use std::{
    hint::black_box,
    os::unix::fs::MetadataExt,
    path::{Path, PathBuf},
    sync::{
        atomic::{AtomicUsize, Ordering},
        Arc, Mutex,
    },
    time::Duration,
};

use clap::Parser;
use crossbeam::channel::{Receiver, Sender};
use helicase::{
    config::advanced::DEFAULT_CONFIG,
    input::{FromFile, FromSlice},
    Config, ParserOptions,
};
use memmap2::Mmap;
use paraseq::{
    prelude::{ParallelProcessor, ParallelReader},
    ProcessError, Record, DEFAULT_MAX_RECORDS,
};

#[derive(clap::Parser)]
struct Args {
    #[clap(long, short, default_value_t = 0)]
    threads: usize,

    #[clap(default_value = "human-genome.fa")]
    path: PathBuf,
}

fn main() {
    let mut args = Args::parse();

    let size = std::fs::File::open(&args.path)
        .unwrap()
        .metadata()
        .unwrap()
        .size();

    println!("threads\tus\tneedletail\tparaseq");

    let print = |name: &str, e: Duration| {
        let secs = e.as_secs_f32();
        let thpt = size as f32 / e.as_nanos() as f32;
        eprintln!("{name:<20}  {secs:>5.2}s {thpt:>5.2}GB/s",);
        print!("\t{thpt:>5.2}");
    };

    for t in [6, 3, 1] {
        print!("{t}");
        eprintln!("Threads: {t}");
        args.threads = t;

        let start = std::time::Instant::now();
        mmap(&args);
        print("mmap trivial", start.elapsed());

        let start = std::time::Instant::now();
        helicase_reader::<DEFAULT_CONFIG>(&args);
        print("helicase reader dft", start.elapsed());

        let start = std::time::Instant::now();
        helicase_mmap::<DEFAULT_CONFIG>(&args);
        print("mmap helicase dft", start.elapsed());

        if t == 1 {
            let start = std::time::Instant::now();
            needletail(&args.path);
            print("Needletail", start.elapsed());

            let start = std::time::Instant::now();
            paraseq(&args.path, args.threads);
            print("Paraseq", start.elapsed());

            let start = std::time::Instant::now();
            helicase_default(&args.path);
            print("helicase dft", start.elapsed());

            let start = std::time::Instant::now();
            helicase_ignore_all(&args.path);
            print("helicase ignore", start.elapsed());

            let start = std::time::Instant::now();
            helicase_header_only(&args.path);
            print("helicase hdr", start.elapsed());

            let start = std::time::Instant::now();
            helicase_string(&args.path);
            print("helicase str", start.elapsed());

            let start = std::time::Instant::now();
            helicase_columnar(&args.path);
            print("helicase col", start.elapsed());

            let start = std::time::Instant::now();
            helicase_packed(&args.path);
            print("helicase pck", start.elapsed());
        }
        println!();
    }
}

fn mmap(args: &Args) {
    let mmap: Mmap;
    let data;
    {
        let file = std::fs::File::open(&args.path).unwrap();
        mmap = unsafe { Mmap::map(&file).unwrap() };
        data = &*mmap;
    }

    let ext = args.path.extension().unwrap();
    let fasta = ext == "fa" || ext == "fasta" || ext == "fna";
    let symbol = if fasta { b'>' } else { b'@' };

    let count = AtomicUsize::new(0);

    std::thread::scope(|scope| {
        let (sender, receiver) = crossbeam::channel::bounded(args.threads);
        let len = data.len().div_ceil(args.threads);

        // Spawn producer threads.
        for i in 0..args.threads {
            let start = (i * len).min(data.len());
            let end = ((i + 1) * len).min(data.len());
            let sender = sender.clone();
            scope.spawn(move || producer(symbol, data, start, end, sender));
        }

        // Spawn consumer threads.
        for _ in 0..args.threads {
            let receiver = receiver.clone();
            let count = &count;
            scope.spawn(|| count.fetch_add(consumer(receiver), Ordering::Relaxed));
        }
    });

    // assert_eq!(count.into_inner(), 2397);
}

fn producer<'d>(
    symbol: u8,
    data: &'d [u8],
    mut start: usize,
    slice_end: usize,
    sender: Sender<Vec<&'d [u8]>>,
) {
    if let Some(pos) = memchr::memchr(symbol, &data[start..]) {
        start += pos;
    } else {
        return;
    }

    let mut batch = vec![];
    let mut batch_len = 0;
    let target_batch_len = 1 << 20; // 1MB

    loop {
        let end = match memchr::memchr(symbol, &data[start + 1..]) {
            Some(x) => start + 1 + x,
            None => slice_end,
        };
        batch.push(&data[start..end]);
        batch_len += end - start;
        if batch_len >= target_batch_len {
            sender.send(batch).unwrap();
            batch = vec![];
            batch_len = 0;
        }

        if end >= slice_end {
            break;
        }
        start = end;
    }

    if batch.len() > 0 {
        sender.send(batch).unwrap();
    }
}

fn consumer(receiver: Receiver<Vec<&[u8]>>) -> usize {
    let mut count = 0;
    for batch in receiver.iter() {
        for record in batch {
            count += 1;
            let first_newline = memchr::memchr(b'\n', record).unwrap();
            black_box(first_newline);
        }
    }
    count
}

fn helicase_reader<const CONFIG: Config>(args: &Args) {
    let reader =
        helicase::paraseq_reader::ParallelHelicaseReader::<'_, CONFIG>::new(&args.path, 16);

    #[derive(Clone)]
    struct Processor<'a> {
        global_count: &'a AtomicUsize,
        count: usize,
    }
    impl<'a, 'b, const CONFIG: Config>
        paraseq::parallel::ParallelProcessor<&'b helicase::FastxParser<'b, CONFIG>>
        for Processor<'a>
    {
        fn process_record(
            &mut self,
            _record: &helicase::FastxParser<'_, CONFIG>,
        ) -> Result<(), ProcessError> {
            self.count += 1;
            Ok(())
        }
        fn on_thread_complete(&mut self) -> Result<(), ProcessError> {
            self.global_count.fetch_add(self.count, Ordering::Relaxed);
            Ok(())
        }
    }

    let global_count = AtomicUsize::new(0);
    let mut processor = Processor {
        global_count: &global_count,
        count: 0,
    };

    reader
        .process_parallel(&mut processor, args.threads)
        .unwrap();
}

fn helicase_mmap<const CONFIG: Config>(args: &Args) {
    let mmap: Mmap;
    let data;
    {
        let file = std::fs::File::open(&args.path).unwrap();
        mmap = unsafe { Mmap::map(&file).unwrap() };
        data = &*mmap;
    }

    let ext = args.path.extension().unwrap();
    let fasta = ext == "fa" || ext == "fasta" || ext == "fna";
    let start = if fasta { b'>' } else { b'@' };

    let count = AtomicUsize::new(0);

    std::thread::scope(|scope| {
        let len = data.len().div_ceil(args.threads);

        let mut splits = (0..=args.threads)
            .map(|i| (i * len).min(data.len()))
            .collect::<Vec<usize>>();
        for x in &mut splits {
            if *x == 0 {
                continue;
            }
            if *x >= data.len() {
                continue;
            }
            // find first > or @ preceded by \n after x
            loop {
                if let Some(pos) = memchr::memchr(start, &data[*x..]) {
                    if data[*x + pos - 1] == b'\n' {
                        *x += pos;
                        break;
                    }
                    *x += pos + 1;
                } else {
                    *x = data.len();
                    break;
                }
            }
        }

        for i in 0..args.threads {
            let start = splits[i];
            let end = splits[i + 1];
            let count = &count;
            scope.spawn(move || {
                let x = helicase::FastxParser::<CONFIG>::from_slice(&data[start..end]);
                let local_count = x.count();
                count.fetch_add(local_count, Ordering::Relaxed);
            });
        }
    });
}

fn helicase_default(path: &Path) -> usize {
    let x = helicase::FastxParser::<DEFAULT_CONFIG>::from_file(path).unwrap();
    x.count()
}

fn helicase_header_only(path: &Path) -> usize {
    const CONFIG: Config = ParserOptions::default().ignore_dna().config();
    let x = helicase::FastxParser::<CONFIG>::from_file(path).unwrap();
    x.count()
}

fn helicase_ignore_all(path: &Path) -> usize {
    const CONFIG: Config = ParserOptions::default()
        .ignore_headers()
        .ignore_dna()
        .config();
    let x = helicase::FastxParser::<CONFIG>::from_file(path).unwrap();
    x.count()
}

fn helicase_string(path: &Path) -> usize {
    const CONFIG: Config = ParserOptions::default()
        .ignore_headers()
        .dna_string()
        .config();
    let x = helicase::FastxParser::<CONFIG>::from_file(path).unwrap();
    x.count()
}

fn helicase_columnar(path: &Path) -> usize {
    const CONFIG: Config = ParserOptions::default()
        .ignore_headers()
        .dna_columnar()
        .config();
    let x = helicase::FastxParser::<CONFIG>::from_file(path).unwrap();
    x.count()
}

fn helicase_packed(path: &Path) -> usize {
    const CONFIG: Config = ParserOptions::default()
        .ignore_headers()
        .dna_packed()
        .config();
    let x = helicase::FastxParser::<CONFIG>::from_file(path).unwrap();
    x.count()
}

fn needletail(path: &Path) -> usize {
    let mut count = 0;
    let mut reader = needletail::parse_fastx_file(path).unwrap();
    while let Some(_record) = reader.next() {
        count += 1;
    }
    // assert_eq!(count, 2397);
    count
}

fn paraseq(path: &Path, threads: usize) -> usize {
    let file = std::fs::File::open(&path).unwrap();
    let mut processor = SeqSum::default();

    let batch_size = DEFAULT_MAX_RECORDS;
    // let batch_size = 1;
    let reader = paraseq::fastx::Reader::new_with_batch_size(file, batch_size).unwrap();
    reader.process_parallel(&mut processor, threads).unwrap();

    // assert_eq!(processor.get_num_records(), 2397);
    processor.get_num_records() as usize
}

// see: https://github.com/noamteyssier/paraseq/blob/main/examples/parallel.rs

#[derive(Default, Clone)]
pub struct SeqSum {
    /// Thread local number of records
    pub num_records: u64,
    /// Global number of records
    pub global_num_records: Arc<Mutex<u64>>,
}
impl SeqSum {
    #[must_use]
    pub fn get_num_records(&self) -> u64 {
        *self.global_num_records.lock().unwrap()
    }
}
impl<Rf: Record> ParallelProcessor<Rf> for SeqSum {
    fn process_record(&mut self, record: Rf) -> Result<(), ProcessError> {
        black_box(record);
        self.num_records += 1;
        Ok(())
    }
    fn on_batch_complete(&mut self) -> Result<(), ProcessError> {
        *self.global_num_records.lock().unwrap() += self.num_records;
        self.num_records = 0;
        Ok(())
    }
}
