use std::{
    hint::black_box,
    os::unix::fs::MetadataExt,
    path::{Path, PathBuf},
    sync::{
        atomic::{AtomicUsize, Ordering},
        Arc, Mutex,
    },
};

use clap::Parser;
use crossbeam::channel::{Receiver, Sender};
use memmap2::Mmap;
use paraseq::{
    fasta,
    prelude::{ParallelProcessor, ParallelReader},
    ProcessError, Record,
};

#[derive(clap::Parser)]
struct Args {
    #[clap(long, short, default_value_t = 0)]
    threads: usize,

    #[clap(default_value = "human-genome.fa")]
    path: PathBuf,

    #[clap(long)]
    needletail: bool,
}

fn main() {
    let mut args = Args::parse();
    args.needletail = true;

    let size = std::fs::File::open(&args.path)
        .unwrap()
        .metadata()
        .unwrap()
        .size();

    {
        let start = std::time::Instant::now();
        needletail(&args.path);
        let elapsed = start.elapsed();
        eprintln!(
            "Needletail {:5.2}s {:>5.2}GB/s",
            elapsed.as_secs_f32(),
            size as f32 / elapsed.as_nanos() as f32
        );
    }

    for t in [6, 5, 4, 3, 2, 1] {
        eprintln!("Threads: {t}");
        args.threads = t;

        let start = std::time::Instant::now();
        mmap(&args);
        let elapsed = start.elapsed();
        eprintln!(
            "Mmap       {:>5.2}s {:>5.2}GB/s",
            elapsed.as_secs_f32(),
            size as f32 / elapsed.as_nanos() as f32
        );

        let start = std::time::Instant::now();
        paraseq(&args.path, args.threads);
        let elapsed = start.elapsed();
        eprintln!(
            "Paraseq    {:>5.2}s {:>5.2}GB/s",
            elapsed.as_secs_f32(),
            size as f32 / elapsed.as_nanos() as f32
        );
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

    let count = AtomicUsize::new(0);

    std::thread::scope(|scope| {
        let (sender, receiver) = crossbeam::channel::bounded(args.threads);
        let len = data.len().div_ceil(args.threads);

        // Spawn producer threads.
        for i in 0..args.threads {
            let start = (i * len).min(data.len());
            let end = ((i + 1) * len).min(data.len());
            let sender = sender.clone();
            scope.spawn(move || producer(data, start, end, sender));
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

fn producer<'d>(data: &'d [u8], mut start: usize, slice_end: usize, sender: Sender<Vec<&'d [u8]>>) {
    if let Some(pos) = memchr::memchr(b'>', &data[start..]) {
        start += pos;
    } else {
        return;
    }

    let mut batch = vec![];
    let mut batch_len = 0;
    let target_batch_len = 1 << 20; // 1MB

    loop {
        let end = match memchr::memchr(b'>', &data[start + 1..]) {
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
    let processor = SeqSum::default();

    // let batch_size = DEFAULT_MAX_RECORDS;
    let batch_size = 1;
    let reader = fasta::Reader::with_batch_size(file, batch_size).unwrap();
    reader.process_parallel(processor.clone(), threads).unwrap();

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
impl ParallelProcessor for SeqSum {
    fn process_record<Rf: Record>(&mut self, record: Rf) -> Result<(), ProcessError> {
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
