use rayon::prelude::{
	IndexedParallelIterator, IntoParallelIterator, IntoParallelRefIterator, ParallelIterator,
};
use std::fmt;
use std::sync::{Arc, Mutex};
use std::{io::BufRead, path::PathBuf, sync::atomic::AtomicUsize};
#[derive(Debug)]
pub struct FastaFile {
	pub file_path: PathBuf,
	pub fasta_entries: Vec<FastaEntry>,
}
impl FastaFile {
	pub fn new(file_path: PathBuf) -> FastaFile {
		let mut new_fasta_file = FastaFile {
			file_path,
			fasta_entries: Vec::new(),
		};
		let ext_file = std::fs::File::open(&new_fasta_file.file_path).unwrap();
		let reader = std::io::BufReader::new(ext_file);
		for line in reader.lines().map(|l| l.unwrap()) {
			if line.starts_with('>') {
				new_fasta_file
					.fasta_entries
					.push(FastaEntry::new(line, Sequence::new(SequenceType::DNA)))
			} else {
				line.chars().for_each(|base| {
					new_fasta_file
						.fasta_entries
						.last_mut()
						.unwrap()
						.sequence
						.add_nucleotide(Nucleotide::from_char(base).expect("Invalid base found"))
				});
			}
		}
		new_fasta_file
	}
	pub fn get_highest_gc_percentage(&self) -> (String, f64) {
		let highest_gc_percentages = Arc::new(Mutex::new(vec![]));
		self.fasta_entries.par_iter().for_each(|entry| {
			let mut local_highest_gc_percentge: f64 = 0.0;
			let mut local_highest_gc_header: String = String::new();
			let gc_percentage = entry.sequence.get_gc_percentage();
			if gc_percentage > local_highest_gc_percentge {
				local_highest_gc_percentge = gc_percentage;
				local_highest_gc_header = entry.header.clone();
			}
			highest_gc_percentages
				.lock()
				.unwrap()
				.push((local_highest_gc_header, local_highest_gc_percentge));
		});
		let highest_gc_percentages = highest_gc_percentages.lock().unwrap();
		let (global_max_heder, global_max_gc_percentage) = highest_gc_percentages
			.iter()
			.max_by(|(_, gc1), (_, gc2)| gc1.partial_cmp(gc2).unwrap())
			.unwrap();
		(global_max_heder.clone(), global_max_gc_percentage.clone())
	}
}
#[derive(Debug)]
pub struct FastaEntry {
	pub header: String,
	pub sequence: Sequence,
}
impl FastaEntry {
	pub fn new(header: String, sequence: Sequence) -> FastaEntry {
		FastaEntry { header, sequence }
	}
}
#[derive(Debug)]
pub struct Sequence {
	pub sequence_type: SequenceType,
	pub sequence: Vec<Nucleotide>,
}
impl Sequence {
	fn new(sequence_type: SequenceType) -> Sequence {
		Sequence {
			sequence_type,
			sequence: Vec::new(),
		}
	}
	pub fn sequence_from_file(file_path: PathBuf, sequence_type: SequenceType) -> Sequence {
		let mut new_sequence = Sequence::new(sequence_type);
		let ext_file = std::fs::File::open(file_path).unwrap();
		let reader = std::io::BufReader::new(ext_file);
		for line in reader.lines().map(|l| l.unwrap()) {
			if !line.starts_with('>') {
				line.chars().for_each(|base| {
					new_sequence
						.add_nucleotide(Nucleotide::from_char(base).expect("Invalid base found"))
				});
			}
		}
		new_sequence
	}
	pub fn sequence_from_string(sequence: String, sequence_type: SequenceType) -> Sequence {
		let mut new_sequence = Sequence::new(sequence_type);
		sequence.chars().for_each(|base| {
			new_sequence.add_nucleotide(Nucleotide::from_char(base).expect("Invalid base found"))
		});
		new_sequence
	}
	fn add_nucleotide(&mut self, nucleotide: Nucleotide) {
		self.sequence.push(nucleotide);
	}
	pub fn get_count_bases(&self) -> (usize, usize, usize, usize) {
		match self.sequence_type {
			SequenceType::DNA => {
				let a_count = AtomicUsize::new(0);
				let c_count = AtomicUsize::new(0);
				let g_count = AtomicUsize::new(0);
				let t_count = AtomicUsize::new(0);
				self.sequence.par_iter().for_each(|base| match base {
					Nucleotide::A => {
						a_count.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
					}
					Nucleotide::C => {
						c_count.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
					}
					Nucleotide::G => {
						g_count.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
					}
					Nucleotide::T => {
						t_count.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
					}
					_ => (),
				});
				(
					a_count.load(std::sync::atomic::Ordering::Relaxed),
					c_count.load(std::sync::atomic::Ordering::Relaxed),
					g_count.load(std::sync::atomic::Ordering::Relaxed),
					t_count.load(std::sync::atomic::Ordering::Relaxed),
				)
			}
			SequenceType::RNA => {
				let a_count = AtomicUsize::new(0);
				let c_count = AtomicUsize::new(0);
				let g_count = AtomicUsize::new(0);
				let u_count = AtomicUsize::new(0);
				self.sequence.par_iter().for_each(|base| match base {
					Nucleotide::A => {
						a_count.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
					}
					Nucleotide::C => {
						c_count.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
					}
					Nucleotide::G => {
						g_count.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
					}
					Nucleotide::U => {
						u_count.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
					}
					_ => (),
				});
				(
					a_count.load(std::sync::atomic::Ordering::Relaxed),
					c_count.load(std::sync::atomic::Ordering::Relaxed),
					g_count.load(std::sync::atomic::Ordering::Relaxed),
					u_count.load(std::sync::atomic::Ordering::Relaxed),
				)
			}
		}
	}
	pub fn transcribe(&self) -> Sequence {
		let new_sequence: Vec<Nucleotide> = match self.sequence_type {
			SequenceType::DNA => self
				.sequence
				.par_iter()
				.map(|base| match base {
					Nucleotide::T => Nucleotide::U,
					_ => base.clone(),
				})
				.collect(),
			SequenceType::RNA => self
				.sequence
				.par_iter()
				.map(|base| match base {
					Nucleotide::U => Nucleotide::T,
					_ => base.clone(),
				})
				.collect(),
		};
		Sequence {
			sequence_type: match self.sequence_type {
				SequenceType::DNA => SequenceType::RNA,
				SequenceType::RNA => SequenceType::DNA,
			},
			sequence: new_sequence,
		}
	}
	pub fn reverse_complement(&self) -> Sequence {
		let new_sequence: Vec<Nucleotide> = match self.sequence_type {
			SequenceType::DNA => self
				.sequence
				.par_iter()
				.map(|base| match base {
					Nucleotide::T => Nucleotide::A,
					Nucleotide::A => Nucleotide::T,
					Nucleotide::C => Nucleotide::G,
					Nucleotide::G => Nucleotide::C,
					_ => base.clone(),
				})
				.rev()
				.collect(),
			SequenceType::RNA => self
				.sequence
				.par_iter()
				.map(|base| match base {
					Nucleotide::U => Nucleotide::A,
					Nucleotide::A => Nucleotide::U,
					Nucleotide::C => Nucleotide::G,
					Nucleotide::G => Nucleotide::C,
					_ => base.clone(),
				})
				.rev()
				.collect(),
		};
		Sequence {
			sequence_type: match self.sequence_type {
				SequenceType::DNA => SequenceType::DNA,
				SequenceType::RNA => SequenceType::RNA,
			},
			sequence: new_sequence,
		}
	}
	pub fn get_gc_percentage(&self) -> f64 {
		let mut gc_count = 0;
		self.sequence.iter().for_each(|base| match base {
			Nucleotide::G | Nucleotide::C => gc_count += 1,
			_ => (),
		});
		gc_count as f64 / self.sequence.len() as f64
	}
	pub fn get_k_mer_starts(&self, k_mer: &Sequence) -> Vec<usize> {
		let n = self.sequence.len();
		let k = k_mer.sequence.len();
		(0..n - k + 1)
			.into_par_iter()
			.filter(|&i| self.sequence[i..i + k] == k_mer.sequence)
			.map(|i| i + 1)
			.collect()
	}
}
impl fmt::Display for Sequence {
	fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
		let mut s = String::new();
		self.sequence.iter().for_each(|base| match base {
			Nucleotide::A => s.push('A'),
			Nucleotide::C => s.push('C'),
			Nucleotide::G => s.push('G'),
			Nucleotide::T => s.push('T'),
			Nucleotide::U => s.push('U'),
		});
		write!(f, "{}", s)
	}
}
#[derive(Eq, PartialEq, Hash, Debug, Clone, Copy)]
pub enum Nucleotide {
	A,
	C,
	G,
	T,
	U,
}
impl Nucleotide {
	fn from_char(c: char) -> Option<Nucleotide> {
		match c.to_ascii_uppercase() {
			'A' => Some(Nucleotide::A),
			'C' => Some(Nucleotide::C),
			'G' => Some(Nucleotide::G),
			'T' => Some(Nucleotide::T),
			'U' => Some(Nucleotide::U),
			_ => None,
		}
	}
}
#[derive(Debug, Clone, Copy)]
pub enum SequenceType {
	DNA,
	RNA,
}
