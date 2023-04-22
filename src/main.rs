use std::io;
pub mod bio_related_structs;
use bio_related_structs::FastaFile;
use bio_related_structs::Sequence;
use bio_related_structs::SequenceType;
fn main() {
    loop {
        println!(
            "Welcome to rusalind, a Rust implementation of the Rosalind bioinformatics problems."
        );
        println!("Please choose from the menu below the function you would like to run.");
        println!("----------------------------------------");
        println!("1. Counting DNA Nucleotides");
        println!("----------------------------------------");
        println!("2. Counting RNA Nucleotides");
        println!("----------------------------------------");
        println!("3. Transcribing DNA into RNA");
        println!("----------------------------------------");
        println!("4. Reverse Complementing a Strand of DNA");
        println!("----------------------------------------");
        println!("5. Computing GC Content");
        println!("----------------------------------------");
        println!("6. Get the highest GC percentage sequence\nfrom a fasta file");
        println!("----------------------------------------");
        println!("7. Find the positions of a motif in a \nsequence");
        println!("----------------------------------------");
        println!("q. Exit");
        println!("----------------------------------------");
        let mut choice = String::new();
        io::stdin()
            .read_line(&mut choice)
            .expect("Failed to read line");

        match choice.trim() {
            "1" => {
                println!("Please enter the DNA sequence you would like to count");
                let mut dna_sequence = String::new();
                io::stdin()
                    .read_line(&mut dna_sequence)
                    .expect("Failed to read line");
                let dna_sequence = dna_sequence.trim();
                let dna_sequence =
                    Sequence::sequence_from_string(dna_sequence.to_string(), SequenceType::DNA);
                let count_result: (usize, usize, usize, usize) = dna_sequence.get_count_bases();
                println!(
                    "A:{},C:{},G:{},T:{}",
                    count_result.0, count_result.1, count_result.2, count_result.3
                );
            }
            "2" => {
                println!("Please enter the RNA sequence you would like to count");
                let mut rna_sequence = String::new();
                io::stdin()
                    .read_line(&mut rna_sequence)
                    .expect("Failed to read line");
                let rna_sequence = rna_sequence.trim();
                let rna_sequence =
                    Sequence::sequence_from_string(rna_sequence.to_string(), SequenceType::RNA);
                let count_result: (usize, usize, usize, usize) = rna_sequence.get_count_bases();
                println!(
                    "A:{},C:{},G:{},U:{}",
                    count_result.0, count_result.1, count_result.2, count_result.3
                );
            }
            "3" => {
                println!("Please enter the DNA sequence you would like to transcribe");
                let mut dna_sequence = String::new();
                io::stdin()
                    .read_line(&mut dna_sequence)
                    .expect("Failed to read line");
                let dna_sequence = dna_sequence.trim();
                let dna_sequence =
                    Sequence::sequence_from_string(dna_sequence.to_string(), SequenceType::DNA);
                let rna_sequence = dna_sequence.transcribe();
                println!("The RNA sequence is: {}", rna_sequence);
            }
            "4" => {
                println!("Please enter the DNA sequence you would like to complement");
                let mut dna_sequence = String::new();
                io::stdin()
                    .read_line(&mut dna_sequence)
                    .expect("Failed to read line");
                let dna_sequence = dna_sequence.trim();
                let dna_sequence =
                    Sequence::sequence_from_string(dna_sequence.to_string(), SequenceType::DNA);
                let complement_sequence = dna_sequence.reverse_complement();
                println!(
                    "The reverse complement sequence is: {}",
                    complement_sequence
                );
            }
            "5" => {
                println!("Please enter the DNA sequence you would like to compute GC content");
                let mut dna_sequence = String::new();
                io::stdin()
                    .read_line(&mut dna_sequence)
                    .expect("Failed to read line");
                let dna_sequence = dna_sequence.trim();
                let dna_sequence =
                    Sequence::sequence_from_string(dna_sequence.to_string(), SequenceType::DNA);
                let gc_content = dna_sequence.get_gc_percentage();
                println!("The GC content is: {}", gc_content);
            }
            "6" => {
                println!("Please enter the path to the fasta file you would like to get the highest GC percentage sequence from");
                let mut fasta_file_path = String::new();
                io::stdin()
                    .read_line(&mut fasta_file_path)
                    .expect("Failed to read line");
                let fasta_file_path = fasta_file_path.trim().to_string();
                let new_fasta_file: FastaFile = FastaFile::new(fasta_file_path.to_string().into());
                let (highest_gc_content_header, highest_gc_percentage) =
                    new_fasta_file.get_highest_gc_percentage();
                println!(
                    "The highest GC content header is: {}\nWith a GC percentage of: {:.6}",
                    highest_gc_content_header,
                    highest_gc_percentage * 100.0
                );
            }
            "7" => {
                println!("Please enter the sequence you would like to find the motif in");
                let mut sequence = String::new();
                io::stdin()
                    .read_line(&mut sequence)
                    .expect("Failed to read line");
                let sequence = sequence.trim();
                let sequence =
                    Sequence::sequence_from_string(sequence.to_string(), SequenceType::DNA);
                println!("Please enter the motif you would like to find");
                let mut motif = String::new();
                io::stdin()
                    .read_line(&mut motif)
                    .expect("Failed to read line");
                let motif = motif.trim();
                let motif = Sequence::sequence_from_string(motif.to_string(), SequenceType::DNA);
                let motif_positions = sequence.get_k_mer_starts(&motif);
                println!("The motif positions are: {:?}", motif_positions);
            }
            "q" => break,
            _ => println!("Please enter a valid option"),
        }
    }
}
