use std::io::{self, BufReader, BufRead, Error,  Write};
use std::env;
use std::fs::File;

mod barcode;

fn main() -> Result<(), Error> {
    let args: Vec<String> = env::args().collect();
    if args.len() != 5 {
        panic!("Run with incorrect parameters ...");
    } 
    // Input file: barcode pool as input
    let path_seq_pool = &args[1];
    // Output file: induced barcodes
    let path_indluced_barcode = &args[2];
    // Output file: barcodes in biological sample
    let path_biological_sample_barcode = &args[3];
    // Output file: sequencing result
    let path_seq_result = &args[4];

    // Read the barcode pool
    let mut seq_pool: barcode::BarcodePool = vec![];
    let f = File::open(path_seq_pool)?;
    let f = BufReader::new(f);
    let mut lines = f.lines();
    lines.next();
    for l in lines {
        let l: String = l.unwrap().clone();
        let seq = l.split(';').last().unwrap();
        seq_pool.push(String::from(seq));
    }

    // Induction
    let tissue_param = barcode::InitParam { bipotent_n : 10, bas_n : 10, lum_n : 10 };
    let mut tissue = barcode::init_tissue(&tissue_param, &seq_pool);
    barcode::write_tissue_barcode(&tissue, path_indluced_barcode)?;

    // Grow tissue
    for _ in 0..15 {
        tissue = barcode::grow_tissue(tissue);
    }
    barcode::write_tissue_barcode(&tissue, path_biological_sample_barcode)?;
    let mut dna_template = barcode::extract_dna(tissue);
    
    // PCR
    for i in 0..20 {
        // PCR mutation rate: https://www.hindawi.com/journals/mbi/2014/287430/tab1/
        barcode::pcr(&mut dna_template, 0.0001);
        io::stderr().write_fmt(format_args!("PCR cycle {} ...\n", i))?;
    }
    barcode::write_pcr_barcode(&dna_template, path_seq_result)?;

    io::stderr().write_fmt(format_args!("{} PCR product.", dna_template.len()))?;
    Ok(())
}


