#[macro_use]
extern crate clap;
use clap::App;
use std::io::{BufReader, BufRead, Error,  Write, stderr};
use std::fs::File;


mod barcode;

fn main() -> Result<(), Error> {

    // Load parameters
    let yaml = load_yaml!("cli.yml");
    let matches = App::from_yaml(yaml).get_matches();
    
    let path_seq_pool = matches.value_of("seq_pool").unwrap();
    let cell_cycle: u32 = matches.value_of("cell_cycle").unwrap_or("10").parse().unwrap_or_else(|_| panic!("Illegal cell-cycle input"));
    let pcr_cycle: u32 = matches.value_of("pcr_cycle").unwrap_or("15").parse().unwrap_or_else(|_| panic!("Illegal pcr-cycle input"));
    let pcr_mutation_rate: f32 = matches.value_of("pcr_mutation_rate").unwrap_or("0.00001").parse().unwrap_or_else(|_| panic!("Illegal pcr-mutation-rate input"));
    let tissue_param: Vec<u32> = matches.value_of("tissue_param").unwrap_or("10:10:10").to_owned()
        .split(':').map(|s| s.parse().unwrap_or_else(|_| panic!("Illegal tissue cell number input"))).collect();
    let tissue_param = barcode::InitParam { bipotent_n: tissue_param[0], bas_n: tissue_param[1], lum_n: tissue_param[2] };
    let path_indluced_barcode = matches.value_of("path_indluced_barcode").unwrap();
    let path_biological_sample_barcode = matches.value_of("path_biological_sample_barcode").unwrap();
    let path_seq_result = matches.value_of("path_seq_result").unwrap();

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
    let mut tissue = barcode::init_tissue(&tissue_param, &seq_pool);
    barcode::write_tissue_barcode(&tissue, path_indluced_barcode)?;

    // Grow tissue
    for _ in 0..cell_cycle {
        tissue = barcode::grow_tissue(tissue);
    }
    barcode::write_tissue_barcode(&tissue, path_biological_sample_barcode)?;
    let mut dna_template = barcode::extract_dna(tissue);
    
    // PCR
    for i in 0..pcr_cycle {
        // PCR mutation rate: https://www.hindawi.com/journals/mbi/2014/287430/tab1/
        barcode::pcr(&mut dna_template, pcr_mutation_rate);
        stderr().write_fmt(format_args!("PCR cycle {} ...\n", i))?;
    }
    barcode::write_pcr_barcode(&dna_template, path_seq_result)?;

    stderr().write_fmt(format_args!("{} PCR product.\n", dna_template.len()))?;
    // debug code: to print template result, START
    // for i in 0..10 {
        // let key = dna_template.keys().nth(i).unwrap();
        // stderr().write_fmt(format_args!("  {:?}\n", dna_template.get_key_value(key)))?;
    // }
    // END
    Ok(())
}
