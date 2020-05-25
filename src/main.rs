#[macro_use]
extern crate clap;
use clap::App;
use std::io::{BufReader, BufRead, Error,  Write, stderr};
use std::fs::File;

mod pcr;


fn main() -> Result<(), Error> {

    // Load parameters
    let yaml = load_yaml!("cli.yml");
    let matches = App::from_yaml(yaml).get_matches();

    let template: String = matches.value_of("template").unwrap().to_string();
    let pcr_cycle: u32 = matches.value_of("pcr_cycle").unwrap_or("15").parse().unwrap_or_else(|_| panic!("Illegal pcr-cycle input"));
    let pcr_mutation_rate: f64 = matches.value_of("pcr_mutation_rate").unwrap_or("0.00001").parse().unwrap_or_else(|_| panic!("Illegal pcr-mutation-rate input"));
    let pcr_efficiency_mean: f64 = matches.value_of("pcr_efficiency_mean").unwrap_or("0.5").parse().unwrap_or_else(|_| panic!("Illegal pcr-efficiency input"));
    let pcr_efficiency_sd: f64 = matches.value_of("pcr_efficiency_sd").unwrap_or("1").parse().unwrap_or_else(|_| panic!("Illegal pcr-efficiency sd input"));
    let prefix_output = matches.value_of("prefix_output").unwrap().to_string();
    let path_seq = prefix_output.clone() + "_seq_barcode.txt";

    let mut dna_template = pcr::read_template(template, pcr_efficiency_mean, pcr_efficiency_sd)?;
    // PCR

    for i in 0..pcr_cycle {
        pcr::pcr(&mut dna_template, pcr_mutation_rate, pcr_efficiency_mean, pcr_efficiency_sd);
        stderr().write_fmt(format_args!("PCR cycle {} ...\n", i))?;
    }
    pcr::write_pcr_barcode(dna_template, &path_seq)?;
    // stderr().write_fmt(format_args!("{} PCR product.\n", dna_template.len()))?;
    // stderr().write_fmt(format_args!("{} PCR product.\n", dna_template.len()))?;
    // debug code: to print template result, START
    // for i in 0..10 {
        // let key = dna_template.keys().nth(i).unwrap();
        // stderr().write_fmt(format_args!("  {:?}\n", dna_template.get_key_value(key)))?;
    // }
    // END
    Ok(())
}
