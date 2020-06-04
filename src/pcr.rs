use rand::seq::SliceRandom;
use rand::rngs::ThreadRng;
use rand::distributions::{Binomial, Normal, Bernoulli, Distribution};
use std::fmt;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufRead, Error,  Write, stderr};

// use rand::distributions::{Bernoulli, Distribution};

pub struct InitParam {
    pub bipotent_n: u32,
    pub lum_n: u32,
    pub bas_n: u32,
}

#[derive(Debug, PartialEq, Clone)]
pub enum CellType {
    BiPotent,
    UniPotentLum,
    UniPotentBas,
}

impl fmt::Display for CellType {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let output =  match self {
            CellType::BiPotent => "BiPotent",
            CellType::UniPotentBas => "Basal",
            CellType::UniPotentLum => "Luminal"
        };
        write!(f, "{}", output)
    }
}

#[derive(Debug, PartialEq, Clone)]
pub struct Cell <'a> {
    pub cell_type: CellType,
    pub barcode: &'a String,
}

impl <'a> Cell <'a> {
    fn new(cell_type: CellType, barcode: &'a String) -> Self {
        Cell {
            cell_type,
            barcode,
        }
    }

    fn cycle_cell (self) -> Tissue <'a> {
        match self.cell_type {
            CellType::BiPotent => vec![Cell {cell_type: CellType::UniPotentLum, barcode: self.barcode}, Cell {cell_type: CellType::UniPotentBas, barcode: self.barcode}],
            _ => vec![self.clone(), self],
        }
    }
}


pub type Tissue <'a> = Vec<Cell <'a>>;


pub type BarcodePool = Vec<String>;

trait BarcodePoolMethods {
    fn pick_seq(&self, rng: &mut ThreadRng) -> &String;
}

impl BarcodePoolMethods for BarcodePool {
    fn pick_seq(&self, rng: &mut ThreadRng) -> &String {
        self.choose(rng).unwrap()
    }
}

pub static NT: [u8; 4] = ['A' as u8, 'T' as u8, 'C' as u8, 'G' as u8];

fn gen_pcr_mutation(seq: &String, rng: &mut ThreadRng) -> String {
    let mut pcr_product = seq.clone().into_bytes();
    let nt_mutation: &mut u8 = pcr_product.choose_mut(rng).unwrap();
    let mutation_to: &u8 = NT.choose(rng).unwrap();
    *nt_mutation = *mutation_to;
    String::from_utf8(pcr_product).unwrap()
}

// fn cal_pcr_mut(distribution: Bernoulli, rng: &mut ThreadRng) -> bool {
// distribution.sample(rng)
// }

pub fn pcr(dna_template: &mut DNATemp, mutation_rate: f64, pcr_efficiency_mean: f64, pcr_efficiency_sd: f64) {
    if pcr_efficiency_sd == 0f64 {
        pcr_model1(dna_template, mutation_rate, pcr_efficiency_mean);
    } else {
        pcr_model6(dna_template, mutation_rate, pcr_efficiency_mean, pcr_efficiency_sd);
    }
}

pub fn pcr_model1(dna_template: &mut DNATemp, mutation_rate: f64, pcr_efficiency: f64) {
    // let mut new_seqs: Vec<String> = Vec::with_capacity(1000_000_000);
    let mut new_seqs: Vec<String> = vec![];
    let mut rng1 = rand::thread_rng();
    let mut rng2 = rand::thread_rng();
    let mut rng3 = rand::thread_rng();

    // calculate pcr mutation error number by Binomail distribution.
    let mut bin = |n| { Binomial::new(n, mutation_rate).sample(&mut rng1) };

    // pcr success?
    let mut bern = || { Bernoulli::new(pcr_efficiency).unwrap().sample(&mut rng2) };

    for (k, v) in dna_template.0.iter_mut() {
        let mutation_count: u64 = bin(v.0);
        for _ in 0..mutation_count as u64 {
            new_seqs.push(gen_pcr_mutation(k, &mut rng3));
        }

        if bern() {
            *v = (v.0 * 2u64 - mutation_count, v.1);
        }
    }

    for seq in new_seqs {
        if let Some(v) = dna_template.0.get_mut(&seq) {
            v.0 += 1;
        } else {
            dna_template.0.insert(seq, (1, pcr_efficiency));
        }
    }
}

// fn bern(pcr_efficiency: f64, rng: &mut ThreadRng) -> bool {
    // Bernoulli::new(pcr_efficiency).unwrap().sample(rng)
// }

pub fn pcr_model6(dna_template: &mut DNATemp, mutation_rate: f64, pcr_efficiency_mean: f64, pcr_efficiency_sd: f64) {
    let mut new_seqs: Vec<String> = vec![];
    let mut rng1 = rand::thread_rng();
    let mut rng2 = rand::thread_rng();
    let mut rng3 = rand::thread_rng();
    let mut rng4 = rand::thread_rng();

    // calculate pcr mutation error number by Binomail distribution.
    let mut bin = |n, length| { Binomial::new(n, mutation_rate * length).sample(&mut rng1) };

    // pcr success or not
    let mut bern = |pcr_efficiency: f64| { Bernoulli::new(pcr_efficiency).unwrap().sample(&mut rng2) };

    // pcr efficiency
    let pcr_efficiency_dist = Normal::new(pcr_efficiency_mean, pcr_efficiency_sd);
    let mut get_pcr_efficiency = || { pcr_efficiency_dist.sample(&mut rng4).abs().min(1.0_f64) };


    for (k, v) in dna_template.0.iter_mut() {
        let mutation_count: u64 = bin(v.0, k.len() as f64);
        for _ in 0..mutation_count as u64 {
            new_seqs.push(gen_pcr_mutation(k, &mut rng3));
        }

        if bern(v.1) {
            *v = (v.0 * 2u64 - mutation_count, v.1);
        }
    }

    for seq in new_seqs {
        if let Some(v) = dna_template.0.get_mut(&seq) {
            v.0 += 1;
        } else {
            let pcr_efficiency = dna_template.0.insert(seq, (1, get_pcr_efficiency()));
        }
    }
}

pub struct DNATemp (HashMap<String, (u64, f64)>);

impl DNATemp {
    fn new() -> Self {
        // DNATemp(HashMap::with_capacity(100_000_000))
        DNATemp(HashMap::new())
    }
}

pub fn read_template(file: String, pcr_efficiency_mean: f64, pcr_efficiency_sd: f64) -> Result<DNATemp, Error> {

    let mut rng4 = rand::thread_rng();
    let pcr_efficiency_dist = Normal::new(pcr_efficiency_mean, pcr_efficiency_sd);
    let mut get_pcr_efficiency = || { pcr_efficiency_dist.sample(&mut rng4).abs().min(1.0_f64) };

    let mut dna_template: DNATemp = DNATemp::new(); 

    let f = File::open(file)?;
    let f = BufReader::new(f);
    let mut lines = f.lines();
    lines.next();
    for l in lines {
        // let l: String = l.unwrap().clone();
        // let seq = l.split(';').last().unwrap();
        let seq = l.unwrap();

        if let Some(v) = dna_template.0.get_mut(&seq) {
            *v = (v.0 + 1, v.1);
        } else {
            dna_template.0.insert(seq.to_owned(), (1, get_pcr_efficiency()));
        }
    }

    Ok(dna_template)
}

pub fn init_tissue <'a> (init_par: &InitParam, barcodepool: &'a BarcodePool) -> Tissue <'a> {
    let mut tissue: Tissue = Tissue::new();
    let mut rng = rand::thread_rng();
    for _ in 0..(init_par.bipotent_n) {
        tissue.push(Cell::new(CellType::BiPotent, barcodepool.pick_seq(&mut rng)));
    }
    for _ in 0..(init_par.lum_n) {
        tissue.push(Cell::new(CellType::UniPotentLum, barcodepool.pick_seq(&mut rng)));
    }
    for _ in 0..(init_par.bas_n) {
        tissue.push(Cell::new(CellType::UniPotentBas, barcodepool.pick_seq(&mut rng)));
    }
    tissue
}

pub fn grow_tissue(tissue: Tissue) -> Tissue {
    let mut new_tissue: Tissue = Tissue::new();
    for cell in tissue {
        new_tissue.append(&mut cell.cycle_cell());
    } 
    new_tissue
}

pub fn write_pcr_barcode(pcr_template: DNATemp, output_path: &str) -> Result<(), Error> {
    let f = File::create(output_path);
    if let Ok(mut f) = f {
        f.write_all(b"barcode\tcount\n")?;
        for (k, v) in pcr_template.0 {
            f.write_fmt(format_args!("{}\t{}\n", k, v.0))?;
        }
    } else {
        let error_msg = format!("Can not find {}.", &output_path);
        panic!(error_msg);
    }
    Ok(())
}

pub fn write_tissue_barcode(tissue: &Tissue, output_path: &str) -> Result<(), Error> {
    let f = File::create(output_path);
    if let Ok(mut f) = f {
        f.write_all(b"barcode\tcell_type\n")?;
        for cell in tissue {
            f.write_fmt(format_args!("{}\t{}\n", cell.barcode, cell.cell_type))?;
        }
    } else {
        let error_msg = format!("Can not find {}.", &output_path);
        panic!(error_msg);
    }
    Ok(())
}

#[test]
    fn cell_test() {
        let seq_pool: BarcodePool = vec!["ATCA".to_owned()];
        let tissue_param = InitParam { bipotent_n : 1, bas_n : 0, lum_n : 0 };
        assert_eq!(init_tissue(&tissue_param, &seq_pool), vec![Cell { cell_type: CellType::BiPotent, barcode: &String::from("ATCA") }]);
    }

// #[test]
// fn extract_dna_test() {
//     let seq_pool: BarcodePool = vec!["ATCA".to_owned(), "TTCAG".to_owned()];
//     let tissue_param = InitParam { bipotent_n : 10, bas_n : 10, lum_n : 10 };
//     let mut tissue = init_tissue(&tissue_param, &seq_pool);
//     for _ in 0..5 {
//         tissue = grow_tissue(tissue);
//     }
//     let mut dna_template = extract_dna(&tissue, CellType::UniPotentBas);
//     println!("{:?} DNATemp generated.", dna_template);
//     for _ in 0..10 {
//         pcr(&mut dna_template, 0.0001);
//     }
//     println!("{:?} PCR product.", dna_template);
// }
