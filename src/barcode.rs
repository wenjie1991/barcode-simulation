use rand::seq::SliceRandom;
use rand::rngs::ThreadRng;
use std::fmt;
use std::collections::HashMap;
use std::fs::File;
use std::io::{Error, Write};

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

pub fn pcr(dna_template: &mut DNATemp, mutation_rate: f32) {
    let mut new_seqs: Vec<String> = Vec::new();
    // let dist: Bernoulli = Bernoulli::new(0.001).unwrap();
    let mut rng = rand::thread_rng();

    for (k, v) in dna_template.iter_mut() {
        let mutation_count: f32 = *v as f32 * mutation_rate; 
        for _ in 1..mutation_count as u64 {
            new_seqs.push(gen_pcr_mutation(k, &mut rng));
        }
        *v = *v * 2 as u64 - mutation_count as u64;
    }

    for seq in new_seqs {
        if let Some(v) = dna_template.get_mut(&seq) {
            *v += 1
        } else {
            dna_template.insert(seq, 1);
        }
    }
}

pub type DNATemp = HashMap<String, u64>;

pub fn extract_dna(tissue: Tissue) -> DNATemp {
    let mut dna_template: DNATemp = DNATemp::new(); // DNATemp::new();
    for cell in tissue {
        if let Some(v) = dna_template.get_mut(cell.barcode) {
            *v += 1;
        } else {
            dna_template.insert(cell.barcode.to_owned(), 1);
        }
    }
    dna_template
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

pub fn write_pcr_barcode(pcr_template: &DNATemp, output_path: &str) -> Result<(), Error> {
    let f = File::create(output_path);
    if let Ok(mut f) = f {
        f.write_all(b"barcode\tcount\n")?;
        for (k, v) in pcr_template {
            f.write_fmt(format_args!("{}\t{}\n", k, v))?;
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

#[test]
fn extract_dna_test() {
    let seq_pool: BarcodePool = vec!["ATCA".to_owned(), "TTCAG".to_owned()];
    let tissue_param = InitParam { bipotent_n : 10, bas_n : 10, lum_n : 10 };
    let mut tissue = init_tissue(&tissue_param, &seq_pool);
    for _ in 0..5 {
        tissue = grow_tissue(tissue);
    }
    let mut dna_template = extract_dna(tissue);
    println!("{:?} DNATemp generated.", dna_template);
    for _ in 0..10 {
        pcr(&mut dna_template, 0.0001);
    }
    println!("{:?} PCR product.", dna_template);
}
