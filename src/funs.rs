use rust_htslib::bcf::header::*;
use rust_htslib::bcf::{Read, Reader};
use std::collections::HashMap;
//use std::convert::TryFrom;
//use std::f64::consts::E;
use std::fs::File;
use std::io::{BufRead, BufReader};

//run the whole salami
pub fn make_freqs(vcf_file: &str, pop_key: &str, frequency: bool) {
    let the_pop_key = process_popkey(pop_key);
    let locus_map = freq_map(&the_pop_key);

    let mut locus_keys: Vec<String> = locus_map.keys().cloned().collect();
    locus_keys.sort();

    //print header
    print!("contig position ");
    for key in locus_keys {
        print!("{} ", key,);
    }
    println!();
    let reader = Reader::from_path(&vcf_file).expect("Error opening file.");
    process_vcf(reader, the_pop_key, frequency);
}

//make hashmap according to which individuals belong to which pops
//first column is the 0-index of where individuals are along the vcf header.
//second is individuals that are in the VCF
//thirs  is the population IDS
pub fn process_popkey(file_name: &str) -> HashMap<i32, (String, String)> {
    let f = File::open(file_name).expect("Unable to open file");
    let reader = BufReader::new(f);

    let mut pop_map = HashMap::new();
    for line in reader.lines() {
        let line = line.expect("Unable to read line");
        let pop_ind = line_splitter(&line);
        //convert the str to String, then to an integer
        let pop_id = pop_ind[0].to_string().parse::<i32>().unwrap();
        let ind_id = (pop_ind[1].to_owned(), pop_ind[2].to_owned());
        pop_map.entry(pop_id.to_owned()).or_insert(ind_id);
    }
    pop_map
}

//build new hashmap for pop allele freqs at a locus
pub fn freq_map(
    pop_map: &HashMap<i32, (String, String)>, //the hashmap for the pop and ind. indices
) -> HashMap<String, Vec<f32>> {
    let mut locus_map = HashMap::new();
    for (_key, value) in pop_map {
        locus_map
            .entry(value.1.to_owned()) //index 1 has the pop IDs
            .or_insert(vec![0.0, 0.0]); //start allele counts at 0 for ref and alt
    }
    locus_map
}

//count of ref and alternate alleles for each pop
//(This could use some refactoring)
pub fn process_vcf(mut reader: Reader, pop_map: HashMap<i32, (String, String)>, frequency: bool) {
    let header = reader.header().clone();
    for record_result in reader.records() {
        let record = record_result.expect("Fail to read record");

        let contig_num = header
            .rid2name(record.rid().expect("fail to read reference id"))
            .expect("fail to read ref name");
        let contig = std::str::from_utf8(contig_num).unwrap().to_owned();
        let pos = record.pos();

        let gp_record = record
            .format(b"GP")
            .float()
            .expect("Couldn't retrieve GP field");

        let locus_map = locus_freqs(&pop_map, gp_record);

        //sort the vector of population keys to print
        let mut locus_keys: Vec<String> = locus_map.keys().cloned().collect();
        locus_keys.sort();

        //print output as ref allele frequency OR counts of ref and alt alleles
        print!("{} {} ", contig, pos,);
        if frequency {
            for key in locus_keys {
                print!(
                    " {}",
                    locus_map[&key][0] / (locus_map[&key][1] + locus_map[&key][0]),
                );
            }
            println!();
        } else {
            for key in locus_keys {
                print!(" {:.0},{:.0}", locus_map[&key][0], locus_map[&key][1],);
            }
            println!();
        }
    }
}

pub fn locus_freqs(
    pop_map: &HashMap<i32, (String, String)>,
    q_vec: rust_htslib::bcf::record::BufferBacked<
        std::vec::Vec<&[f32]>,
        rust_htslib::bcf::record::Buffer,
    >,
) -> HashMap<String, Vec<f32>> {
    //turn here down into a function?
    let mut locus_map = freq_map(&pop_map);

    for pop in q_vec.iter().enumerate() {
        let pop_gp = p_vec(pop.1);
        //only bi-allelic sites!!
        if pop_gp.len() == 3 {
            //println!("{} {:?}", pop.0, pop_gp);
            //get the pop. of the current ind.
            let ind_idx = &(pop_map[&(pop.0 as i32)]).1;
            let freq_up = locus_map
                .entry(ind_idx.to_string())
                .or_insert(vec![0.0, 0.0]);
            //count alleles while accounting for genotype probs
            freq_up[0] += 2.0 * pop_gp[0] + 1.0 * pop_gp[1];
            freq_up[1] += 2.0 * pop_gp[2] + 1.0 * pop_gp[1];
        }
    }
    locus_map
}

pub fn locate_gp(format: &str) -> usize {
    //let a = "GT:PL:DP:GP:AD:GQ";
    let index = format
        .split(":")
        .collect::<Vec<&str>>()
        .iter()
        .position(|&r| r == "GP")
        .unwrap();
    index
}

//calculate the vector of genotype probs, normalize to sum to 1 (why don't they already?!?!)
pub fn p_vec(q_vec: &[f32]) -> Vec<f32> {
    //Q = -10*log10(P)
    //P = 10^(-Q/10)
    //let q_split: Vec<&str> = q_string.split(",").collect();
    if q_vec.len() == 3 {
        let ps: Vec<f32> = q_vec.iter().map(|x| get_p(*x)).collect();
        let p_sum: f32 = ps.iter().sum();
        let ps_norm = ps.iter().map(|x| x / p_sum).collect();
        ps_norm
    } else {
        vec![0.0, 0.0, 0.0]
    }
}

//convert the phred-scaled genotype probs to raw probs
pub fn get_p(q: f32) -> f32 {
    let base: f32 = 10.0;

    let p: f32 = base.powf(-q / 10.0);

    //DROPPING PROBS > 1
    if p > 1.0 {
        0.0
    } else {
        p
    }
}

//split str on white space into a str vector
pub fn line_splitter(line: &str) -> Vec<&str> {
    let line_vec: Vec<&str> = line.split_whitespace().collect();
    line_vec
}

//split the genotype format string on colon, return the genotype probs
//default index position is 4
pub fn get_qstring(line: &str, gp_pos: usize) -> &str {
    let line_vec: Vec<&str> = line.split(':').collect();
    line_vec[gp_pos]
}

/*
UNIT TESTS
*/

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_zero_vec() {
        let zero_probs = vec![0.0, 0.0, 0.0];
        let sum_probs: f32 = p_vec(&zero_probs).iter().sum();
        let one = 1.0f32;
        assert_eq!(sum_probs, one)
    }

    #[test]
    fn test_p_vec() {
        let probs = vec![4.0, 3.0, 6.0];
        let sum_probs: f32 = p_vec(&probs).iter().sum();
        assert_eq!(sum_probs, 1.0)
    }

    #[test]
    #[should_panic(expected = "Unable to open file")]
    fn missing_vcf() {
        make_freqs("", "", true);
    }

    #[test]
    #[should_panic(expected = "Unable to open file")]
    fn missing_popkey() {
        process_popkey("");
    }

    #[test]
    #[should_panic(expected = "VCF site has the incorrect number of fields: 1")]
    fn bad_vcf() {
        let test_vcf = "example_data/fake.vcf.gz";
        let testkey = "example_data/pop_key.txt";
        let the_pop_key = process_popkey(testkey);
        let reader = Reader::from_path(test_vcf).expect("Error opening file.");
        //let reader = BufReader::new(BGZFReader::new(f));
        process_vcf(reader, the_pop_key, true);
    }
}
