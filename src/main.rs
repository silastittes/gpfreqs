mod funs;
use clap::{App, Arg};
use funs::make_freqs;

fn main() {
    let matches = App::new("gpfreqs")
        .version("0.0.1")
        .author("Silas Tittes <silas.tittes@gmail.com>")
        .about("Use genotype probabilities in a VCF to calculate designated population allele frequencies.")
        //.arg(Arg::new("config")
        //    .short('c')
        //    .long("config")
        //    .value_name("FILE")
        //    .about("Sets a custom config file")
        //    .takes_value(true))
        .arg(
            Arg::new("frequency")
                .short('f')
                .takes_value(false)
                .required(false)
                .about("Returns reference allele frequency rather than ref alt counts."),
        )
        .arg(
            Arg::new("gp_index")
                .short('g')
                .takes_value(true)
                .required(true) 
                .about("The 0-index position of the three genotype probabilites (GP) in the FORMAT field (eigth column) of the input VCF.
For example if FORMAT is, 'GT:PL:DP:AD:GP:GQ', input would be -g 4."),
        )
        .arg(
            Arg::new("vcf")
                .short('v')
                .takes_value(true)
                .required(true)
                .about("Path to the vcf input file."),
            //.index(1),
        )
        .arg(
            Arg::new("popkey")
                .short('p')
                .required(true)
                .takes_value(true)
                .about(
                    "File containing population information. 
First three columns must be:
- a zero-based index of each individuals position in the vcf
(index starts at the first sample, skipping the first 10 fields).
- the name of each individual as it appears in the vcf file.
- an ID for which population each individual belongs to. 
                    
Must be whitespace separated and without header names.
for example, a file could be a sample as:
                    
0 individual1 pop1
",
                ),
        )
        .get_matches();

    //let vcf_name = "example_data/small.vcf.gz";
    let vcf_name = matches.value_of("vcf").unwrap();
    //let pop_name = "example_data/pop_key.txt";
    let pop_name = matches.value_of("popkey").unwrap();

    let gp_index = matches.value_of("gp_index").unwrap().to_string().parse::<usize>().unwrap();

    if matches.is_present("frequency") {
        make_freqs(vcf_name, pop_name, true, gp_index)
    } else {
        make_freqs(vcf_name, pop_name, false, gp_index)
    };
}
