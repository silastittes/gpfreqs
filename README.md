## `gpfreqs` – Calculate designated population allele counts (or frequencies) based on genotype probabilities in a VCF.

# Installation

`gofreqs` is written in Rust. While I may eventually turn it into a crate, installation must be done locally for now. 

[Assuming Rust in installed](https://www.rust-lang.org/tools/install), installing `gpfreqs` is as simple as 

```
git clone https://github.com/silastittes/gpfreqs.git
cd gpfreqs/
cargo build --release
```

The compiled executable will be available as `target/release/gpfreqs`

# Documentation

```
gpfreqs 0.0.1

Silas Tittes <silas.tittes@gmail.com>

Use genotype probabilities in a VCF to calculate designated population allele frequencies.

USAGE:
    gpfreqs [FLAGS] -g <gp_index> -v <vcf> -p <popkey>

FLAGS:
    -f               Returns reference allele frequency rather than ref alt counts.
    -h, --help       Print help information
    -V, --version    Print version information

OPTIONS:
    -g <gp_index>        The 0-index position of the three genotype probabilites (GP) in the FORMAT
                         field (eigth column) of the input VCF.
                                         For example if FORMAT is, 'GT:PL:DP:AD:GP:GQ', input would
                         be -g 4.
    -p <popkey>          File containing population information.
                         First three columns must be:
                         - a zero-based index of each individuals position in the vcf
                         (index starts at the first sample, skipping the first 10 fields).
                         - the name of each individual as it appears in the vcf file.
                         - an ID for which population each individual belongs to.
                         
                         Must be whitespace separated and without header names.
                         for example, a file could be a sample as:
                         
                         0 individual1 pop1
    -v <vcf>             Path to the vcf input file.
```

The vcf file can (and should) be compressed.

