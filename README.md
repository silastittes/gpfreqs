## `gpfreqs` -- Calculate designated population allele counts (or frequencies) based on genotype probabilities in a VCF.

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
Use genotype probabilities in a VCF to calculate population genetic statistics.

USAGE:
    gpfreqs [FLAGS] -v <vcf> -p <popkey>

FLAGS:
    -f               Returns reference allele frequency rather than ref alt counts.
    -h, --help       Print help information
    -V, --version    Print version information

OPTIONS:
    -p <popkey>        File containing population information.
                       First three columns must be:
                       - a zero-based index of each individuals position in the vcf
                       (index starts at the first sample, skipping the first 10 fields).
                       - the name of each individual as it appears in the vcf file.
                       - an ID for which population each individual belongs to.
                       
                       Must be whitespace separated and without header names.
                       for example, a file could be a sample as:
                       
                       0 individual1 pop1
    -v <vcf>           Path to the vcf input file.
```


# Example

```
target/release/gpfreqs -v example_data/small.vcf.gz -p example_data/pop_key.txt -f
```

The vcf file can (and should) be compressed.

