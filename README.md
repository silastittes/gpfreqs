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
gpfreqs 0.0.4

Silas Tittes <silas.tittes@gmail.com>

Use genotype probabilities in a VCF to calculate designated population allele frequencies.

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
    -v <vcf>           Path to the vcf input file. Can gzipped (File should end in .gz) or
                       uncompressed.
```

# Example

```
target/release/gpfreqs -f -v example_data/small.vcf.gz -p example_data/pop_key.txt | less -S
```

This should return

```
contig position pop1 pop2 pop3 
Super-Scaffold_48 109  0.5 0.48 0.5
Super-Scaffold_48 177  0.45918366 0.48 0.43939394
Super-Scaffold_48 193  0.4489796 0.44 0.43939394
Super-Scaffold_48 210  NaN NaN NaN
Super-Scaffold_48 219  0.47959185 0.5 0.4848485
Super-Scaffold_48 285  0.3265306 0.4 0.3939394
Super-Scaffold_48 290  0.45918366 0.5 0.5
Super-Scaffold_48 325  0.39795917 0.42 0.37878788
Super-Scaffold_48 352  0.68601024 0.6292135 0.700375
```

The input VCF file can (and should) be compressed.

