# gvcf_to_vcf

This program will convert a GVCF file to VCF format. The GVCF format batches homozygous reference calls with a set minimum DP value. This script will replace the GVCF batch with individual calls for each base. 

GVCF includes a value called "<NON_REF>" to represent allele calls that aren't reference calls. These may or may not be named ALT variants. This script will remove these and adjust FORMAT or INFO fields that depend on the number of ALT calls.

Additionally, in some GVCF files, GT is not calculated. This script will optionally re-call GT values based upon AD values. The GT caller is very simple and based on minimum count and AF thresholds. By default, at least one read with a minimum AF of 10% is required. This means that if the total depth is less than 11, any single read will be enough to call a genotype.

## Installation

This script expects to run in the context of virtualenv container (venv). In particular, the script requires the installation of the `pysam` library. The library is used to read reference values from a genomic FASTA file. The script *can* run in any environment that has this library installed. It *should* be run in the context of a venv, so the library installation is separated from the global Python installation.

Example setup:

    $ python3 -m venv .venv
    $ .venv/bin/pip install -r requirements.txt
    $ .venv/bin/python3 gvcf_to_vcf.py {args} input.g.vcf.gz ref.fasta | bgzip > new.vcf.gz
