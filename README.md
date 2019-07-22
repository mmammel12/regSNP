# regsnp-intron

regsnp-intron predicts the disease-causing probability of intronic single nucleotide variants (iSNVs) based on both genomic and protein structural features.

## Prerequisites

**Python (>= 2.7.11):**
(Python 3 is not currently supported.)

The following Python libraries are also required. They will be automatically installed if you use pip (see [Installation](#Installation)).

- Pandas (== 0.17.1),
- pysam (== 0.15.2),
- pymongo (>= 3.8.0)

## Installation

1. To install, you need to install all the required Python libraries first:

```bash
pip install pandas==0.17.1
pip install pysam==0.15.2
pip install pymongo
```

2. Then run the following command to clone the repository.

```bash
git clone https://github.com/mmammel12/regSNP.git
```
3. You will also need a human genome file, run the following commands to download then extract it

```bash
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gunzip hg19.fa.gz
```

2. Make sure you have the settings file set up correctly (see [Configuration](#Configuration))

## Configuration

1. Modify the data in settings/settings.json file. Type `regsnp_intron --help` to find the location of default settings.json file. You can also provide customized settings.json file with `-s` (see [Usage](#Usage)):

```json
{
  "dbURI": "MONGO_DB_URI",
  "dbUsername": "MONGO_DB_USERNAME",
  "dbPassword": "MONGO_DB_PASSWORD",
  "hg_dir": "PATH/TO/HUMAN_GENOME_FILE"
}
```

## Usage

from the regSNP directory run the following commands:

```bash
cd regsnp_intron
python regsnp_intron.py [-h] [-s SFNAME] [-f] ifname out_dir
```

Given a list of intronic SNVs, predict the disease-causing probability based
on genomic and protein structural features.

positional arguments:
  ifname                input SNV file. Contains four columns: chrom, pos, ref, alt.
  out_dir               directory contains output files

optional arguments:
  -h, --help            show this help message and exit
  -s SFNAME, --sfname SFNAME
                        JSON file containing settings. Default setting file
                        located at: regsnp_intron/settings/settings.json
  -f, --force           overwrite existing directory


## Output

The following files will be generated under the output directory:

- snp.prediction.txt: tab-delimited text file containing prediction results and all the features for iSNVs.
- invalid.txt: invalid input from input file, will not exist if all data is valid
- tmp: temporary folder containing all the intermediate results (can be deleted).

snp.prediction.txt contains the following columns:

```
chrom: Chromosome
pos: Position.
ref: Reference allele.
alt: Alternative allele.
disease: Categorical prediction.
prob: Disease-causing probability [D, PD, B]. Higher score indicates higher probability of being pathologic.
splicing_site: Indicates on/off splicing site. Splicing sites are defined as +7bp from donor site and -13bp from acceptor site.
features: The rest of columns contain all the genomic and protein structural features around each iSNV.
```
