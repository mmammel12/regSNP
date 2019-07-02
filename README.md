# regsnp-intron
regsnp-intron predicts the disease-causing probability of intronic single nucleotide variants (iSNVs) based on both genomic and protein structural features.

## Prerequisites
**ANNOVAR (>= 2016Feb01):**
Follow the instructions at <http://annovar.openbioinformatics.org/en/latest> to install, and prepare Ensembl gene annotation.
```bash
tar -xf annovar.latest.tar.gz
cd annovar
perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar ensGene humandb/
```
**BEDTools (>= 2.25.0):**
Follow the instructions at <http://bedtools.readthedocs.io/en/latest> install, and make sure the programs are in your PATH.

**Python (>= 2.7.11):**
Installing libraries such as Numpy and Scipy can be a little difficult for inexperienced users. We highly recommend installing [Anaconda](https://docs.continuum.io/anaconda). Anaconda conveniently installs Python and other commonly used packages for scientific computing and data science. (Python 3 is not currently supported.)

The following Python libraries are also required. They will be automatically installed if you use pip (see [Installation](#Installation)).

* Numpy (>= 1.10.4),
* Scipy (>= 0.17.0),
* Pandas (>= 0.17.1),
* Scikit-learn (>= 0.17),
* pysam (>= 0.8.4),
* bx-python (0.7.3),
* pybedtools (>= 0.7.6)

## Installation
1. The easiest way is to install with pip (recommended). This will also install all the required Python libraries: 
```bash
pip install regsnp-intron
```
2. To manually install, you need to install all the required Python libraries first. Then, use the following command:
```bash
git clone https://github.com/linhai86/regsnp_intron.git
cd regsnp_intron
python setup.py install
```

## Configuration
1. Download the database and annotation files for human genome (hg19):
```bash
wget to_be_added
unzip db.zip
```
2. Modify the paths in settings/settings.json file. Type `regsnp_intron --help` to find the location of default settings.json file. You can also provide customized settings.json file with `-s` (see [Usage](#Usage)):
```json
{
    "annovar_path": "/path/to/annovar",
    "annovar_db_path": "/path/to/annovar/humandb",
    "db_dir": "/path/to/db"
}
```

## Usage
```bash
usage: regsnp_intron [-h] [-s SFNAME] [-f] ifname out_dir

Given a list of intronic SNVs, predict the disease-causing probability based
on genomic and protein structural features.

positional arguments:
  ifname                input SNV file. Contains four columns: chrom, pos,
                        ref, alt.
  out_dir               directory contains output files

optional arguments:
  -h, --help            show this help message and exit
  -s SFNAME, --sfname SFNAME
                        JSON file containing settings. Default setting file
                        locate at: regsnp_intron/settings/settings.json
  -f, --force           overwrite existing directory

```

## Output
The following files will be generated under the output directory:

* snp.prediction.txt: tab-delimited text file containing prediction results and all the features for iSNVs.
* snp.features.txt: tab-delimited text file containing all the features for iSNVs (can be deleted).
* tmp: temporary folder containing all the intermediate results (can be deleted).

snp.prediction.txt contains the following columns:
```
chrom: Chromosome
pos: Position.
ref: Reference allele.
alt: Alternative allele.
disease: Categorical prediction.
prob: Disease-causing probability [0, 1]. Higher score indicates higher probability of being pathologic.
splicing_site: Indicates on/off splicing site. Splicing sites are defined as +7bp from donor site and -13bp from acceptor site.
features: The rest of columns contain all the genomic and protein structural features around each iSNV.
. 
.
.
```

## Citation
Please cite:
```
To be added
```
