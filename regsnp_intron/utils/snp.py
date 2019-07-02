#!/usr/bin/env python
import argparse
import logging
import os.path

import pandas as pd
import pysam

class SNP(object):
    def __init__(self, in_fname):
        self.in_fname = os.path.expanduser(in_fname)
        self.logger = logging.getLogger(__name__)

    def sort(self, out_fname):
        in_f = pd.read_csv(self.in_fname, sep='\s+', names=['#chrom', 'pos', 'ref', 'alt'], skipinitialspace=True)
        in_f.sort_values(by=['#chrom', 'pos']).to_csv(out_fname, sep='\t', header=False, index=False)

    def is_valid(self, ref_fname):
        """
        An SNV input file is valid if:
            1. contains four columns: chrom, pos, ref, alt
            2. chrom and pos are in reference genome.
            3. ref and alt are single nucleotide allele
        :param ref_fname: reference fasta file with .fai index.
        :return: a boolean value.
        """
        if os.path.getsize(self.in_fname) <= 0:
            self.logger.error('Input file is empty.')
            return False
        with open(self.in_fname) as f, pysam.FastaFile(ref_fname) as ref:
            references = dict(zip(ref.references, ref.lengths))
            alleles = {'A', 'C', 'G', 'T'}
            for line in f:
                cols = line.strip().split()
                if len(cols) != 4:
                    self.logger.error('Invalid file format: {0}. '
                                      'File should contain tab-delimited four columns: chrom, pos, ref, alt'.
                                      format(line.strip()))
                    return False
                chrom, pos, ref_allele, alt_allele = cols
                if chrom not in references:
                    self.logger.error('Chromosome cannot be found in the reference: {0}.'.format(line.strip()))
                    return False
                if ref_allele not in alleles or alt_allele not in alleles:
                    self.logger.error('Ref and alt should be in {{A, C, G, T}}: {0}'.format(line.strip()))
                    return False
                try:
                    pos = int(pos)
                except ValueError:
                    self.logger.error('Coordinate must be integer: {0}.'.format(line.strip()))
                    return False
                if pos > references[chrom] or pos < 0:
                    self.logger.error('Coordinate exceed the range of chromosome: {0}.'.
                                      format(line.rstrip()))
                    return False
                if ref_allele != ref.fetch(reference=chrom, start=(pos - 1), end=pos).upper():
                    if alt_allele != ref.fetch(reference=chrom, start=(pos - 1), end=pos).upper():
                        self.logger.error('Input reference allele is not consistent with reference genome: {0}.'.
                                          format(line.strip()))
                        return False
                    else:
                        self.logger.warning('Input reference allele and alternative allele are switched: {0}.'
                                            .format(line.strip()))
        return True

    def switch_alleles(self, in_fname, ref_fname, out_fname):
        """
        Switch ref_allele and alt_allele in input file if the alt_allele matches reference genome.
        :param in_fname: input snp
        :param ref_fname: reference fasta file with .fai index
        :param out_fname: output snp
        :return:
        """
        with open(in_fname) as in_f, pysam.FastaFile(ref_fname) as ref, open(out_fname, 'w') as out_f:
            for line in in_f:
                cols = line.strip().split()
                chrom, pos, ref_allele, alt_allele = cols
                pos = int(pos)
                if ref_allele != ref.fetch(reference=chrom, start=(pos - 1), end=pos).upper() and \
                    alt_allele == ref.fetch(reference=chrom, start=(pos - 1), end=pos).upper():
                    tmp = ref_allele
                    ref_allele = alt_allele
                    alt_allele = tmp
                out_f.write('\t'.join(map(str, [chrom, pos, ref_allele, alt_allele])))
                out_f.write('\n')


def main():
    parser = argparse.ArgumentParser(description='''
            Given input SNP file, prepare for feature calculator.''')
    parser.add_argument('ifname',
            help='input snp file, contains four columns: chrom, pos, ref, alt')
    parser.add_argument('ofname',
            help='output snp file')
    args = parser.parse_args()
    snp = SNP(args.ifname)
    snp.sort(args.ofname)

if __name__ == '__main__':
    main()
