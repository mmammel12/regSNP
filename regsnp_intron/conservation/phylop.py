#!/usr/bin/env python
import argparse
import os.path

import numpy as np
from bx.bbi.bigwig_file import BigWigFile


class Phylop(object):
    def __init__(self, bw_fname):
        """
        :param bw_fname: Phylop 100way bigwig file name.
        """
        self.bw_handle = open(os.path.expanduser(bw_fname))
        self.bw = BigWigFile(self.bw_handle)

    def get(self, chrom, start, end, flanking=0):
        """
        :param chrom: chr1, chr2, etc.
        :param start: 0-based.
        :param end: 1-based.
        :param flanking: length of flanking sequence on each side.
        """
        return np.nanmean(self.bw.get_as_array(chrom, start - flanking, end + flanking))

    def calculate(self, fname, out_fname):
        """
        :param fname: SNP BED.
        :param out_fname: output file.
        """
        with open(fname) as bed_f, open(out_fname, 'w') as out_f:
            out_f.write('\t'.join(self._build_header()) + '\n')
            for line in bed_f:
                cols = line.rstrip().split('\t')
                chrom, start, end = cols[:3]
                start = int(start)
                end = int(end)
                scores = [self.get(chrom, start, end), self.get(chrom, start, end, 3), self.get(chrom, start, end, 7)]
                out_f.write('\t'.join(map(str, cols + scores)) + '\n')

    def close(self):
        self.bw_handle.close()

    def _build_header(self):
        header = ['#chrom_snp', 'start_snp', 'end_snp', 'ref', 'alt', 'feature', 'gene_id',
                  'chrom', 'start', 'end', 'name', 'score', 'strand', 'distance']
        header += ['phylop1',
                   'phylop3',
                   'phylop7']
        return header


def main():
    parser = argparse.ArgumentParser(description='''
            Given SNP Bed file and Phylop 100 way bigwig file, calculate the
            Phylop score around SNP''')
    parser.add_argument('pfname',
            help='Phylop 100 way bigwig file')
    parser.add_argument('bfname',
            help='BED file')
    parser.add_argument('ofname',
            help='output file')
    args = parser.parse_args()
    phylop = Phylop(args.pfname)
    phylop.calculate(args.bfname, args.ofname)
    phylop.close()

if __name__ == '__main__':
    main()
