#!/usr/bin/env python
import argparse
import os.path

from pybedtools import BedTool


class ClosestExon(object):
    def __init__(self, efname):
        self.efname = os.path.expanduser(efname)

    def get_closest_exon(self, in_fname, out_fname, idx_bstart=8, idx_distance=13):
        """
        :param in_fname: input SNV files in BED-like format.
        Default assume 7 columns: chrom, start (0-based), end (1-based), ref, alt, region (intronic/exonic), gene_id
        :param out_fname:
        :param idx_bstart: column index of start position from file b in output. '-1' indicates no closest exon found.
        :param idx_distance: column index of distance in output. '0' indicates overlap.
        """
        in_f = BedTool(in_fname)
        #  Annotate the closest exons and filter out the SNVs that overlap with exons and SNVs without closest exons.
        BedTool(in_f.sort().closest(self.efname, D='b', t='first')).\
            filter(lambda x: x[idx_bstart] != '-1' and x[idx_distance] != '0').\
            moveto(out_fname)


def main():
    parser = argparse.ArgumentParser(description='''
            Given BED file and exon annotation, find the closest exon.''')
    parser.add_argument('efname',
            help='exon bed file')
    parser.add_argument('ifname',
            help='input bed file')
    parser.add_argument('ofname',
            help='output bed file')
    args = parser.parse_args()
    closest_exon = ClosestExon(args.efname)
    closest_exon.get_closest_exon(args.ifname, args.ofname)

if __name__ == '__main__':
    main()
