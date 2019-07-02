#!/usr/bin/env python
import argparse
import logging
import os.path


class VCF(object):
    def __init__(self, in_fname):
        self.in_fname = os.path.expanduser(in_fname)
        self.logger = logging.getLogger(__name__)

    def convert_to_txt(self, out_fname, pass_filter=True):
        """
        convert vcf to 4-column txt input file for regsnp_intron.
        :param out_fname: output txt file name
        :param pass_filter: whether filter out variants that failed the filter
        """
        with open(self.in_fname) as vcf, open(os.path.expanduser(out_fname), 'w') as out_f:
            for line in vcf:
                if line.find('#') == 0:
                    continue
                cols = line.rstrip().split('\t')
                chrom, pos, rs_id, ref, alts, qual, filters, info, vcf_format = cols[:9]
                alts = alts.split(',')
                filters = filters.split(',')
                if not pass_filter or not filters or filters == ['PASS']:
                    for alt in alts:
                        if len(ref) == 1 and len(alt) == 1:
                            out_f.write('\t'.join(map(str, [chrom, pos, ref, alt])) + '\n')


def main():
    parser = argparse.ArgumentParser(description='''
            Convert a vcf file to a 4-column txt file: chrom, pos, ref, alt.''')
    parser.add_argument('ifname',
            help='input vcf file')
    parser.add_argument('ofname',
            help='output txt file')
    args = parser.parse_args()
    vcf = VCF(args.ifname)
    vcf.convert_to_txt(args.ofname)

if __name__ == '__main__':
    main()
