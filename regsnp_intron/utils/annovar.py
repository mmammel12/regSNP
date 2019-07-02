#!/usr/bin/env python
import argparse
import os.path
import re
import subprocess


class Annovar(object):
    def __init__(self, annovar_path, db_path):
        self.annovar_path = os.path.expanduser(annovar_path)
        self.db_path = os.path.expanduser(db_path)

    def annotate(self, in_fname, out_fname):
        in_fname = self._prepare_input(in_fname, out_fname + '.avinput')
        subprocess.check_call([os.path.join(self.annovar_path, 'table_annovar.pl'),
                               in_fname,
                               self.db_path,
                               '-buildver',
                               'hg19',
                               '-out',
                               out_fname,
                               '-remove',
                               '-protocol',
                               'ensGene',
                               '-operation',
                               'g',
                               '-nastring',
                               '.'])
        with open(out_fname + '.hg19_multianno.txt') as in_f, open(out_fname + '.intronic', 'w') as out_f:
            header = in_f.readline()
            for line in in_f:
                cols = line.rstrip().split('\t')
                if re.search(r'intronic|splicing', cols[5]) and not re.search(r'exon', cols[5]):
                    out_f.write('\t'.join(map(str, [cols[0], int(cols[1]) - 1] + cols[2:7])) + '\n')

    def _prepare_input(self, in_fname, out_fname):
        in_fname = os.path.expanduser(in_fname)
        out_fname = os.path.expanduser(out_fname)
        with open(in_fname) as in_f, open(out_fname, 'w') as out_f:
            for line in in_f:
                chrom, pos, ref, alt = line.rstrip().split('\t')
                out_f.write('\t'.join([chrom, pos, pos, ref, alt]) + '\n')
        return out_fname


def main():
    parser = argparse.ArgumentParser(description='''
            Given SNP input, invoke Annovar and generate ensGene annotation.''')
    parser.add_argument('annovar_path',
            help='path to annovar home folder')
    parser.add_argument('db_path',
            help='path to annovar db')
    parser.add_argument('ifname',
            help='input file name')
    parser.add_argument('ofname',
            help='output file name')
    args = parser.parse_args()
    annovar = Annovar(args.annovar_path, args.db_path)
    annovar.annotate(args.ifname, args.ofname)

if __name__ == '__main__':
    main()
