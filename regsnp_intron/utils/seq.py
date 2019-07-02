#!/usr/bin/env python
import argparse
import logging
import os.path
import string

import pysam


class FlankingSeq(object):
    def __init__(self, rfname, length):
        self.ref = pysam.Fastafile(os.path.expanduser(rfname))
        self.length = length
        self.logger = logging.getLogger(__name__)

    def fetch_flanking_seq(self, in_fname, out_fname, otype='tab'):
        with open(in_fname) as in_f, open(out_fname, 'w') as out_f:
            if otype == 'tab':
                out_f.write('\t'.join(['#chrom', 'pos', 'ref', 'alt', 'ref_seq', 'alt_seq']) + '\n')
            for line in in_f:
                cols = line.rstrip().split('\t')
                chrom, start, pos, ref, alt = cols[:5]
                strand = cols[12]
                pos = int(pos)
                seq, seq_alt = self._fetch_seq(chrom, pos, ref, alt)
                assert strand == '+' or strand == '-'
                if strand == '-':
                    seq = self._reverse_complement(seq)
                    seq_alt = self._reverse_complement(seq_alt)
                if otype == 'tab':
                    out_f.write('\t'.join(map(str, [chrom, pos, ref, alt, seq, seq_alt])) + '\n')
                elif otype == 'fasta':
                    out_f.write('>' + ':'.join([chrom, str(pos), ref, alt]) + '\n')
                    out_f.write(seq + '\n')
                    out_f.write('>' + ':'.join([chrom, str(pos), ref, alt]) + '\n')
                    out_f.write(seq_alt + '\n')

    def _fetch_seq(self, chrom, pos, ref_allele, alt_allele):
        pos = pos
    
        start = (pos - 1) - self.length  # 0-based. "chr:(pos - 1) ~ pos" is the coordinate of the SNP.
        end = pos + self.length
        seq = self.ref.fetch(reference=chrom, start=start, end=end).lower()
    
        if seq[self.length:self.length + 1].lower() != ref_allele.lower():
            self.logger.error('The input reference allele does not match the reference fasta:\n')
            self.logger.error('\t'.join([chrom, str(pos), ref_allele, alt_allele,
                                         seq[self.length:self.length + 1]]) + '\n')

        # Replace the center of the sequence with the alternative allele.
        seq_alt = seq[:self.length] + alt_allele + seq[self.length + 1:]
    
        return seq, seq_alt

    def _reverse_complement(self, seq):
        trans_tab = string.maketrans('ACGTacgt', 'TGCAtgca')
        seq_comp = seq.translate(trans_tab)
        seq_rev_comp = seq_comp[::-1]
        return seq_rev_comp

    def close(self):
        self.ref.close()


def main():
    parser = argparse.ArgumentParser(description='''
            Given SNP file, fetch flanking sequence.''')
    parser.add_argument('rfname',
            help='reference fasta file. The .fai index should be in the same folder.')
    parser.add_argument('bfname',
            help='SNP input file')
    parser.add_argument('ofname',
            help='output file')
    parser.add_argument('-l', '--length', type=int, default=20,
            help='flanking sequence length')
    args = parser.parse_args()
    seq = FlankingSeq(args.rfname, args.length)
    seq.fetch_flanking_seq(args.bfname, args.ofname)
    seq.close()

if __name__ == "__main__":
    main()
