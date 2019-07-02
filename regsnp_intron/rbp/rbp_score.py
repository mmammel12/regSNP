#!/usr/bin/env python
import collections
import os.path

import rbp


class RBPScore(object):
    """
    Store a collection of PSSM files. Given a file contains the input sequences, calculate the matching scores for all
    the PSSMs.
    """
    def __init__(self, pssm_path, pssm_list_fname):
        self.pssm_path = os.path.expanduser(pssm_path)  # the path of PSSM files
        self.pssm_list = self.parse_pssm_list(pssm_list_fname)  # a list of valid PSSM files
        self.pssm_num = len(self.pssm_list)
        self.pssms = self.generate_pssms(self.pssm_path, self.pssm_list)

    def parse_pssm_list(self, fname):
        """
        parse the list of valid PSSM file.
        :param fname: each row contains a valid PSSM file name.
        :return: a list of PSSM file name.
        """
        pssm_list = []
        with open(os.path.expanduser(fname)) as f:
            for line in f:
                pssm_list.append(line.rstrip())
        return pssm_list

    def generate_pssms(self, pssm_path, pssm_list):
        """

        :param pssm_path: the path of PSSM files.
        :param pssm_list: a list of valid PSSM files.
        :return: a collections.OrderedDict of RBPs.
        """
        pssms = collections.OrderedDict()
        for pssm_fname in pssm_list:
            pssms[pssm_fname] = rbp.RBP(os.path.join(pssm_path, pssm_fname))
        return pssms

    def cal_matching_score(self, seq_fname, out_fname, keep_ncol=6):
        """
        given a sequence file that contains reference sequences and alternative sequences, calculate the max matching
        score of each RBP for each sequence.
        :param seq_fname: a file whose first six columns are: chrom, pos, ref, alt, ref_seq (the center is the SNP),
        alt_seq (the center is the SNP).
        :param out_fname: output file name
        :param keep_ncol: keep only the first keep_ncol columns
        :return:
        """
        with open(os.path.expanduser(seq_fname)) as seq_f, open(os.path.expanduser(out_fname), 'w') as out_f:
            header = seq_f.readline().rstrip().split('\t')
            header = header[:keep_ncol]
            header = header + [x + '_ref' for x in self.pssm_list] + [x + '_alt' for x in self.pssm_list]
            out_f.write('\t'.join(header) + '\n')
            for line in seq_f:
                cols = line.rstrip().split('\t')
                chrom, pos, ref, alt, ref_seq, alt_seq = cols[:keep_ncol]
                ref_scores = []
                alt_scores = []
                for pssm in self.pssms:
                    ref_scores.append(self.pssms[pssm].match(ref_seq))
                    alt_scores.append(self.pssms[pssm].match(alt_seq))
                out_f.write('\t'.join(map(str,
                                          [chrom, pos, ref, alt, ref_seq, alt_seq] + ref_scores + alt_scores)) + '\n')


def main():
    pass

if __name__ == '__main__':
    main()
