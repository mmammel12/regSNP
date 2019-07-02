#!/usr/bin/env python
import os.path

import numpy as np


class RBP(object):
    """
    Store PSSM of single RBP. Given a sequence, calculate the max matching score.
    """
    def __init__(self, pssm_fname):
        self.name = os.path.basename(pssm_fname)  # pssm file name
        self.pssm = self.parse_pssm(pssm_fname)   # pssm numpy ndarray. A 4 x m matrix where m is the motif length.
                                                  # Row 0,1,2,3 correspond to A,C,G,T.
        self.len = self.pssm.shape[1]  # motif length

    def parse_pssm(self, fname):
        """
        parse PSSM matrix file
        :param fname: a file contains a 4 x m PSSM matrix. Each row corresponds to A,C,G,T, and m is the motif length.
        :return: a 4 x m numpy ndarray.
        """
        matrix = np.genfromtxt(fname)
        return matrix

    def match(self, seq, cover_center=True):
        """
        calculate the max matching score.
        :param seq: input DNA sequence. Should be the same strand as RNA. The length should be longer than RBP motif.
        :param cover_center: indicate the motif must cover the center of sequence.
        :return: the max matching score.
        """
        seq = seq.upper()
        seq_len = len(seq)

        if cover_center:
            assert seq_len > self.len * 2, 'sequence length ({0}) should be longer than 2 * motif length ({1}).'\
                .format(seq_len, self.len)
            assert seq_len % 2 == 1, 'sequence should be up_flanking + mid + down_flanking, the length should be even.'
            center_idx = seq_len // 2
            seq = seq[(center_idx - self.len + 1):(center_idx + self.len)]
            seq_len = len(seq)

        assert seq_len >= self.len, 'sequence length {0} is shorter than motif length {1}.'.format(seq_len, self.len)
        max_score = float('-inf')

        for i in range(0, seq_len - self.len + 1):  # Use a sliding window.
            sub_seq = seq[i:(i + self.len)]
            score = self._match(sub_seq)
            max_score = max(max_score, score)

        return max_score

    def _match(self, seq, pssm=None):
        """
        calculate the matching score.
        :param seq: input DNA sequence. Should be the same strand as RNA. The length should be equal to RBP motif.
        :param pssm: a 4 x m numpy ndarray represents PSSM. Each row corresponds to A,C,G,T.
        :return: the matching score
        """
        seq_len = len(seq)
        assert seq_len == self.len, 'sequence length {0} should equal to motif length {1}.'.format(seq_len, self.len)
        letters = {'A':0, 'C':1, 'G':2, 'T':3}  # the row index of each letter in PSSM.
        score = 0

        if not pssm:
            pssm = self.pssm

        for i, v in enumerate(seq):
            if v in letters:
                current_score = pssm[letters[v], i]
            else:
                current_score = 0.0  # 'N'
            score += current_score

        return score


def main():
    pass

if __name__ == '__main__':
    main()
