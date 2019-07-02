#!/usr/bin/env python
import argparse
import numpy
import pysam
import string


class JunctionStrength(object):
    def __init__(self, donor_ic, acceptor_ic):
        self.donor_ic_matrix = self.parse_matrix(donor_ic)
        self.acceptor_ic_matrix = self.parse_matrix(acceptor_ic)

    def parse_matrix(self, fname):
        matrix = []
        with open(fname) as f:
            header = f.readline()
            for line in f:
                cols = line.rstrip().split('\t')
                matrix.append(map(float, cols[1:]))
        matrix = numpy.array(matrix)
        return matrix.T
    
    def cal_donor_ic(self, seq):
        assert len(seq) == self.donor_ic_matrix.shape[1]
        seq = seq.upper()
        letter = {'G': 0, 'A': 1, 'C': 2, 'T': 3}
        score = 0
        for i, v in enumerate(seq):
            score += self.donor_ic_matrix[letter[v], i]
        return score

    def cal_acceptor_ic(self, seq):
        assert len(seq) == self.acceptor_ic_matrix.shape[1]
        seq = seq.upper()
        letter = {'G': 0, 'A': 1, 'C': 2, 'T': 3}
        score = 0
        for i, v in enumerate(seq):
            score += self.acceptor_ic_matrix[letter[v], i]
        return score

    def cal_junction_strength(self, ref_fname, bed_fname, out_fname):
        with open(bed_fname) as bed_f, open(out_fname, 'w') as out_f:
            out_f.write('\t'.join(self._build_header()) + '\n')
            ref = pysam.Fastafile(ref_fname)
            for line in bed_f:
                cols = line.rstrip().split('\t')
                chrom_snp, start_snp, end_snp, ref_snp, alt_snp, feature, gene_id,\
                    chrom, start, end, transcript_id, score, strand, distance = cols
                start_snp, end_snp, start, end, distance = map(int, [start_snp, end_snp, start, end, distance])
                if strand == '+':
                    astart = start - 13  # acceptor: 3' end of intron
                    aend = start + 1
                    dstart = end - 3  # donor: 5' end of intron
                    dend = end + 7
                elif strand == '-':
                    astart = end - 1
                    aend = end + 13
                    dstart = start - 7
                    dend = start + 3
                aseq = ref.fetch(chrom, astart, aend).lower()
                dseq = ref.fetch(chrom, dstart, dend).lower()
    
                aseq_snp = None
                dseq_snp = None
                if -13 <= distance < 0:  # snp in intron region of acceptor site
                    assert aseq[(end_snp - astart - 1)]
                    aseq_snp = aseq[:(end_snp - astart - 1)] + alt_snp + aseq[(end_snp - astart):]
                    assert len(aseq) == len(aseq_snp), 'Unequal length: {0}, {1}'.format(len(aseq), len(aseq_snp))
                if 0 < distance <= 7:  # snp in intron region of donor site
                    assert dseq[(end_snp - dstart - 1)]
                    dseq_snp = dseq[:(end_snp - dstart - 1)] + alt_snp + dseq[(end_snp - dstart):]
                    assert len(dseq) == len(dseq_snp), 'Unequal length: {0}, {1}'.format(len(dseq), len(dseq_snp))
    
                if strand == '-':
                    aseq = self.reverse_complement(aseq)
                    dseq = self.reverse_complement(dseq)
                    aseq_snp = self.reverse_complement(aseq_snp)
                    dseq_snp = self.reverse_complement(dseq_snp)
    
                aic = 'NA'
                dic = 'NA'
                aic_change = 'NA'
                dic_change = 'NA'
                if distance < 0 and aseq.find('N') == -1:
                    aic = self.cal_acceptor_ic(aseq)
                if distance > 0 and dseq.find('N') == -1:
                    dic = self.cal_donor_ic(dseq)
                if aseq_snp:
                    if aseq_snp.find('N') == -1:
                        aic_change = self.cal_acceptor_ic(aseq_snp) - aic
                if dseq_snp:
                    if dseq_snp.find('N') == -1:
                        dic_change = self.cal_donor_ic(dseq_snp) - dic
                out_f.write('\t'.join(map(str, cols + [astart, aend, aseq, dstart, dend, dseq,
                                                       aic, dic, aic_change, dic_change])) + '\n')
            ref.close()
    
    def reverse_complement(self, seq):
        if seq:
            trans_tab = string.maketrans('ACGTacgt', 'TGCAtgca')
            seq_comp = seq.translate(trans_tab)
            seq_rev_comp = seq_comp[::-1]
            return seq_rev_comp
        return None

    def _build_header(self):
        header = ['#chrom_snp', 'start_snp', 'end_snp', 'ref', 'alt', 'feature', 'gene_id',
                  'chrom', 'start', 'end', 'name', 'score', 'strand', 'distance']
        header += ['acceptor_start',
                   'acceptor_end',
                   'acceptor_seq',
                   'donor_start',
                   'donor_end',
                   'donor_seq',
                   'aic',
                   'dic',
                   'aic_change',
                   'dic_change']
        return header


def main():
    parser = argparse.ArgumentParser(description='''
            Given AS strength matrix, bed file and reference fasta file, calculate
            the acceptor/donor site strength.''')
    parser.add_argument('ref_fname',
            help='reference fasta')
    parser.add_argument('bed_fname',
            help='bed file including skipped exon')
    parser.add_argument('ofname',
            help='output file name')
    parser.add_argument('--donor_ic')
    parser.add_argument('--acceptor_ic')
    args = parser.parse_args()
    junction_strength = JunctionStrength(args.donor_ic, args.acceptor_ic)
    junction_strength.cal_junction_strength(args.ref_fname, args.bed_fname, args.ofname)

if __name__ == '__main__':
    main()
