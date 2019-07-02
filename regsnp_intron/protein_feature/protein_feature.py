#!/usr/bin/env python
import argparse
import os.path
import re

import fetch_db
from spined import SpineD
from ptm import PTM
from pfam import Pfam
from gene_pred_ext import GenePredExt


class ProteinFeature(object):
    def __init__(self, db_fname, gene_pred_fname):
        self.db = fetch_db.DB(os.path.expanduser(db_fname))
        self.spined = SpineD(self.db)
        self.ptm = PTM(self.db)
        self.pfam = Pfam(self.db)

        self.gene_pred = GenePredExt(os.path.expanduser(gene_pred_fname))

    def calculate_protein_feature(self, bfname, out_fname):
        with open(os.path.expanduser(bfname)) as bed_f, open(os.path.expanduser(out_fname), 'w') as out_f:
            header = self._build_header()
            out_f.write('\t'.join(header) + '\n')
            for line in bed_f:
                cols = line.rstrip().split('\t')
                estart = int(cols[8])
                eend = int(cols[9])
                transcript_id = re.search(r'(\w+)_exon', cols[10]).group(1)
                if transcript_id in self.gene_pred.transcripts:
                    pstart, pend = self.gene_pred.get_protein_coord(transcript_id, estart, eend)
                    disorder, ss, asa = self.spined.cal_spined(transcript_id, pstart, pend)
                    ptm = self.ptm.cal_ptm(transcript_id, pstart, pend)
                    pfam = self.pfam.cal_pfam(transcript_id, pstart, pend)
                out_f.write('\t'.join(map(str, cols + list(disorder) + list(ss) + list(asa) + [ptm, pfam])) + '\n')

    def _build_header(self):
        header = ['#chrom_snp', 'start_snp', 'end_snp', 'ref', 'alt', 'feature', 'gene_id',
                  'chrom', 'start', 'end', 'name', 'score', 'strand', 'distance']
        header += ['min_disorder',
                   'max_disorder',
                   'mean_disorder',
                   'mean_disorder_struct',
                   'mean_disorder_dis',
                   'switch_num',
                   'min_disorder_len',
                   'max_disorder_len',
                   'mean_disorder_len',
                   'min_struct_len',
                   'max_struct_len',
                   'mean_struct_len']
        header += ['min_ss',
                   'max_ss',
                   'mean_ss',
                   'min_alpha',
                   'max_alpha',
                   'mean_alpha',
                   'min_beta',
                   'max_beta',
                   'mean_beta',
                   'min_coil',
                   'max_coil',
                   'mean_coil']
        header += ['min_asa',
                   'max_asa',
                   'mean_asa']
        header += ['ptm']
        header += ['pfam']
        return header


def main():
    parser = argparse.ArgumentParser(description='''
            Parse GenePred table (Extended) and extract features.''')
    parser.add_argument('ensembldb',
            help='ensembl sqlite db file, containing protein structure info.')
    parser.add_argument('gfname',
            help='GenePred table (Extended) file name, from UCSC table browser.')
    parser.add_argument('bfname',
            help='bedtools closest region distances output')
    parser.add_argument('ofname',
            help='output file name')
    args = parser.parse_args()
    protein_feature = ProteinFeature(args.ensembldb, args.gfname)
    protein_feature.calculate_protein_feature(args.bfname, args.ofname)

if __name__ == '__main__':
    main()
