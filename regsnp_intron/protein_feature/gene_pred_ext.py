#!/usr/bin/env python
import collections
import logging
import os.path


class GenePredExt(object):
    def __init__(self, fname):
        self.fname = os.path.expanduser(fname)
        self.logger = logging.getLogger(__name__)
        self.transcripts = self.extract_transcripts()

    def parse_file(self):
        with open(self.fname) as f:
            for line in f:
                if line.find('#') == 0:
                    continue
                cols = line.rstrip().split('\t')
                yield cols

    def extract_transcripts(self):
        gene_pred = self.parse_file()
        transcripts = collections.defaultdict(dict)
        for query_bin, transcript_id, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount,\
            exonStarts, exonEnds, score, gene_id, cdsStartStat, cdsEndStat, exonFrames in gene_pred:
            if cdsStartStat != 'none' or cdsEndStat != 'none':  # if it's a protein coding gene.
                txStart, txEnd, cdsStart, cdsEnd, exonCount = map(int, [txStart, txEnd, cdsStart, cdsEnd, exonCount])
                exonStarts = map(int, exonStarts.rstrip(',').split(','))
                exonEnds = map(int, exonEnds.rstrip(',').split(','))
                exonFrames = map(int, exonFrames.rstrip(',').split(','))
                transcripts[transcript_id]['transcript_id'] = transcript_id
                transcripts[transcript_id]['gene_id'] = gene_id
                transcripts[transcript_id]['chrom'] = chrom
                transcripts[transcript_id]['strand'] = strand
                transcripts[transcript_id]['txStart'] = txStart
                transcripts[transcript_id]['txEnd'] = txEnd
                transcripts[transcript_id]['cdsStart'] = cdsStart
                transcripts[transcript_id]['cdsEnd'] = cdsEnd
                transcripts[transcript_id]['exonCount'] = exonCount
                transcripts[transcript_id]['exons'] = zip(exonStarts, exonEnds)
                transcripts[transcript_id]['length'] = sum(x[1] - x[0] for x in zip(exonStarts, exonEnds))
                transcripts[transcript_id]['exonFrames'] = exonFrames
        self.logger.debug('Total number of protein coding transcripts: {0}.\n'.format(len(transcripts)))
        return transcripts

    def get_protein_coord(self, transcript_id, gstart, gend):
        """
        :param transcript_id: transcript id.
        :param gstart: genomic start of exon, 0-based.
        :param gend: genomic end of exon, 1-based.
        :return: start and end of protein coordinates.
        """
        pstart = None
        pend = None
        if transcript_id in self.transcripts:
            transcript = self.transcripts[transcript_id]
            strand = transcript['strand']
            cdsStart = transcript['cdsStart']
            cdsEnd = transcript['cdsEnd']
            exons = transcript['exons']
    
            g_cds_len = 0  # the CDS length of query exon
            g_cds_start = cdsStart if gstart < cdsStart < gend else gstart
            g_cds_end = cdsEnd if gstart < cdsEnd < gend else gend
            g_cds_len = g_cds_end - g_cds_start
    
            offset = 0
            if gstart < cdsEnd and gend > cdsStart:
                if strand == '+':
                    for estart, eend in exons:
                        cstart = estart
                        cend = eend
                        if cdsStart > eend:
                            continue
                        if estart <= cdsStart < eend:
                            cstart = cdsStart
                        if gstart != estart or gend != eend:
                            offset += cend - cstart
                        else:
                            break
                    pstart = offset
                    pend = offset + g_cds_len
                if strand == '-':
                    for estart, eend in reversed(exons):
                        if cdsEnd < estart:
                            continue
                        cstart = estart
                        cend = eend
                        if estart < cdsEnd <= eend:
                            cend = cdsEnd
                        if gstart != estart or gend != eend:
                            offset += cend - cstart
                        else:
                            break
                    pstart = offset
                    pend = offset + g_cds_len
                pstart /= 3
                pend /= 3
            else:
                self.logger.debug('Exon {0}-{1} is out of the CDS region {2}-{3} of transcript: {4}.\n'.
                                  format(gstart, gend, cdsStart, cdsEnd, transcript_id))
            return pstart, pend


if __name__ == '__main__':
    gene_pred = GenePredExt('../../db/hg19_ensGene.txt')
    print(gene_pred.get_protein_coord('ENST00000379236', 1147321, 1147518))
