#!/usr/bin/env python
import os.path
import sqlite3


class DB(object):
    def __init__(self, fname):
        self.fname = os.path.expanduser(fname)  # sqlite3 db file name

    def query_structure(self, transcript_id):
        conn = sqlite3.connect(self.fname)
        conn.row_factory = sqlite3.Row
        c = conn.cursor()
        transcript_id_tuple = (transcript_id, )
        c.execute('select * from structure where transcript_id == ?', transcript_id_tuple)
        r = c.fetchall()
        assert len(r) == 1 or r == [], '{0} contains multiple spineD records.\n'.format(transcript_id)
        structure_info = {}
        if r:
            structure_info['transcript_id'] = r[0]['transcript_id']
            structure_info['aa'] = r[0]['aa'].split(',')
            structure_info['beta_sheet'] = map(float, r[0]['beta_sheet'].split(','))
            structure_info['random_coil'] = map(float, r[0]['random_coil'].split(','))
            structure_info['alpha_helix'] = map(float, r[0]['alpha_helix'].split(','))
            structure_info['asa'] = map(float, r[0]['asa'].split(','))
            structure_info['disorder'] = map(float, r[0]['disorder'].split(','))
            structure_info['length'] = len(structure_info['aa'])
        conn.close()
        return structure_info

    def query_pfam(self, transcript_id):
        conn = sqlite3.connect(self.fname)
        conn.row_factory = sqlite3.Row
        c = conn.cursor()
        transcript_id_tuple = (transcript_id, )
        c.execute('select * from pfam where transcript_id == ?', transcript_id_tuple)
        r = c.fetchall()
        pfam_info = []
        if r:
            for record in r:
                pfam_info.append((record['transcript_id'],
                                  record['start'],
                                  record['end'],
                                  record['family'],
                                  record['name'],
                                  record['clan']))
        conn.close()
        return pfam_info

    def query_ptm(self, transcript_id):
        conn = sqlite3.connect(self.fname)
        conn.row_factory = sqlite3.Row
        c = conn.cursor()
        transcript_id_tuple = (transcript_id, )
        c.execute('select * from ptm where transcript_id == ?', transcript_id_tuple)
        r = c.fetchall()
        ptm_info = []
        if r:
            for record in r:
                ptm_info.append((record['transcript_id'],
                                 record['uniprot_id'],
                                 record['position'],
                                 record['modification']))
        conn.close()
        return ptm_info


def main():
    db = DB("../db/ensembl.db")
    structure = db.query_structure('ENST00000597436')
    print(structure)
    pfam = db.query_pfam('ENST00000409508')
    print(pfam)
    ptm = db.query_ptm('ENST00000409508')
    print(ptm)

if __name__ == '__main__':
    main()
