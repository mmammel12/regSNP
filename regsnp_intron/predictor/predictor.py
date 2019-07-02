#!/usr/bin/env python
import argparse
import json
import logging
import os.path

import numpy as np
import pandas as pd
from sklearn.externals import joblib


class Predictor(object):
    def __init__(self, on_ss_imp, off_ss_imp, on_ss_clf, off_ss_clf, on_ss_fpr_tpr, off_ss_fpr_tpr):
        self.on_ss_imp = joblib.load(on_ss_imp)
        self.off_ss_imp = joblib.load(off_ss_imp)
        self.on_ss_clf = joblib.load(on_ss_clf)
        self.off_ss_clf = joblib.load(off_ss_clf)
        self.on_ss_fpr_tpr = pd.read_csv(on_ss_fpr_tpr, sep='\t')
        self.off_ss_fpr_tpr = pd.read_csv(off_ss_fpr_tpr, sep='\t')
        self.logger = logging.getLogger(__name__)

    def predict(self, ifname, ofname):
        """
        Predict diseasing-causing probability of all intronic SNVs
        :param ifname: m x n input file. m is the number of SNVs, n is the number of features.
        :param ofname: output file prefix.
        """
        self.logger.info('Loading all the features.')
        data = pd.read_csv(ifname, sep='\t', header=0)
        data_on_ss = data[(data['distance'] >= -13) & (data['distance'] <= 7) & (data['distance'] != 0)]
        data_off_ss = data[(data['distance'] < -13) | (data['distance'] > 7)]
        if not data_on_ss.empty:
            self.logger.info('Predicting on splicing site SNVs.')
            preds, scores = self._predict_on_ss(data_on_ss)
            fpr, tpr, disease = self._cal_fpr_tpr(scores, True)
            data_on_ss.insert(4, 'splicing_site', ['on'] * len(preds))
            data_on_ss.insert(4, 'fpr', fpr)
            data_on_ss.insert(4, 'tpr', tpr)
            data_on_ss.insert(4, 'prob', scores[:, 1])
            data_on_ss.insert(4, 'disease', disease)
        else:
            data_on_ss.insert(4, 'splicing_site', '')
            data_on_ss.insert(4, 'fpr', 0.0)
            data_on_ss.insert(4, 'tpr', 0.0)
            data_on_ss.insert(4, 'prob', 0.0)
            data_on_ss.insert(4, 'disease', '')
        if not data_off_ss.empty:
            self.logger.info('Predicting off splicing site SNVs')
            preds, scores = self._predict_off_ss(data_off_ss)
            fpr, tpr, disease = self._cal_fpr_tpr(scores, False)
            data_off_ss.insert(4, 'splicing_site', ['off'] * len(preds))
            data_off_ss.insert(4, 'fpr', fpr)
            data_off_ss.insert(4, 'tpr', tpr)
            data_off_ss.insert(4, 'prob', scores[:, 1])
            data_off_ss.insert(4, 'disease', disease)
        else:
            data_off_ss.insert(4, 'splicing_site', '')
            data_off_ss.insert(4, 'fpr', 0.0)
            data_off_ss.insert(4, 'tpr', 0.0)
            data_off_ss.insert(4, 'prob', 0.0)
            data_off_ss.insert(4, 'disease', '')
        self.logger.info('Generate final output.')
        result = pd.concat([data_on_ss, data_off_ss]).\
            sort_values(by=['#chrom', 'pos'])
        result.to_csv(ofname + '.txt', sep='\t', na_rep='NA', index=False)
        records = result.loc[:, ['#chrom', 'pos', 'ref', 'alt', 'disease', 'prob', 'tpr', 'fpr', 'splicing_site', 'name', 'strand']].\
            to_dict(orient='records')
        json.dump({"data": records}, open(ofname + '.json', 'w'))

    def _predict_on_ss(self, data):
        """
        Predict diseasing-causing probability of on-splicing-site SNVs
        :param data: m x n pandas data frame. m is the number of SNVs, n is the number of features.
        """
        data = data.drop(['#chrom', 'pos', 'ref', 'alt', 'name', 'strand', 'distance'], axis=1)
        data = self.on_ss_imp.transform(data)
        preds = self.on_ss_clf.predict(data)
        scores = self.on_ss_clf.predict_proba(data)
        return preds, scores

    def _predict_off_ss(self, data):
        """
        Predict diseasing-causing probability of off-splicing-site SNVs
        :param data: m x n pandas data frame. m is the number of SNVs, n is the number of features.
        """
        data = data.drop(['#chrom', 'pos', 'ref', 'alt', 'name', 'strand', 'distance', 'aic_change', 'dic_change'],
                         axis=1)
        data = self.off_ss_imp.transform(data)
        preds = self.off_ss_clf.predict(data)
        scores = self.off_ss_clf.predict_proba(data)
        return preds, scores

    def _cal_fpr_tpr(self, scores, on_ss):
        """

        :param scores: predicted scores
        :param on_ss: boolean indicate whether it's on_ss or off_ss
        :return: fpr, tpr and disease
        """
        fpr_tpr = None
        if on_ss:
            fpr_tpr = self.on_ss_fpr_tpr
        else:
            fpr_tpr = self.off_ss_fpr_tpr
        fpr = [fpr_tpr.fpr[(abs(x - fpr_tpr.pvalue)).idxmin()] for x in scores[:,1]]
        tpr = [fpr_tpr.tpr[(abs(x - fpr_tpr.pvalue)).idxmin()] for x in scores[:,1]]
        cutoff = [0.0, 0.05, 0.1, 1.1]
        disease_category = ['D', 'PD', 'B']
        disease = np.digitize(fpr, cutoff)
        disease = [disease_category[i - 1] for i in disease]
        return fpr, tpr, disease


def main():
    parser = argparse.ArgumentParser(description='''
            Given model files and SNV features generated by feature calculator,
            predict the disease-causing probability.''')
    parser.add_argument('db_path',
            help='path to models')
    parser.add_argument('sfname',
            help='SNV file with features generated by feature calculator')
    parser.add_argument('ofname',
            help='output file prefix')
    args = parser.parse_args()
    predictor = Predictor(os.path.join(args.db_path, 'on_ss_imp.pkl'),
                          os.path.join(args.db_path, 'off_ss_imp.pkl'),
                          os.path.join(args.db_path, 'on_ss_clf.pkl'),
                          os.path.join(args.db_path, 'off_ss_clf.pkl'),
                          os.path.join(args.db_path, 'on_ss_fpr_tpr.txt'),
                          os.path.join(args.db_path, 'off_ss_fpr_tpr.txt'))
    predictor.predict(args.sfname, args.ofname)

if __name__ == '__main__':
    main()
