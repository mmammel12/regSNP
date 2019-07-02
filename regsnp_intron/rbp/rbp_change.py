#!/usr/bin/env python
import collections
import os.path
import numpy as np
import scipy.stats
import scipy.optimize

import rbp_score


class RBPChange(object):
    def __init__(self, pssm_path, pssm_list_fname, ms_fname):
        """
        RBP binding score change.
        :param pssm_path: the path of PSSM files
        :param pssm_list_fname: a list of valid PSSM files
        :param ms_fname: the file contains binding/nonbinding mean and sd for all PSSMs
        :return:
        """
        self.rbps = rbp_score.RBPScore(pssm_path, pssm_list_fname)
        self.mean_var = self.parse_mean_var(ms_fname)
        self.ics = self.cal_ic()
        self.beta_params = self.cal_beta_params()

    def parse_mean_var(self, fname):
        """
        Parse pre-computed mean and variance of matching score for binding and non-binding events.
        :param fname: contains five columns: pssm_fname, binding_mean, binding_sd, nonbinding_mean, nonbinding_sd
        :return: a collections.OrderedDict of tuple, each element is {pssm_fname: [binding_mean, binding_sd,
        nonbinding_mean, nonbinding_sd]}
        """
        mean_var = collections.OrderedDict()
        with open(os.path.expanduser(fname)) as f:
            header = f.readline()
            for line in f:
                motif, bmean, bsd, nbmean, nbsd = line.rstrip().split('\t')
                mean_var[motif] = (map(float, [bmean, bsd, nbmean, nbsd]))
        return mean_var

    def cal_ic(self, background=[0.25, 0.25, 0.25, 0.25]):
        ics = collections.OrderedDict()  # information contents for all the PSSMs
        for rbp in self.rbps.pssms:
            ics[rbp] = self._cal_ic(self.rbps.pssms[rbp].pssm, background)
        return ics

    def _cal_ic(self, pssm, background):
        """
        Calculate information content.
        :param pssm: 4 x m numpy ndarray of PSSM.
        :param background: background frequency
        :return: information content.
        """
        background = np.array(background)
        background = np.repeat(background, pssm.shape[1]).reshape(4,-1)  # Convert background to a matrix.
        freq = (2 ** pssm) * background
        ic = np.sum(freq * pssm)
        return ic

    def cal_beta_params(self):
        """
        Calculate beta distribution parameters a, b for all PSSMs.
        :return:
        """
        beta_params = collections.OrderedDict() # beta paramaters a, b for all the PSSMs
        for rbp in self.ics:
            mode = 1 / 2 ** self.ics[rbp]
            quantile = 0.005
            qtl_val = mode / 10.0
            beta_params[rbp] = self._cal_beta_params(mode, quantile, qtl_val)
        return beta_params

    def _cal_beta_params(self, mode, quantile, qtl_val, lower=1.1, upper=10):
        """
        Calculate beta distribution parameters a, b based on mode and quantile value atl_val.
        :param mode: set to 1/2 ** ic. For beta distribution, mode = (a - 1) / (a + b - 2).
        :param quantile:
        :param qtl_val:
        :param lower: lower end points of the interval to be searched for the root.
        :param upper: upper end points of the interval to be searched for the root.
        :return:
        """
        def beta_prior(x, mode, quantile, qtl_val):
            return scipy.stats.beta.cdf(qtl_val, x, (x - 1) / mode - (x - 1) + 1) - quantile
        a = scipy.optimize.brentq(beta_prior, a=lower, b=upper, args=(mode, quantile, qtl_val))
        b = (a - 1) / mode - (a - 1) + 1  # relationship between alpha, beta and mode
        return a, b

    def cal_change(self, score_fname, out_fname, keep_ncol=6):
        """
        Given matching score file, calculate log2(Odds Ratio) and posterior probability of binding change.
        :param score_fname: output of RBP_score.cal_matching_score()
        :param out_fname:
        :param keep_ncol:
        :return:
        """
        with open(os.path.expanduser(score_fname)) as score_f, open(os.path.expanduser(out_fname), 'w') as out_f:
            header = score_f.readline().rstrip().split('\t')
            header = header[:keep_ncol]
            header = header + [x + '_logOR' for x in self.rbps.pssm_list] + [x + '_pvalue' for x in self.rbps.pssm_list]
            out_f.write('\t'.join(header) + '\n')
            for line in score_f:
                cols = line.rstrip().split('\t')
                chrom, pos, ref, alt, ref_seq, alt_seq = cols[:keep_ncol]
                ref_scores = map(float, cols[keep_ncol:(keep_ncol + self.rbps.pssm_num)])
                alt_scores = map(float, cols[(keep_ncol + self.rbps.pssm_num):])
                assert len(ref_scores) == len(alt_scores)

                log_ors = []
                pvalues = []

                for i, rbp in enumerate(self.rbps.pssm_list):
                    bmean, bsd, nbmean, nbsd = self.mean_var[rbp]
                    a, b = self.beta_params[rbp]
                    pvalue = self._cal_prob(ref_scores[i], alt_scores[i], bmean, bsd, nbmean, nbsd, a, b)
                    pvalues.append(pvalue)
                    log_or = self._cal_logor(ref_scores[i], alt_scores[i], bmean, bsd, nbmean, nbsd)
                    log_ors.append(log_or)
                out_f.write('\t'.join(map(str, [chrom, pos, ref, alt, ref_seq, alt_seq] + log_ors + pvalues)) + '\n')

    def _cal_prob(self, ref_score, alt_score, bmean, bsd, nbmean, nbsd, a, b):
        """
        Calculate posterior probability of binding change when the prior follows a beta distribution specified by a, b.
        :param ref_score:
        :param alt_score:
        :param bmean:
        :param bsd:
        :param nbmean:
        :param nbsd:
        :param a:
        :param b:
        :return:
        """
        if max(ref_score, alt_score) < bmean - 3 * bsd:
            return 0
        else:
            b_ref = scipy.stats.norm.cdf(ref_score, bmean, bsd)
            b_alt = scipy.stats.norm.cdf(alt_score, bmean, bsd)
            nb_ref = 1 - scipy.stats.norm.cdf(ref_score, nbmean, nbsd)
            nb_alt = 1 - scipy.stats.norm.cdf(alt_score, nbmean, nbsd)

            def beta_distri(x):
                return x * (1 - x) * (b_ref * nb_alt + nb_ref * b_alt ) / \
                       ((x * b_ref + (1 - x) * nb_ref) * (x * b_alt + (1 - x) * nb_alt)) * \
                       scipy.stats.beta.pdf(x, a, b)

            prob, err = scipy.integrate.quad(beta_distri, 0, 1, epsrel=10 ** -10)
            return prob

    def _cal_logor(self, ref_score, alt_score, bmean, bsd, nbmean, nbsd):
        """
        Calculate log2(Odds Ratio) of binding change bewteen ref and alt.
        :param ref_score:
        :param alt_score:
        :param bmean:
        :param bsd:
        :param nbmean:
        :param nbsd:
        :return:
        """
        b_ref = scipy.stats.norm.cdf(ref_score, bmean, bsd)
        b_alt = scipy.stats.norm.cdf(alt_score, bmean, bsd)
        nb_ref = 1 - scipy.stats.norm.cdf(ref_score, nbmean, nbsd)
        nb_alt = 1 - scipy.stats.norm.cdf(alt_score, nbmean, nbsd)
        odds_ref = b_ref / nb_ref
        odds_alt = b_alt / nb_alt
        return np.log2(odds_alt / odds_ref)


def main():
    pass

if __name__ == '__main__':
    main()
