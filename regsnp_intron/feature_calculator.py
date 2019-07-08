#!/usr/bin/env python
import argparse
import json
import logging
import os
import os.path

import pandas as pd
import pymongo
from pymongo import MongoClient
import json
from bson.json_util import dumps
import operator
import csv

from utils.vcf import VCF
from utils.seq import FlankingSeq
from utils.snp import SNP
from utils.annovar import Annovar
from utils.closest_exon import ClosestExon
from rbp.rbp_change import RBPChange
from protein_feature.protein_feature import ProteinFeature
from junc_score.junction_strength import JunctionStrength
from conservation.phylop import Phylop


class FeatureCalculator(object):
    def __init__(self, settings, ifname, out_dir, iformat="txt"):
        self.settings = settings
        self.db_dir = os.path.expanduser(settings["db_dir"])
        self.ifname = os.path.expanduser(ifname)
        self.out_dir = os.path.expanduser(out_dir)
        self.iformat = iformat  # input format: txt or vcf
        self.logger = logging.getLogger(__name__)

    def calculate_feature(self):
        # return boolean
        needPredictor = False

        out_dir_tmp = os.path.join(self.out_dir, "tmp")
        os.mkdir(out_dir_tmp)

        # Convert input vcf to txt
        if self.iformat == "vcf":
            vcf_input = VCF(self.ifname)
            txt_input = os.path.join(out_dir_tmp, "snp_input.txt")
            vcf_input.convert_to_txt(txt_input)
            self.ifname = txt_input

        # Check input format
        self.logger.info("Checking input file format.")
        snp = SNP(self.ifname)
        ref_name = os.path.join(self.db_dir, "hg19/hg19.fa")
        if not snp.is_valid(ref_name):
            raise RuntimeError("Invalid input file format.")

        # Sort input
        self.logger.info("Sorting input file.")
        snp.sort(os.path.join(out_dir_tmp, "snp.sorted"))

        # Switch alleles
        self.logger.info("Switching alleles.")
        snp.switch_alleles(
            os.path.join(out_dir_tmp, "snp.sorted"),
            ref_name,
            os.path.join(out_dir_tmp, "snp.switched"),
        )

        # pull data from db
        needCompute = self._queryDB()

        # if self._queryDB returned true, continue with computation
        if needCompute:
            # Annotate SNVs
            self.logger.info("Annotating SNVs with ANNOVAR.")
            annovar_path = os.path.expanduser(self.settings["annovar_path"])
            annovar_db_path = os.path.expanduser(self.settings["annovar_db_path"])
            annovar = Annovar(annovar_path, annovar_db_path)
            annovar.annotate(
                os.path.join(out_dir_tmp, "snp.switched"),
                os.path.join(out_dir_tmp, "snp"),
            )

            # Calculate distance to closest protein coding exons and extract intronic SNVs
            self.logger.info("Calculating distance to closest protein coding exons.")
            protein_coding_exon_fname = os.path.join(
                self.db_dir, "hg19_ensGene_exon.bed"
            )
            closest_exon = ClosestExon(protein_coding_exon_fname)
            closest_exon.get_closest_exon(
                os.path.join(out_dir_tmp, "snp.intronic"),
                os.path.join(out_dir_tmp, "snp.distance"),
            )

            if os.path.getsize(os.path.join(out_dir_tmp, "snp.distance")) <= 0:
                self.logger.error("No input SNVs are in valid intronic regions.")
                raise RuntimeError("No input SNVs are in valid intronic regions.")

            # Extract flanking sequence
            self.logger.info("Fetching flanking sequence.")
            seq = FlankingSeq(ref_name, 20)
            seq.fetch_flanking_seq(
                os.path.join(out_dir_tmp, "snp.distance"),
                os.path.join(out_dir_tmp, "snp.seq"),
            )
            seq.fetch_flanking_seq(
                os.path.join(out_dir_tmp, "snp.distance"),
                os.path.join(out_dir_tmp, "snp.fa"),
                otype="fasta",
            )
            seq.close()

            # Calculate RBP binding change
            self.logger.info("Calculating RBP binding change.")
            pssm_path = os.path.join(self.db_dir, "motif/pwm")
            pssm_list_fname = os.path.join(self.db_dir, "motif/pwm_valid.txt")
            ms_fname = os.path.join(self.db_dir, "motif/binding_score_mean_sd.txt")
            rbp = RBPChange(pssm_path, pssm_list_fname, ms_fname)
            rbp.rbps.cal_matching_score(
                os.path.join(out_dir_tmp, "snp.seq"),
                os.path.join(out_dir_tmp, "snp.rbp_score"),
            )
            rbp.cal_change(
                os.path.join(out_dir_tmp, "snp.rbp_score"),
                os.path.join(out_dir_tmp, "snp.rbp_change"),
            )

            # Calculate protein structural features (disorder, secondary structure, ASA, Pfam, PTM)
            self.logger.info("Extracting protein structural features.")
            db_fname = os.path.join(self.db_dir, "ensembl.db")
            gene_pred_fname = os.path.join(self.db_dir, "hg19_ensGene.txt")
            protein_feature = ProteinFeature(db_fname, gene_pred_fname)
            protein_feature.calculate_protein_feature(
                os.path.join(out_dir_tmp, "snp.distance"),
                os.path.join(out_dir_tmp, "snp.protein_feature"),
            )

            # Calculate junction score
            self.logger.info("Calculating junction strength change.")
            donor_ic_fname = os.path.join(self.db_dir, "motif/donorsite.pssm")
            acceptor_ic_fname = os.path.join(self.db_dir, "motif/acceptorsite.pssm")
            junction_strength = JunctionStrength(donor_ic_fname, acceptor_ic_fname)
            junction_strength.cal_junction_strength(
                ref_name,
                os.path.join(out_dir_tmp, "snp.distance"),
                os.path.join(out_dir_tmp, "snp.junc"),
            )

            # Calculate phylop conservation score
            self.logger.info("Calculating conservation score.")
            phylop_fname = os.path.join(
                self.db_dir, "phylop/hg19.100way.phyloP100way.bw"
            )
            phylop = Phylop(phylop_fname)
            phylop.calculate(
                os.path.join(out_dir_tmp, "snp.distance"),
                os.path.join(out_dir_tmp, "snp.phylop"),
            )
            phylop.close()

            self._merge_features()
        # snp.switched doesn't exist
        else:
            # TODO change naming in _queryDB and uncomment this
            # rename prediction.txt to snp.prediction.txt
            # os.rename(
            #     os.path.join(self.out_dir, "prediction.txt"),
            #     os.path.join(out_dir_tmp, "snp.prediction.txt"),
            # )
            # rename prediction.json to snp.prediciton.json
            # os.rename(
            #     os.path.join(self.out_dir, "prediction.json"),
            #     os.path.join(out_dir_tmp, "snp.prediction.json"),
            # )
            # predictor not necessary, return false
            needPredictor = False
        return needPredictor

    def _queryDB(self):
        # return boolean
        # TODO change this to optionally allow computing
        needCalculate = False
        # tmp directory path
        out_dir_tmp = os.path.join(self.out_dir, "tmp")
        # create tempSwitched to rewrite snp.switched
        tempSwitched = ""
        # create output file string
        output = ""
        # create temp json dictionary
        json_str = "{ data: [ "
        # create connection to mongoD serverAdminB
        client = MongoClient(
            "mongodb+srv://cluster0-souoy.gcp.mongodb.net/test",
            username="serverAdmin",
            password="s3cr3tpass",
        )
        # get the DB
        db = client.muriDB
        # get the collection (basically a table)
        items = db.muriCol
        # parse snp.switched
        inFile = os.path.join(out_dir_tmp, "snp.switched")
        with open(inFile) as in_f:
            for line in in_f:
                cols = line.rstrip().split("\t")
                # build query dictionary
                query = {
                    "#chrom": cols[0],
                    "pos": cols[1],
                    "ref": cols[2],
                    "alt": cols[3],
                }
                # query for matching data
                item = items.find_one(query)
                # if data not in db
                if item == None:
                    # write line to tempSwitched, tab delimited
                    for element in cols:
                        tempSwitched += element + "\t"
                    # remove last \t from tempSwitched and replace with \n
                    tempSwitched = tempSwitched[:-2] + "\n"
                    # tempSwitched will be used to rewrite snp.switched
                # else, append all data to output file string, tab delimited
                else:
                    for key, value in item.iteritems():
                        if key != "_id":
                            output += value + "\t"
                    # remove last \t from line and replace with \n
                    output = output[:-2] + "\n"
                    # write data as JSON
                    json_str += dumps(item) + ","
        json_str += " ] }"
        json_str.encode("utf-8")
        # if tempSwitched is empty
        if len(tempSwitched) == 0:
            # return true
            needCalculate = False
        # else, overwrite snp.switched with tempSwitched
        else:
            with open(inFile, "w") as switched:
                switched.write(tempSwitched)
        # create file called prediction.txt and prediction.json in out_dir
        outFile = os.path.join(self.out_dir, "snp.prediction.txt")
        outJSONFile = os.path.join(self.out_dir, "snp.prediction.json")
        with open(outFile, "w") as out_f, open(outJSONFile, "w") as out_json_f:
            # write output file string to snp.prediction.txt and snp.prediction.json
            out_f.write(output)
            out_json_f.write(json_str)
        return needCalculate

    def _merge_features(self):
        self.logger.info("Merging all the features.")
        out_dir_tmp = os.path.join(self.out_dir, "tmp")
        rbp_change = pd.read_csv(
            os.path.join(out_dir_tmp, "snp.rbp_change"), sep="\t", header=0
        )
        snps = rbp_change.loc[:, ["#chrom", "pos", "ref", "alt"]]
        rbp_change.drop(
            ["#chrom", "pos", "ref", "alt", "ref_seq", "alt_seq"], axis=1, inplace=True
        )

        column_to_drop = [
            "#chrom_snp",
            "start_snp",
            "end_snp",
            "ref",
            "alt",
            "feature",
            "gene_id",
            "chrom",
            "start",
            "end",
            "score",
        ]
        protein_feature = pd.read_csv(
            os.path.join(out_dir_tmp, "snp.protein_feature"), sep="\t", header=0
        )
        protein_feature.drop(column_to_drop, axis=1, inplace=True)

        junction_strength = pd.read_csv(
            os.path.join(out_dir_tmp, "snp.junc"),
            sep="\t",
            header=0,
            usecols=["aic", "dic", "aic_change", "dic_change"],
        )

        phylop = pd.read_csv(
            os.path.join(out_dir_tmp, "snp.phylop"),
            sep="\t",
            header=0,
            usecols=["phylop1", "phylop3", "phylop7"],
        )

        result = pd.concat(
            [snps, rbp_change, protein_feature, junction_strength, phylop], axis=1
        )
        result.to_csv(
            os.path.join(self.out_dir, "snp.features.txt"),
            sep="\t",
            index=False,
            na_rep="NA",
        )

        pred = os.path.join(self.out_dir, "features.txt")
        features = os.path.join(self.out_dir, "snp.features.txt")
        # temp string to hold prediciton.txt data
        pred_string = ""
        # read in data of both files
        with open(pred) as dbFile, open(features) as snpFile:
            for line in dbFile:
                pred_string += line
            for line in snpFile:
                pred_string += line
        # write combined data to snp.features.txt
        with open(features, "w") as out_f:
            out_f.write(pred_string)
        data = csv.reader(features, delimiter="\t")
        # sort snp.features.txt
        data = pd.read_csv(features)
        data.sort_index(axis=0)
        data.to_csv(features, sep="\t", header=False, index=False)


def main():
    parser = argparse.ArgumentParser(
        description="""
            Given input SNP file, calculate features for classifier."""
    )
    parser.add_argument("sfname", help="JSON file containing settings.")
    parser.add_argument(
        "ifname", help="input SNP file. Contains four columns: chrom, pos, ref, alt."
    )
    parser.add_argument("out_dir", help="directory contains output files")
    args = parser.parse_args()

    settings = json.load(open(args.sfname))
    feature_calculator = FeatureCalculator(settings, args.ifname, args.out_dir)
    feature_calculator.calculate_feature()


if __name__ == "__main__":
    main()
