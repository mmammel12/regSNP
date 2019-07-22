#!/usr/bin/env python2.7
import argparse
import json
import logging
import os
import os.path

import pymongo
import json
import datetime
from bson.json_util import dumps
import csv

from utils.vcf import VCF
from utils.snp import SNP


class FeatureCalculator(object):
    def __init__(self, settings, ifname, out_dir, iformat="txt"):
        self.settings = settings
        self.db_URI = os.path.expanduser(settings["dbURI"])
        self.db_user = os.path.expanduser(settings["dbUsername"])
        self.db_pass = os.path.expanduser(settings["dbPassword"])
        self.ifname = os.path.expanduser(ifname)
        self.out_dir = os.path.expanduser(out_dir)
        self.hg_dir = os.path.expanduser(settings["hg_dir"])
        self.iformat = iformat  # input format: txt or vcf
        self.logger = logging.getLogger(__name__)

    def calculate_feature(self):
        out_dir_tmp = os.path.join(self.out_dir, "tmp")
        os.mkdir(out_dir_tmp)

        # Convert input vcf to txt
        if self.iformat == "vcf":
            vcf_input = VCF(self.ifname)
            txt_input = os.path.join(out_dir_tmp, "snp_input.txt")
            vcf_input.convert_to_txt(txt_input)
            self.ifname = txt_input

        # change delimiter to '\t' if necessary
        if self.iformat == "txt" or self.iformat == "csv":
            # determine delimiter
            with open(self.ifname) as in_f:
                delim = csv.Sniffer().sniff(in_f.read(1024))
                in_f.seek(0)
                # if delimiter is not '\t'
                if delim.delimiter != "\t":
                    # change delimiter to '\t'
                    reader = csv.reader(in_f, delimiter=delim.delimiter)
                    fixed_file = os.path.join(
                        self.out_dir, "snp_fixed_input.txt")
                    with open(fixed_file, "w") as out_f:
                        writer = csv.writer(out_f, delimiter="\t")
                        writer.writerows(reader)
                        # set fixed_file to self.ifname
                        self.ifname = fixed_file

        # Check input format
        self.logger.info("Checking input file format.")
        snp = SNP(self.ifname)

        # Sort input
        self.logger.info("Sorting input file.")
        snp.sort(os.path.join(out_dir_tmp, "snp.sorted"))

        # Switch alleles
        self.logger.info("Switching alleles.")
        snp.switch_alleles(
            os.path.join(out_dir_tmp, "snp.sorted"),
            self.hg_dir,
            os.path.join(out_dir_tmp, "snp.switched"),
        )

        # pull data from db
        self._queryDB()

    def _queryDB(self):
        self.logger.info("Querying database")
        # tmp directory path
        out_dir_tmp = os.path.join(self.out_dir, "tmp")
        # create headers string
        headers = ""
        # create temp json dictionary
        json_str = '{"data":['
        # create string to hold invalid lines
        invalid_str = ""
        # create connection to mongoDB
        client = pymongo.MongoClient(
            self.db_URI, username=self.db_user, password=self.db_pass
        )
        # get the DB
        db = client.muriDB
        # get the collection
        items = db.muriCol
        # get the queries collection
        queries = db.queries
        # create list to hold each query result
        resultsList = []
        # create counters for lines and errorLines
        lines = 0
        errorLines = 0
        # parse snp.switched
        inFile = os.path.join(out_dir_tmp, "snp.switched")
        with open(inFile) as in_f:
            needHeader = True
            for line in in_f:
                # increment lines
                lines += 1
                # get columns
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
                # add current date to query
                query["date"] = datetime.datetime.utcnow()
                # add query to queries collection
                queries.insert(query)
                # if data not in db
                if item == None:
                    # line isn't valid
                    # increment errorLines
                    errorLines += 1
                    # write to error log with css
                    errorMessage = "<div class='invalid'>Error: {0} {1} {2} {3} is not a valid combination. Line will not be included in results</div>".format(
                        cols[0], cols[1], cols[2], cols[3]
                    )
                    self.logger.info(errorMessage)
                    # write invalid data to invalid_str
                    invalid_str += "{0} {1} {2} {3}\n".format(
                        cols[0], cols[1], cols[2], cols[3]
                    )
                # else, append all data to output file string, tab delimited
                else:
                    # delete _id, it is not needed
                    del item["_id"]

                    # write header if needed
                    if needHeader:
                        headerList = [x for x in item.keys()]
                        resultsList.append(headerList)
                        needHeader = False
                    itemList = [x for x in item.values()]
                    resultsList.append(itemList)

                    # create dict with only necessary data for json
                    simple_json = {
                        "#chrom": item["#chrom"],
                        "pos": item["pos"],
                        "alt": item["alt"],
                        "ref": item["ref"],
                        "disease": item["disease"],
                        "splicing_site": item["splicing_site"],
                        "tpr": item["tpr"],
                        "fpr": item["fpr"],
                        "prob": item["prob"],
                        "name": item["name"],
                        "strand": item["strand"],
                    }
                    # fix strand type conversion
                    if simple_json["strand"] == "-0":
                        simple_json["strand"] = "-"
                    elif simple_json["strand"] == "0":
                        simple_json["strand"] = "+"
                    # write data as JSON
                    json_str += dumps(simple_json) + ","
        # if there is valid output
        if errorLines < lines:
            # end json_str
            json_str = json_str[:-1] + "]}"
            # remove unicode artifacts
            json_str.replace("u'", "'")
            # create file called snp.prediction.txt and snp.prediction.json and invalid.txt in out_dir
            outFile = os.path.join(self.out_dir, "snp.prediction.txt")
            outJSONFile = os.path.join(self.out_dir, "snp.prediction.json")
            invalidFile = os.path.join(self.out_dir, "invalid.txt")

            # create list of indices
            indices = []
            # sore resultsList by strand
            strandIndex = resultsList[0].index("strand")
            indices.append(strandIndex)
            resultsList[1:] = sorted(
                resultsList[1:], key=lambda x: x[strandIndex])
            # sore resultsList by name
            nameIndex = resultsList[0].index("name")
            indices.append(nameIndex)
            resultsList[1:] = sorted(
                resultsList[1:], key=lambda x: x[nameIndex])
            # sort resultsList by splicing_site
            ssIndex = resultsList[0].index("splicing_site")
            indices.append(ssIndex)
            resultsList[1:] = sorted(resultsList[1:], key=lambda x: x[ssIndex])
            # sore resultsList by fpr
            fprIndex = resultsList[0].index("fpr")
            indices.append(fprIndex)
            resultsList[1:] = sorted(
                resultsList[1:], key=lambda x: x[fprIndex])
            # sore resultsList by tpr
            tprIndex = resultsList[0].index("tpr")
            indices.append(tprIndex)
            resultsList[1:] = sorted(
                resultsList[1:], key=lambda x: x[tprIndex])
            # sore resultsList by prob
            probIndex = resultsList[0].index("prob")
            indices.append(probIndex)
            resultsList[1:] = sorted(
                resultsList[1:], key=lambda x: x[probIndex])
            # sore resultsList by disease
            diseaseIndex = resultsList[0].index("disease")
            indices.append(diseaseIndex)
            resultsList[1:] = sorted(
                resultsList[1:], key=lambda x: x[diseaseIndex])
            # sore resultsList by alt
            altIndex = resultsList[0].index("alt")
            indices.append(altIndex)
            resultsList[1:] = sorted(
                resultsList[1:], key=lambda x: x[altIndex])
            # sore resultsList by ref
            refIndex = resultsList[0].index("ref")
            indices.append(refIndex)
            resultsList[1:] = sorted(
                resultsList[1:], key=lambda x: x[refIndex])
            # sore resultsList by pos
            posIndex = resultsList[0].index("pos")
            indices.append(posIndex)
            resultsList[1:] = sorted(
                resultsList[1:], key=lambda x: x[posIndex])
            # sore resultsList by #chrom
            chromIndex = resultsList[0].index("#chrom")
            indices.append(chromIndex)
            resultsList[1:] = sorted(
                resultsList[1:], key=lambda x: x[chromIndex])

            # create indice order
            order = []
            for i in reversed(indices):
                # add index to order
                order.append(i)
                # add header string to headers
                headers += resultsList[0][i] + "\t"
            # add the rest of the headers after the sorted ones
            for i, value in enumerate(resultsList[0]):
                if i not in indices:
                    headers += value + "\t"
            # remove final \t and replace with \n
            headers = headers[:-2] + "\n"

            # add the rest of the indices to order
            for i in range(len(resultsList[0])):
                if i not in order:
                    order.append(i)

            if len(invalid_str) > 0:
                # invalid lines exits, write to file
                with open(invalidFile, "w") as invalid_f:
                    invalid_f.write(invalid_str)

            with open(outFile, "w") as out_f, open(outJSONFile, "w") as out_json_f:
                # write output file strings to snp.prediction.txt, snp.prediction.json
                out_f.write(headers)
                for i in resultsList[1:]:
                    for j in order:
                        # fix strand type conversion
                        if j != strandIndex:
                            out_f.write(i[j] + "\t")
                        else:
                            if i[j] == "-0":
                                out_f.write("-\t")
                            elif i[j] == "0":
                                out_f.write("+\t")
                    out_f.write("\n")
                out_json_f.write(json_str)


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
