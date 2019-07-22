#!/usr/bin/env python
import argparse
import json
import logging
import os.path
import shutil
import sys


from feature_calculator import FeatureCalculator


def main():
    default_setting_path = os.path.realpath(
        os.path.join(
            os.path.dirname(os.path.realpath(__file__)
                            ), "../regsnp_intron/settings/settings.json"
        )
    )
    parser = argparse.ArgumentParser(
        description="""
            Given a list of intronic SNVs, predict the disease-causing probability
            based on genomic and protein structural features."""
    )
    parser.add_argument(
        "-s",
        "--sfname",
        default=default_setting_path,
        help="JSON file containing settings. Default setting file locate at: {0}".format(
            default_setting_path
        ),
    )
    parser.add_argument(
        "-f", "--force", action="store_true", help="overwrite existing directory"
    )
    parser.add_argument(
        "--iformat", choices=["txt", "vcf", "csv"], help="input format type"
    )
    parser.add_argument(
        "ifname", help="input SNV file. Contains four columns: chrom, pos, ref, alt."
    )
    parser.add_argument("out_dir", help="directory contains output files")
    args = parser.parse_args()

    out_dir = os.path.expanduser(args.out_dir)
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    elif args.force:
        shutil.rmtree(out_dir, ignore_errors=True)
        os.mkdir(out_dir)
    else:
        sys.stderr.write(
            "Output directory exists: {0}. Use '-f/--force' to overwrite.\n".format(
                out_dir
            )
        )
        sys.exit(1)

    logging.basicConfig(
        format="%(asctime)s - %(levelname)s - %(name)s - %(message)s",
        level=logging.INFO,
    )
    settings = json.load(open(args.sfname))
    try:
        feature_calculator = FeatureCalculator(
            settings, args.ifname, args.out_dir, args.iformat
        )
        feature_calculator.calculate_feature()
    except RuntimeError:
        shutil.rmtree(out_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
