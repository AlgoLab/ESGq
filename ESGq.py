#!/usr/bin/env python3

import sys
import os
import argparse
import logging
import subprocess

from esg import build_gfa
from psi import compute_psi
from dpsi import compute_dpsi

FORMAT = "[%(asctime)s] %(message)s"
logging.basicConfig(stream=sys.stderr, format=FORMAT, level=logging.INFO)


def run_suppa(GTF, WD, events):
    log_path = os.path.join(WD, "suppa2.log")
    cmd = [
        "suppa.py",
        "generateEvents",
        "-i",
        GTF,
        "-f",
        "ioe",
        "-o",
        WD + "/",
        "-e",
    ] + events
    p = subprocess.run(cmd, stderr=open(log_path, "w"))
    return p.returncode, log_path


def run_vgindex(GFA, INDEX_PREFIX, THREADS):
    log_path = INDEX_PREFIX + ".autoindex.log"
    cmd = [
        "vg",
        "autoindex",
        "--workflow",
        "giraffe",
        "--gfa",
        GFA,
        "-t",
        str(THREADS),
        "-p",
        INDEX_PREFIX,
    ]
    p = subprocess.run(cmd, stderr=open(log_path, "w"))
    return p.returncode, log_path


def run_giraffe(INDEX_PREFIX, FQ1, FQ2, THREADS, GAF):
    log_path = GAF + ".vg-giraffe.log"
    if FQ2 == None:
        # single-end
        cmd = [
            "vg",
            "giraffe",
            "-Z",
            f"{INDEX_PREFIX}.giraffe.gbz",
            "-m",
            f"{INDEX_PREFIX}.min",
            "-d",
            f"{INDEX_PREFIX}.dist",
            "-f",
            FQ1,
            "-o",
            "GAF",
        ]
    else:
        # paired-end
        cmd = [
            "vg",
            "giraffe",
            "-Z",
            f"{INDEX_PREFIX}.giraffe.gbz",
            "-m",
            f"{INDEX_PREFIX}.min",
            "-d",
            f"{INDEX_PREFIX}.dist",
            "-f",
            FQ1,
            "-f",
            FQ2,
            "-o",
            "GAF",
        ]
    p = subprocess.run(cmd, stdout=open(GAF, "w"), stderr=open(log_path, "w"))
    return p.returncode, log_path


def main(args):
    # try:
    os.makedirs(args.WD, exist_ok=True)
    # except FileExistsError:
    #     logging.critical("Output folder already exits.")
    #     logging.critical("Halting..\n")
    #     sys.exit(1)

    logging.info("Running SUPPA2 to generate events..")
    suppa2_wd = os.path.join(args.WD, "suppa2")
    os.makedirs(suppa2_wd)
    retcode, suppa2log_path = run_suppa(args.GTF, suppa2_wd, ["SE", "SS", "RI"])
    if retcode != 0:
        logging.critical(f"SUPPA2 did not run succesfully (return code {retcode}).")
        logging.critical(f"See {suppa2log_path} for more details.")
        logging.critical("Halting..\n")
        sys.exit(1)
    logging.info("SUPPA2 ran succesfully.")

    logging.info("Builgind ESGs..")
    graph_prefix = os.path.join(args.WD, "ESGs")
    gfa_path = graph_prefix + ".gfa"
    addinfo_path = graph_prefix + ".addinfo"
    build_gfa(args.FA, args.GTF, suppa2_wd, gfa_path, addinfo_path)
    logging.info("Done.")

    logging.info("Indexing ESGs..")
    retcode, vgindexlog_path = run_vgindex(gfa_path, graph_prefix, args.threads)
    if retcode != 0:
        logging.critical(f"Indexing did not run succesfully (return code {retcode}).")
        logging.critical(f"See {vgindexlog_path} for more details.")
        logging.critical("Halting..\n")
        sys.exit(1)
    logging.info("Indexing ran succesfully.")

    # TODO: if paired, all must be paired
    C1 = {}
    for sample in args.C1:
        fq1_path = sample
        fq2_path = None
        if "," in sample:
            fq1_path, fq2_path = sample.split(",")

        # CHECKME: this may not work always
        bn = os.path.basename(fq1_path)
        if bn.endswith("gz"):
            bn = os.path.splitext(os.path.splitext(bn)[0])[0]
        else:
            bn = os.path.splitext(bn)[0]
        if fq2_path != None:
            bn = bn[:-2]

        logging.info(f"Aligning {bn}..")
        gaf_path = os.path.join(args.WD, f"{bn}.gaf")
        retcode, giraffelog_path = run_giraffe(
            graph_prefix, fq1_path, fq2_path, args.threads, gaf_path
        )
        if retcode != 0:
            logging.critical(
                f"Alignment did not run succesfully (return code {retcode})."
            )
            logging.critical(f"See {giraffelog_path} for more details.")
            logging.critical("Halting..\n")
            sys.exit(1)
        logging.info(f"Computing PSI from {bn}..")
        psi_path = gaf_path + ".psi"
        compute_psi(gaf_path, addinfo_path, psi_path)
        C1[bn] = psi_path
    C2 = {}
    for sample in args.C2:
        fq1_path = sample
        fq2_path = None
        if "," in sample:
            fq1_path, fq2_path = sample.split(",")
        bn = os.path.basename(fq1_path)
        if bn.endswith("gz"):
            bn = os.path.splitext(os.path.splitext(bn)[0])[0]
        else:
            bn = os.path.splitext(bn)[0]
        if fq2_path != None:
            bn = bn[:-2]

        logging.info(f"Aligning {bn}..")
        gaf_path = os.path.join(args.WD, f"{bn}.gaf")
        retcode, giraffelog_path = run_giraffe(
            graph_prefix, fq1_path, fq2_path, args.threads, gaf_path
        )
        if retcode != 0:
            logging.critical(
                f"Alignment did not run succesfully (return code {retcode})."
            )
            logging.critical(f"See {giraffelog_path} for more details.")
            logging.critical("Halting..\n")
            sys.exit(1)

        logging.info(f"Computing PSI from {bn}..")
        psi_path = gaf_path + ".psi"
        compute_psi(gaf_path, addinfo_path, psi_path)
        C2[bn] = psi_path

    dpsi_path = os.path.join(args.WD, "events.dpsi")
    logging.info(f"Computing dPSI..")
    compute_dpsi(C1, C2, dpsi_path)
    logging.info(f"dPSI stored in {dpsi_path}")
    logging.info("Done.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="ESGq",
        description="Event Splicing Graph -based Quantification of AS events across conditions",
    )
    parser.add_argument("FA", help="Reference in FASTA format")
    parser.add_argument("GTF", help="Gene annotation in GTF format")
    parser.add_argument("WD", help="Output/Working directory")
    parser.add_argument(
        "-1", dest="C1", nargs="+", help="Samples for condition 1 (can be gzipped)"
    )
    parser.add_argument(
        "-2", dest="C2", nargs="+", help="Samples for condition 2 (can be gzipped)"
    )
    parser.add_argument(
        "-t",
        dest="threads",
        default=1,
        type=int,
        help="Number of threads to use (default: 1)",
    )
    args = parser.parse_args()
    main(args)
