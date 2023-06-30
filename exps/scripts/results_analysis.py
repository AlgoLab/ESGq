#!/usr/bin/env python3

import sys
from matplotlib import pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
from pathlib import Path
import itertools
import pickle
import math


def parse_suppa(suppa_dpsi, pvalue=0.05):
    EVENTS = {"SE": {}, "A3": {}, "A5": {}, "RI": {}}
    for i, line in enumerate(open(suppa_dpsi)):
        if i == 0:
            continue
        idx, dpsi, pv = line.strip("\n").split("\t")
        dpsi, pv = float(dpsi), float(pv)
        if pv > pvalue:
            continue
        gene, rest = idx.split(";")
        etype, chrom, *positions, strand = rest.split(":")
        if not chrom.startswith("chr"):
            chrom = "chr" + chrom
        if etype == "SE":
            ab, cd = positions
            intron1 = tuple(int(x) for x in ab.split("-"))
            intron1 = (intron1[0] + 1, intron1[1] - 1)
            intron2 = tuple(int(x) for x in cd.split("-"))
            intron2 = (intron2[0] + 1, intron2[1] - 1)
            k = f"{chrom}:{intron1[0]}-{intron1[1]}-{intron2[0]}-{intron2[1]}"
            if k in EVENTS[etype]:
                print(f"Duplicate event at {k}", file=sys.stderr)
            EVENTS[etype][k] = dpsi
        elif (etype == "A5" and strand == "+") or (etype == "A3" and strand == "-"):
            ab, cd = positions
            shorter_intron = tuple(int(x) for x in ab.split("-"))
            shorter_intron = (shorter_intron[0] + 1, shorter_intron[1] - 1)
            longer_intron = tuple(int(x) for x in cd.split("-"))
            longer_intron = (longer_intron[0] + 1, longer_intron[1] - 1)
            assert longer_intron[1] == shorter_intron[1]
            k = f"{chrom}:{longer_intron[0]}-{shorter_intron[0]}-{longer_intron[1]}"
            if k in EVENTS[etype]:
                print(f"Duplicate event at {k}", file=sys.stderr)
            l = shorter_intron[0] - longer_intron[0] + 1
            EVENTS[etype][k] = (dpsi, l)
        elif (etype == "A3" and strand == "+") or (etype == "A5" and strand == "-"):
            ab, cd = positions
            shorter_intron = tuple(int(x) for x in ab.split("-"))
            shorter_intron = (shorter_intron[0] + 1, shorter_intron[1] - 1)
            longer_intron = tuple(int(x) for x in cd.split("-"))
            longer_intron = (longer_intron[0] + 1, longer_intron[1] - 1)
            assert longer_intron[0] == shorter_intron[0]
            k = f"{chrom}:{shorter_intron[0]}-{shorter_intron[1]}-{longer_intron[1]}"
            if k in EVENTS[etype]:
                print(f"Duplicate event at {k}", file=sys.stderr)
            l = longer_intron[1] - shorter_intron[1] + 1
            EVENTS[etype][k] = (dpsi, l)
        elif etype == "RI":
            a, bc, d = positions
            a = int(a)
            d = int(d)
            intron = tuple(int(x) for x in bc.split("-"))
            intron = (intron[0] + 1, intron[1] - 1)
            k = f"{chrom}:{a}-{intron[0]}-{intron[1]}-{d}"
            if k in EVENTS[etype]:
                print(f"Duplicate event at {k}", file=sys.stderr)
            EVENTS[etype][k] = dpsi
    return EVENTS


def parse_rmats_se(fpath, pvalue=0.05):
    EVENTS = {}
    for line in open(fpath):
        if line.startswith("ID"):
            continue
        (
            idx,
            gene,
            gene_sym,
            chrom,
            strand,
            ex_s,
            ex_e,
            usex_s,
            usex_e,
            dsex_s,
            dsex_e,
            _idx,
            inc_jc_1,
            sk_jc_1,
            inc_jc_2,
            sk_jc_2,
            inc_len,
            sk_len,
            pv,
            fdr,
            inclvl_1,
            inclvl_2,
            delta_incl,
        ) = line.strip("\n").split("\t")

        ex_s, usex_s, dsex_s = int(ex_s), int(usex_s), int(dsex_s)
        ex_e, usex_e, dsex_e = int(ex_e), int(usex_e), int(dsex_e)
        # Starts are 0-based, ends aren't (or they are exclusive)
        ex_s += 1
        usex_s += 1
        dsex_s += 1

        pv = float(pv)
        if pv > pvalue:
            continue

        intron1 = (int(usex_e) + 1, int(ex_s) - 1)
        intron2 = (int(ex_e) + 1, int(dsex_s) - 1)
        intron3 = (intron1[0], intron2[1])

        k = f"{chrom}:{intron1[0]}-{intron1[1]}-{intron2[0]}-{intron2[1]}"
        if k in EVENTS:
            print(f"Duplicate event at {k}", file=sys.stderr)
        EVENTS[k] = delta_incl
    return EVENTS


def parse_rmats_a3(fpath, pvalue=0.05):
    EVENTS = {}
    for line in open(fpath):
        if line.startswith("ID"):
            continue
        # Exon position are 0-based
        (
            idx,
            gene,
            gene_sym,
            chrom,
            strand,
            lex_s,
            lex_e,
            sex_s,
            sex_e,
            ex_s,
            ex_e,
            _idx,
            inc_jc_1,
            sk_jc_1,
            inc_jc_2,
            sk_jc_2,
            inc_len,
            sk_len,
            pv,
            fdr,
            inclvl_1,
            inclvl_2,
            delta_incl,
        ) = line.strip("\n").split("\t")

        # Starts are 0-based, ends aren't (or they are exclusive)
        ex_s, lex_s, sex_s = int(ex_s), int(lex_s), int(sex_s)
        ex_e, lex_e, sex_e = int(ex_e), int(lex_e), int(sex_e)
        ex_s += 1
        sex_s += 1
        sex_s += 1

        pv = float(pv)
        if pv > pvalue:
            continue

        # strand +
        longer_intron = (int(ex_e) + 1, int(sex_s) - 2)  # CHECKME: why do we need -2?
        shorter_intron = (int(ex_e) + 1, int(lex_s))
        if strand == "-":
            longer_intron = (int(sex_e) + 1, int(ex_s) - 1)
            shorter_intron = (int(lex_e) + 1, int(ex_s) - 1)

        k = ""
        l = 0
        if strand == "+":
            # assert longer_intron[0] == shorter_intron[0]
            k = f"{chrom}:{longer_intron[0]}-{shorter_intron[1]}-{longer_intron[1]}"
            l = longer_intron[1] - shorter_intron[1] + 1
        else:
            # assert longer_intron[1] == shorter_intron[1]
            k = f"{chrom}:{longer_intron[0]}-{shorter_intron[0]}-{longer_intron[1]}"
            l = shorter_intron[0] - longer_intron[0] + 1
        if k in EVENTS:
            print(f"Duplicate event at {k}", file=sys.stderr)
        EVENTS[k] = (delta_incl, l)
    return EVENTS


def parse_rmats_a5(fpath, pvalue=0.05):
    EVENTS = {}
    for line in open(fpath):
        if line.startswith("ID"):
            continue
        # Exon position are 0-based
        (
            idx,
            gene,
            gene_sym,
            chrom,
            strand,
            lex_s,
            lex_e,
            sex_s,
            sex_e,
            ex_s,
            ex_e,
            _idx,
            inc_jc_1,
            sk_jc_1,
            inc_jc_2,
            sk_jc_2,
            inc_len,
            sk_len,
            pv,
            fdr,
            inclvl_1,
            inclvl_2,
            delta_incl,
        ) = line.strip("\n").split("\t")

        # Starts are 0-based, ends aren't (or they are exclusive)
        ex_s, lex_s, sex_s = int(ex_s), int(lex_s), int(sex_s)
        ex_e, lex_e, sex_e = int(ex_e), int(lex_e), int(sex_e)
        ex_s += 1
        sex_s += 1
        sex_s += 1

        pv = float(pv)
        if pv > pvalue:
            continue

        # strand +
        longer_intron = (int(sex_e) + 1, int(ex_s) - 1)
        shorter_intron = (int(lex_e) + 1, int(ex_s) - 1)
        if strand == "-":
            longer_intron = (int(ex_e) + 1, int(sex_s) - 2)
            shorter_intron = (int(ex_e) + 1, int(lex_s))

        k = ""
        l = 0
        if strand == "+":
            k = f"{chrom}:{longer_intron[0]}-{shorter_intron[0]}-{longer_intron[1]}"
            l = shorter_intron[0] - longer_intron[0] + 1
        else:
            k = f"{chrom}:{longer_intron[0]}-{shorter_intron[1]}-{longer_intron[1]}"
            l = longer_intron[1] - shorter_intron[1] + 1
        if k in EVENTS:
            print(f"Duplicate event at {k}", file=sys.stderr)
        EVENTS[k] = (delta_incl, l)
    return EVENTS


def parse_rmats_ri(fpath, pvalue=0.05):
    EVENTS = {}
    for line in open(fpath):
        if line.startswith("ID"):
            continue
        # Exon position are 0-based
        (
            idx,
            gene,
            gene_sym,
            chrom,
            strand,
            ex_s,  # retained exon
            ex_e,
            fex_s,  # first exon
            fex_e,
            sex_s,  # second exon
            sex_e,
            _idx,
            inc_jc_1,
            sk_jc_1,
            inc_jc_2,
            sk_jc_2,
            inc_len,
            sk_len,
            pv,
            fdr,
            inclvl_1,
            inclvl_2,
            delta_incl,
        ) = line.strip("\n").split("\t")

        # Starts are 0-based, ends aren't (or they are exclusive)
        ex_s, fex_s, sex_s = int(ex_s), int(fex_s), int(sex_s)
        ex_e, fex_e, sex_e = int(ex_e), int(fex_e), int(sex_e)
        ex_s += 1
        sex_s += 1
        sex_s += 1

        fex_s += 1

        pv = float(pv)
        if pv > pvalue:
            continue

        assert ex_s == fex_s and ex_e == sex_e
        # strand +
        k = f"{chrom}:{ex_s}-{fex_e+1}-{sex_s-2}-{ex_e}"
        if k in EVENTS:
            print(f"Duplicate event at {k}", file=sys.stderr)
        EVENTS[k] = delta_incl
    return EVENTS


## esg dpsi parser
def parse_esg(esg_dpsi):
    EVENTS = {"SE": {}, "A3": {}, "A5": {}, "RI": {}}
    with open(esg_dpsi, "r") as f:
        for line in f:
            if line.startswith("Event"):
                continue
            idx, dpsi = (
                line.strip("\n").split("\t")[0],
                line.strip("\n").split("\t")[-1],
            )
            gene, etype, chrom, *positions, strand = idx.split("_")
            if not chrom.startswith("chr"):
                chrom = "chr" + chrom
            if etype == "SE":
                ab, cd = positions
                intron1 = tuple(int(x) for x in ab.split("-"))
                intron1 = (intron1[0] + 1, intron1[1] - 1)
                intron2 = tuple(int(x) for x in cd.split("-"))
                intron2 = (intron2[0] + 1, intron2[1] - 1)
                k = f"{chrom}:{intron1[0]}-{intron1[1]}-{intron2[0]}-{intron2[1]}"
                if k in EVENTS[etype]:
                    print(f"Duplicate event at {k}", file=sys.stderr)
                # print(k in suppa_events[etype].keys(), file=sys.stderr)
                EVENTS[etype][k] = dpsi
            elif (etype == "A5" and strand == "+") or (etype == "A3" and strand == "-"):
                ab, cd = positions
                shorter_intron = tuple(int(x) for x in ab.split("-"))
                shorter_intron = (shorter_intron[0] + 1, shorter_intron[1] - 1)
                longer_intron = tuple(int(x) for x in cd.split("-"))
                longer_intron = (longer_intron[0] + 1, longer_intron[1] - 1)
                assert longer_intron[1] == shorter_intron[1]
                k = f"{chrom}:{longer_intron[0]}-{shorter_intron[0]}-{longer_intron[1]}"
                if k in EVENTS[etype]:
                    print(f"Duplicate event at {k}", file=sys.stderr)
                l = shorter_intron[0] - longer_intron[0] + 1
                EVENTS[etype][k] = (dpsi, l)
            elif (etype == "A3" and strand == "+") or (etype == "A5" and strand == "-"):
                ab, cd = positions
                shorter_intron = tuple(int(x) for x in ab.split("-"))
                shorter_intron = (shorter_intron[0] + 1, shorter_intron[1] - 1)
                longer_intron = tuple(int(x) for x in cd.split("-"))
                longer_intron = (longer_intron[0] + 1, longer_intron[1] - 1)
                assert longer_intron[0] == shorter_intron[0]
                k = f"{chrom}:{shorter_intron[0]}-{shorter_intron[1]}-{longer_intron[1]}"
                if k in EVENTS[etype]:
                    print(f"Duplicate event at {k}", file=sys.stderr)
                l = longer_intron[1] - shorter_intron[1] + 1
                EVENTS[etype][k] = (dpsi, l)
            elif etype == "RI":
                a, bc, d = positions
                a = int(a)
                d = int(d)
                intron = tuple(int(x) for x in bc.split("-"))
                intron = (intron[0] + 1, intron[1] - 1)
                k = f"{chrom}:{a}-{intron[0]}-{intron[1]}-{d}"
                if k in EVENTS[etype]:
                    print(f"Duplicate event at {k}", file=sys.stderr)
                EVENTS[etype][k] = dpsi
    return EVENTS


def main():
    data_dir = sys.argv[1]
    pvalue = float(sys.argv[2])

    if data_dir[-1] == "/":
        data_dir = data_dir[:-1]

    k_kmer = ["13", "21", "31"]

    print("k", "EventType", "EventID", "rMATS", "SUPPA2", "ESGq", sep=",")
    
    rmats_prefix = f"{data_dir}/rMATS"
    esgq_dpsi = f"{data_dir}/ESGq/events.dpsi"

    rmats_se = parse_rmats_se(rmats_prefix + "/" + "SE.MATS.JC.txt", pvalue)
    rmats_a3 = parse_rmats_a3(rmats_prefix + "/" + "A3SS.MATS.JC.txt", pvalue)
    rmats_a5 = parse_rmats_a5(rmats_prefix + "/" + "A5SS.MATS.JC.txt", pvalue)
    rmats_ri = parse_rmats_ri(rmats_prefix + "/" + "RI.MATS.JC.txt", pvalue)
    esg_events = parse_esg(esgq_dpsi)

    for km in k_kmer:
        suppa_dpsi = f"{data_dir}/suppa2.k{km}/DIFF.dpsi"

        suppa_events = parse_suppa(suppa_dpsi, pvalue)
        etype = "SE"
        for k in set(rmats_se.keys()) | set(suppa_events[etype].keys()):
            chrom, positions = k.split(":")
            rpsi = float(rmats_se[k]) if k in rmats_se else float("NaN")
            spsi = (
                float(-suppa_events[etype][k])
                if k in suppa_events[etype]
                else float("NaN")
            )
            epsi = (
                float(esg_events[etype][k]) if k in esg_events[etype] else float("NaN")
            )
            print(km, etype, k, rpsi, spsi, epsi, sep=",")

        etype = "A3"
        for k in set(rmats_a3.keys()) | set(suppa_events[etype].keys()):
            chrom, positions = k.split(":")
            rpsi = float(rmats_a3[k][0]) if k in rmats_a3 else float("NaN")
            spsi = (
                float(-suppa_events[etype][k][0])
                if k in suppa_events[etype]
                else float("NaN")
            )
            epsi = (
                float(esg_events[etype][k][0])
                if k in esg_events[etype]
                else float("NaN")
            )
            print(km, etype, k, rpsi, spsi, epsi, sep=",")

        etype = "A5"
        for k in set(rmats_a5.keys()) | set(suppa_events[etype].keys()):
            chrom, positions = k.split(":")
            rpsi = float(rmats_a5[k][0]) if k in rmats_a5 else float("NaN")
            spsi = (
                float(-suppa_events[etype][k][0])
                if k in suppa_events[etype]
                else float("NaN")
            )
            epsi = (
                float(esg_events[etype][k][0])
                if k in esg_events[etype]
                else float("NaN")
            )
            print(km, etype, k, rpsi, spsi, epsi, sep=",")

        etype = "RI"
        for k in set(rmats_ri.keys()) | set(suppa_events[etype].keys()):
            chrom, positions = k.split(":")
            rpsi = float(rmats_ri[k]) if k in rmats_ri else float("NaN")
            spsi = (
                float(-suppa_events[etype][k])
                if k in suppa_events[etype]
                else float("NaN")
            )
            epsi = (
                1 - float(esg_events[etype][k]) if k in esg_events[etype] else float("NaN")
            )
            print(km, etype, k, rpsi, spsi, epsi, sep=",")


if __name__ == "__main__":
    main()
