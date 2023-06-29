import sys
import os
import glob

from Bio import SeqIO
import gffutils

NODE_IDX = 1


def open_gtf(gtf_path):
    try:
        gtf = gffutils.FeatureDB("{}.db".format(gtf_path), keep_order=True)
    except ValueError:
        gtf = gffutils.create_db(
            gtf_path,
            dbfn="{}.db".format(gtf_path),
            force=True,
            keep_order=True,
            disable_infer_genes=True,
            disable_infer_transcripts=True,
            merge_strategy="merge",
            sort_attribute_values=True,
        )
    return gtf


def print_gene_gaf(ref, eidx, Ts, ogfa, ogfa_add):
    global NODE_IDX

    nodes = {}
    for tidx, T in Ts.items():
        tidx = "P" if tidx == 1 else "S"
        for nidx, (s, e) in T:
            if nidx in nodes:
                continue
            nodes[nidx] = []
            nodeseq = ref[s - 1 : e - 1 + 1]
            subnodes = [nodeseq[i : i + 32] for i in range(0, len(nodeseq), 32)]
            for sn in subnodes:
                nodes[nidx].append(NODE_IDX)
                print("S", NODE_IDX, sn, sep="\t", file=ogfa)
                NODE_IDX += 1
        path = []
        last_nidx = -1
        for nidx in T:
            nidx = nidx[0]
            if last_nidx != -1:
                print(
                    "L",
                    nodes[last_nidx][-1],
                    "+",
                    nodes[nidx][0],
                    "+",
                    "*",
                    sep="\t",
                    file=ogfa,
                )
                print(eidx, tidx, nodes[last_nidx][-1], nodes[nidx][0], file=ogfa_add)
            for subnode1, subnode2 in zip(nodes[nidx][:-1], nodes[nidx][1:]):
                print("L", subnode1, "+", subnode2, "+", "*", sep="\t", file=ogfa)
                path.append(subnode1)
            path.append(nodes[nidx][-1])
            last_nidx = nidx
        print(
            "P",
            f"{eidx}_{tidx}",
            ",".join([f"{n}+" for n in path]),
            "*",
            sep="\t",
            file=ogfa,
        )


def analyze_ES(ref, exons, fpath, ogfa, ogfa_add):
    if not os.path.exists(fpath):
        return
    for line in open(fpath):
        if line.startswith("seqname"):
            continue
        idx = line.strip("\n").split("\t")[2]
        gene, rest = idx.split(";")
        _etype, chrom, intron1, intron2, strand = rest.split(":")
        intron1 = (int(intron1.split("-")[0]), int(intron1.split("-")[1]))
        intron2 = (int(intron2.split("-")[0]), int(intron2.split("-")[1]))
        idx = idx.replace(";", "_")
        idx = idx.replace(":", "_")
        pre, mid, post = (
            (float("inf"), intron1[0]),
            (intron1[1], intron2[0]),
            (intron2[1], -1),
        )
        for start, end in exons[chrom][gene]:
            if end == pre[1]:
                if start < pre[0]:
                    pre = (start, end)
            elif start == post[0]:
                if end > post[1]:
                    post = (start, end)

        Ts = {1: [(1, pre), (2, mid), (3, post)], 2: [(1, pre), (3, post)]}
        print_gene_gaf(ref[chrom], idx, Ts, ogfa, ogfa_add)


def analyze_SS(ref, exons, fpath, ogfa, ogfa_add):
    if not os.path.exists(fpath):
        return
    for line in open(fpath):
        if line.startswith("seqname"):
            continue
        idx = line.strip("\n").split("\t")[2]
        gene, rest = idx.split(";")
        _etype, chrom, intron1, intron2, strand = rest.split(":")
        intron1 = (int(intron1.split("-")[0]), int(intron1.split("-")[1]))
        intron2 = (int(intron2.split("-")[0]), int(intron2.split("-")[1]))
        idx = idx.replace(";", "_")
        idx = idx.replace(":", "_")
        mode = (strand == "+" and _etype == "A5") or (strand == "-" and _etype == "A3")

        if mode:
            alt1, alt2, const = (
                (float("inf"), intron1[0]),
                (float("inf"), intron2[0]),
                (intron1[1], -1),
            )
        else:
            alt1, alt2, const = (
                (intron1[1], -1),
                (intron2[1], -1),
                (float("inf"), intron1[0]),
            )
        for start, end in exons[chrom][gene]:
            if mode:
                if start == const[0]:
                    if end > const[1]:
                        const = (start, end)
                elif end == alt1[1]:
                    if start < alt1[0]:
                        alt1 = (start, end)
                elif end == alt2[1]:
                    if start < alt2[0]:
                        alt2 = (start, end)
            else:
                if end == const[1]:
                    if start < const[0]:
                        const = (start, end)
                elif start == alt1[0]:
                    if end > alt1[1]:
                        alt1 = (start, end)
                elif start == alt2[0]:
                    if end > alt2[1]:
                        alt2 = (start, end)
        Ts = {}
        begin, end = 0, 0
        if mode:
            begin = min(alt1[0], alt2[0])
            end = const[1]
            Ts = {1: [(1, alt1), (2, const)], 2: [(3, alt2), (2, const)]}
        else:
            begin = const[0]
            end = max(alt1[1], alt2[1])
            Ts = {1: [(1, const), (2, alt1)], 2: [(1, const), (3, alt2)]}

        print_gene_gaf(ref[chrom], idx, Ts, ogfa, ogfa_add)


def analyze_IR(ref, fpath, ogfa, ogfa_add):
    if not os.path.exists(fpath):
        return
    for line in open(fpath):
        if line.startswith("seqname"):
            continue
        idx = line.strip("\n").split("\t")[2]
        gene, rest = idx.split(";")
        _etype, chrom, exon1_s, intron, exon2_e, strand = rest.split(":")
        exon1 = (int(exon1_s), int(intron.split("-")[0]))
        exon2 = (int(intron.split("-")[1]), int(exon2_e))
        idx = idx.replace(";", "_")
        idx = idx.replace(":", "_")
        # if idx != "FBgn0013323_RI_2L_576081_576404-576514_576550_-":
        #     continue
        begin = exon1[0]
        end = exon2[1]
        Ts = {
            1: [(1, exon1), (2, exon2)],
            2: [(1, exon1), (3, (exon1[1] + 1, exon2[0] - 1)), (2, exon2)],
        }
        print_gene_gaf(ref[chrom], idx, Ts, ogfa, ogfa_add)


def build_gfa(fa_path, gtf_path, suppa2_wd, gfa_path, addinfo_path):
    ref = {}
    for record in SeqIO.parse(fa_path, "fasta"):
        ref[record.id] = str(record.seq)

    gtf = open_gtf(gtf_path)
    # print("Extracting exons..", file=sys.stderr)
    exons = {}
    for exon in gtf.features_of_type("exon"):
        chrom = exon.chrom
        gidx = exon.attributes["gene_id"][0]
        if chrom not in exons:
            exons[chrom] = {}
        if gidx not in exons[chrom]:
            exons[chrom][gidx] = set()
        exons[chrom][gidx].add((exon.start, exon.end))

    ogfa = open(gfa_path, "w")
    ogfa_add = open(addinfo_path, "w")
    # print("Analyzing ES..", file=sys.stderr)
    analyze_ES(ref, exons, os.path.join(suppa2_wd, "_SE_strict.ioe"), ogfa, ogfa_add)
    # print("Analyzing A3..", file=sys.stderr)
    analyze_SS(ref, exons, os.path.join(suppa2_wd, "_A3_strict.ioe"), ogfa, ogfa_add)
    # print("Analyzing A5..", file=sys.stderr)
    analyze_SS(ref, exons, os.path.join(suppa2_wd, "_A5_strict.ioe"), ogfa, ogfa_add)
    # print("Analyzing IR..", file=sys.stderr)
    analyze_IR(ref, os.path.join(suppa2_wd, "_RI_strict.ioe"), ogfa, ogfa_add)
    # print("Analyzing MX..", file=sys.stderr)
    # analyze_MX(exons, os.path.join(indir, "_MX_strict.ioe"))
    ogfa.close()
    ogfa_add.close()


if __name__ == "__main__":
    pass
