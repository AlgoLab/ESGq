#!/usr/bin/env python3

import sys


def compute_psi(gaf_path, addinfo_path, psi_path):
    edges = {}
    for line in open(gaf_path):
        line = line.strip("\n").split("\t")
        if line[2] == "*":
            continue
        path = line[5]
        assert ("<" in path and ">" not in path) or ("<" not in path and ">" in path)
        strand = path[0] == ">"  # True: +
        path = path.replace("<", ">")
        path = path.split(">")[1:]
        if not strand:
            path = path[::-1]
        path = [int(x) for x in path]
        for x, y in zip(path[:-1], path[1:]):
            edges[(x, y)] = edges[(x, y)] + 1 if (x, y) in edges else 1

    opsi = open(psi_path, "w")
    last_eidx = ""
    isoforms = [[], []]
    for line in open(addinfo_path):
        if line.startswith("An") or line.startswith("Ex"):
            # FIXME
            continue
        eidx, t, x, y = line.strip("\n").split(" ")
        x, y = int(x), int(y)
        if last_eidx != "" and eidx != last_eidx:
            # compute
            PSI = 0
            w1, w2, w3 = 0, 0, 0
            if "SE" in last_eidx:
                w1 = edges[isoforms[0][0]] if isoforms[0][0] in edges else 0
                w2 = edges[isoforms[0][1]] if isoforms[0][1] in edges else 0
                w3 = edges[isoforms[1][0]] if isoforms[1][0] in edges else 0
                if w1 + w2 + w3 == 0:
                    PSI = "NaN"
                else:
                    PSI = ((w1 + w2) / 2) / ((w1 + w2) / 2 + w3)
            elif "A3" in last_eidx or "A5" in last_eidx:
                w1 = edges[isoforms[0][0]] if isoforms[0][0] in edges else 0
                w2 = edges[isoforms[1][0]] if isoforms[1][0] in edges else 0
                if w1 + w2 == 0:
                    PSI = "NaN"
                else:
                    PSI = w1 / (w1 + w2)
            elif "RI" in last_eidx:
                w1 = edges[isoforms[0][0]] if isoforms[0][0] in edges else 0
                w2 = edges[isoforms[1][0]] if isoforms[1][0] in edges else 0
                w3 = edges[isoforms[1][1]] if isoforms[1][1] in edges else 0
                if w1 + w2 + w3 == 0:
                    PSI = -1
                else:
                    PSI = 1 - w1 / (w1 + (w2 + w3) / 2)
            assert PSI == "NaN" or (PSI >= -1 and PSI <= 1)
            print(last_eidx, w1, w2, w3, PSI, file=opsi)
            isoforms = [[], []]
        last_eidx = eidx
        isoforms[0 if t == "P" else 1].append((x, y))
    # last event
    PSI = 0
    w1, w2, w3 = 0, 0, 0
    if "SE" in last_eidx:
        w1 = edges[isoforms[0][0]] if isoforms[0][0] in edges else 0
        w2 = edges[isoforms[0][1]] if isoforms[0][1] in edges else 0
        w3 = edges[isoforms[1][0]] if isoforms[1][0] in edges else 0
        if w1 + w2 + w3 == 0:
            PSI = "NaN"
        else:
            PSI = ((w1 + w2) / 2) / ((w1 + w2) / 2 + w3)
    elif "A3" in last_eidx or "A5" in last_eidx:
        w1 = edges[isoforms[0][0]] if isoforms[0][0] in edges else 0
        w2 = edges[isoforms[1][0]] if isoforms[1][0] in edges else 0
        if w1 + w2 == 0:
            PSI = "NaN"
        else:
            PSI = w1 / (w1 + w2)
    elif "RI" in last_eidx:
        w1 = edges[isoforms[0][0]] if isoforms[0][0] in edges else 0
        w2 = edges[isoforms[1][0]] if isoforms[1][0] in edges else 0
        w3 = edges[isoforms[1][1]] if isoforms[1][1] in edges else 0
        if w1 + w2 + w3 == 0:
            PSI = -1
        else:
            PSI = 1 - w1 / (w1 + (w2 + w3) / 2)
    assert PSI == "NaN" or (PSI >= -1 and PSI <= 1)
    print(last_eidx, w1, w2, w3, PSI, file=opsi)
    opsi.close()


if __name__ == "__main__":
    main()
