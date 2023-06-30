#!/usr/bin/env python3

import sys
import pandas as pd
from scipy.stats import pearsonr


def get_corr(df, t1, t2):
    corr, _ = pearsonr(df[t1], df[t2])
    corr = round(corr, 3)
    return corr


def dump_from(csv_path, l, label):
    fulldf = pd.read_csv(csv_path)
    fulldf = fulldf[
        fulldf["ESGq"].notna() & fulldf["rMATS"].notna() & fulldf["SUPPA2"].notna()
    ]
    df = fulldf[fulldf["k"] == 31]
    corr = get_corr(df, "ESGq", "rMATS")
    print(label, l, "\\esgq", "\\rmats", corr, sep=" & ", end=" \\\\\n")
    for k in [13, 21, 31]:
        df = fulldf[fulldf["k"] == k]
        corr = get_corr(df, "ESGq", "SUPPA2")
        print(label, l, "\\esgq", f"\\suppa (k{k})", corr, sep=" & ", end=" \\\\\n")
        corr = get_corr(df, "rMATS", "SUPPA2")
        print(label, l, "\\rmats", f"\\suppa (k{k})", corr, sep=" & ", end=" \\\\\n")


def main():
    csv_path_1 = sys.argv[1]  # 50, PE
    csv_path_2 = sys.argv[2]  # 100, PE
    csv_path_3 = sys.argv[3]  # 150, PE
    csv_path_4 = sys.argv[4]  # 150, SE

    dump_from(csv_path_1, 50, "PE")
    dump_from(csv_path_2, 100, "PE")
    dump_from(csv_path_3, 150, "PE")
    dump_from(csv_path_4, 150, "SE")


if __name__ == "__main__":
    main()
