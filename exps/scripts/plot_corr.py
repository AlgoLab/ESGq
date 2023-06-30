import sys
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr

plt.rcParams.update({"font.size": 13})


def main():
    csv_path = sys.argv[1]
    prefix = sys.argv[2]

    fulldf = pd.read_csv(csv_path)
    fulldf = fulldf[
        fulldf["ESGq"].notna() & fulldf["rMATS"].notna() & fulldf["SUPPA2"].notna()
    ]
    for k in fulldf["k"].unique():
        print(k)
        df = fulldf[fulldf["k"] == k]
        print(k, len(df))
        corr, _ = pearsonr(df["ESGq"], df["rMATS"])
        corr = round(corr, 3)
        g = sns.jointplot(
            data=df,
            x="ESGq",
            y="rMATS",
            hue="EventType",
            kind="scatter",
            xlim=(-1.05, 1.05),
            ylim=(-1.05, 1.05),
        )
        plt.text(s=f"Pearson correlation: {corr}", x=-0.3, y=-1)
        plt.tight_layout()
        plt.savefig(prefix + f".k{k}.ESGq-vs-rMATS.png")
        plt.close()

        corr, _ = pearsonr(df["ESGq"], df["SUPPA2"])
        corr = round(corr, 3)
        sns.jointplot(
            data=df,
            x="ESGq",
            y="SUPPA2",
            hue="EventType",
            kind="scatter",
            xlim=(-1.05, 1.05),
            ylim=(-1.05, 1.05),
            legend=None,
        )
        plt.text(s=f"Pearson correlation: {corr}", x=-0.3, y=-1)
        plt.tight_layout()
        plt.savefig(prefix + f".k{k}.ESGq-vs-SUPPA2.png")
        plt.close()

        corr, _ = pearsonr(df["rMATS"], df["SUPPA2"])
        corr = round(corr, 3)
        sns.jointplot(
            data=df,
            x="rMATS",
            y="SUPPA2",
            hue="EventType",
            kind="scatter",
            xlim=(-1.05, 1.05),
            ylim=(-1.05, 1.05),
            legend=None,
        )
        plt.text(s=f"Pearson correlation: {corr}", x=-0.3, y=-1)
        plt.tight_layout()
        plt.savefig(prefix + f".k{k}.rMATS-vs-SUPPA2.png")


if __name__ == "__main__":
    main()
