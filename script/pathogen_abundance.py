#!/usr/bin/env python
import math
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import FuncFormatter
import numpy as np
from config import *

def round_to_1(x):
   return round(x, -int(math.floor(math.log10(abs(x)))))


def plot(df, metadata, taxa, fig_num):
    print(taxa)
    all_samples = metadata[["subject", "day", "plot_group"]].drop_duplicates()
    df = df.query("taxonomy == @taxa")

    samples_with_data = df[["subject", "day"]].drop_duplicates()

    missing_samples = pd.merge(
       all_samples,
       samples_with_data,
       on=["subject", "day"],
       how="left",
       indicator=True
    ).query('_merge == "left_only"').drop(columns=["_merge"])

    missing_samples["pct"] = 0.0
    missing_samples["taxonomy"] = taxa
    missing_samples = missing_samples[df.columns]
    df = pd.concat([df, missing_samples], ignore_index=True)
    df = df.query("day > -5")

    grid_rows = 6
    grid_cols = 6
    grid_pos = 0
    fig = plt.figure(figsize=(18, 18), dpi=250)
    fig.subplots_adjust(hspace=0.4, wspace=0.1)
    sns.set(font_scale=0.9)

    for group in ["Placebo", "cohort1", "cohort2", "cohort3"]:
        for subject in sorted(df.query("plot_group == @group")["subject"].unique()):
            grid_pos += 1
            subj_df = df.query("subject == @subject").sort_values("day")
            subj_df["day"] = subj_df["day"].astype(str)

            max_pct = subj_df['pct'].max()
            subj_df['pct'] = subj_df['pct'] + 0.0001
            ax = fig.add_subplot(grid_rows, grid_cols, grid_pos)
            p = sns.lineplot(data=subj_df, x="day", y="pct", marker="o", ax=ax)

            ax.set_title(group_name.get(group, '?') + r": $\bf{"+str(subject)+"}$")
            ax.set_yscale("log")
            formatter = FuncFormatter(lambda y, _: round_to_1(y))
            ax.yaxis.set_major_formatter(formatter)

            ax.set_ylim(0.0001, 100)
            if grid_pos % grid_cols == 1:
                ax.set_ylabel(f"% {taxa.replace('_', ' ')}")
            else:
                ax.set_ylabel("")
                ax.set_yticklabels([])

            if grid_pos > grid_cols * (grid_rows-1):
                ax.set_xlabel("Day")
            else:
                ax.set_xlabel("")

            plt.xticks(rotation=45)



    out_fn = f"../figures/pathogens_{taxa.lower()}_supplementary_figure_S{fig_num}.pdf"
    plt.savefig(out_fn, bbox_inches="tight")


def read_mpa_taxa(fn, taxa, metric="frac"):

    if not os.path.exists(fn):
        return None
    with open(fn) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")

            last_tax_part = parts[0].split("|")[-1]
            if last_tax_part != taxa:
                continue

            if metric == "frac":
                return float(parts[2])
            else:
                return int(parts[4])

    return 0

def get_baseline_alpha(data):
    baseline_values = {}

    for row in data.to_dict(orient="records"):
        if row["Day"] in [-1, -2]:
            baseline_values[row["unblinded_subject_id"]] = row["shannon_alpha"]

    return baseline_values

if __name__ == "__main__":
    metadata = pd.read_table("../data/metagenomics_sequencing_metadata.tsv")[["unblinded_subject_id","Day", "plot_group"]].rename(columns={'unblinded_subject_id':'subject','Day':'day'})
    genera = [
       {"name":"Citrobacter", "fig_num": "10"},
       {"name":"Klebsiella", "fig_num": "8"},
       {"name":"Salmonella", "fig_num": "9"},
       {"name":"Staphylococcus", "fig_num": "13"},
    ]
    species = [
       {"name":"Clostridioides_difficile", "fig_num": "11"},
       {"name":"Clostridium_perfringens", "fig_num": "12"},
       {"name":"Enterobacter_cloacae_complex", "fig_num": "15"},
       {"name":"Streptococcus_agalactiae", "fig_num": "14"},
    ]

    df = pd.read_table("../data/genus_abundance_data.tsv")
    for data in genera:
       plot(df, metadata, data["name"], data["fig_num"])

    df = pd.read_table("../data/species_abundance_data.tsv")
    for data in species:
       plot(df, metadata, data["name"], data["fig_num"])


