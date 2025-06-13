#!/usr/bin/env python
import os
import math
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle

from config import *


if __name__ == "__main__":

    amr_df_raw = pd.read_table("../data/strain_amr.tsv")

    amr_df = (amr_df_raw.groupby(["unblinded_subject_id", "plot_group", "Day", "Drug"])["Susceptibility"].value_counts().unstack(fill_value=0).reset_index())
    amr_df["fraction_resistant"] = (
        (amr_df.get("RESISTANT", 0) + 0.5 * amr_df.get("INTERMEDIATE", 0)) /
        (amr_df.get("RESISTANT", 0) + amr_df.get("INTERMEDIATE", 0) + amr_df.get("SUSCEPTIBLE", 0))
    )
    df = amr_df.drop(columns=["INTERMEDIATE","RESISTANT","SUSCEPTIBLE"])

    drugs = df['Drug'].unique()
    plot_groups = df['plot_group'].unique()

    cols = 9
    rows = 4
    fig, axes = plt.subplots(
        rows,
        cols,
        figsize=(16, 9),
        constrained_layout=True,
        gridspec_kw = {'height_ratios': [2,1,1,2]},
    )

    cmap = LinearSegmentedColormap.from_list("white_to_red", ["#FFFFFF", "#A74337"])
    for j, group in enumerate(plot_groups):
        for i, drug in enumerate(drugs):
            ax = axes[j, i]
            subset = df.query("plot_group == @group and Drug == @drug").sort_values(["plot_group","unblinded_subject_id"])
            heatmap_data = subset.pivot(index="unblinded_subject_id", columns="Day", values=subset.columns[-1])
            heatmap_data.loc[r"$\bf{Total}$"] = heatmap_data.mean(skipna=True)

            sns.heatmap(heatmap_data, ax=ax, cmap=cmap, cbar=True if i == len(drugs)-1 and j==0 else False, linewidths=0.3, linecolor='#aaa', vmin=0, vmax=1, cbar_kws={'label': 'Fraction of isolates with resistance'})
            num_rows, num_cols = heatmap_data.shape
            ax.add_patch(Rectangle((0, 0), num_cols, num_rows, fill=False, edgecolor='#777', linewidth=0.6))
            ax.set_facecolor("lightgray")

            ax.tick_params(left = i==0, bottom = j==len(plot_groups)-1)
            if j == 0:
                ax.set_title(drug, fontsize=10)

            if j == len(plot_groups) - 1:
                ax.set_xlabel("Day")
            else:
                ax.set_xlabel("")
                ax.set_xticklabels([])

            if i > 0:
                ax.set_yticklabels([])
                ax.set_ylabel("")
            else:
                ax.set_ylabel(group_name[group])

    plt.savefig("../figures/amr_heatmaps_figure_S20.pdf", bbox_inches="tight")
