#!/usr/bin/env python
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from config import *


# Since seaborn doesn't support stacked barplots, this implements a custom plotter for it
def stacked_barplot(ax, df, labels, colors, title, border=True, show_xlabel=True, show_ylabel=True, show_legend=False):

    sum_df = df.pivot_table(index="day", columns="taxonomy", values="pct", aggfunc="sum", fill_value=0)

    x_positions = range(len(sum_df.index))
    bottom_values = None

    for idx, column in enumerate(labels):
        if column not in sum_df:
            continue

        # Pick a color for the sub-bar
        if isinstance(colors, dict):
            col = colors[column]
        else:
            if column == "Other":
                col = "#dddddd"
            else :
                col = colors[idx]

        if border:
            ax.bar(x_positions, sum_df[column], label=column, bottom=bottom_values, color=col)
        else:
            ax.bar(x_positions, sum_df[column], label=column, bottom=bottom_values, color=col, linewidth=0)

        # Update bottom values for stacking
        bottom_values = sum_df[column] if bottom_values is None else bottom_values + sum_df[column]

    # Set the x-ticks and their labels
    ax.set_xticks(x_positions)
    xtick_labels = sum_df.index
    ax.set_xticklabels(xtick_labels, rotation=90)

    # Always show 0-100% on the y-axis
    ax.set_ylim(0,100)

    #if grid_pos > grid_cols*(grid_rows-1):
    if show_xlabel:
        ax.set_xlabel("Day")
    #if grid_pos % grid_cols == 1:
    if show_ylabel:
        ax.set_ylabel("Percentage")
    else:
        ax.set_yticklabels([])

    ax.set_title(title)

    if show_legend:
        ax.legend(title="Taxonomy", bbox_to_anchor=(1.03, 1), loc='upper left')


def plot_means(level, out_filename):
    assert level == "p"

    # Read abundance data
    df = pd.read_table(settings[level]['data_file'])

    # Initialize figure
    fig = plt.figure(figsize=(13, 3), dpi=250)
    fig.subplots_adjust(hspace=0.3, wspace=0.1)
    sns.set(rc={'axes.facecolor':'f4f4f4'}, font_scale=1.1)

    grid_rows, grid_cols, grid_pos = 1, 4, 0
    for group in df["plot_group"].unique():
        grid_pos += 1
        group_df = df.query("plot_group == @group")
        mean_df = group_df.groupby(['taxonomy', 'day'])['pct'].mean().reset_index()

        top_taxa = ["Actinobacteria","Bacteroidetes","Firmicutes","Proteobacteria", "Synergistetes", "Verrucomicrobia"]
        top_taxa_df = mean_df[mean_df["taxonomy"].isin(top_taxa)]

        # Sum up non-top species into Other
        other_taxa_df = (
            mean_df[~mean_df["taxonomy"].isin(top_taxa)]
            .groupby(["day"], as_index=False)
            .agg({"pct": "sum"})
            .assign(taxonomy="Other")
        )
        result_df = pd.concat([top_taxa_df, other_taxa_df], ignore_index=True)
        taxa = top_taxa + ["Other"]

        ax = fig.add_subplot(grid_rows, grid_cols, grid_pos)
        stacked_barplot(ax, result_df, taxa, phyla_colors, title = group_name[group], show_ylabel = grid_pos == 1)

    # Add the shared legend
    legend_handles = [mpatches.Patch(color=color, label=phylum) for phylum, color in phyla_colors.items()]
    fig.legend(
        handles=legend_handles,
        title="Phylum",
        loc='center left',
        bbox_to_anchor=(0.92, 0.65)
    )

    plt.savefig(out_filename, bbox_inches="tight")


def plot_individual(level, out_filename):

    assert level in ["s", "p"]

    # Read abundance data
    df = pd.read_table(settings[level]['data_file'])

    # Initialize figure
    fig = plt.figure(figsize=(43, 42), dpi=250)
    fig.subplots_adjust(hspace=settings[level]['hspace'], wspace=settings[level]['wspace'])
    sns.set(font_scale=settings[level]['font_scale'])

    grid_rows, grid_cols, grid_pos = 6, 6, 0
    for group in ["Placebo", "cohort1", "cohort2", "cohort3"]:
        for subject in df.query("plot_group == @group")["subject"].unique():
            grid_pos += 1

            subj_df = df.query("subject == @subject")
            subj_df.loc[:, 'taxonomy'] = subj_df['taxonomy'].str.replace('_', ' ', regex=False)

            # Select top taxonomies
            if level == "p":
                top_taxa = ["Actinobacteria","Bacteroidetes","Firmicutes","Proteobacteria", "Synergistetes", "Verrucomicrobia"]
            else:
                top_taxa = (subj_df.groupby("taxonomy")["pct"].mean().nlargest(20).index)
            top_taxa_df = subj_df[subj_df["taxonomy"].isin(top_taxa)]

            # Sum up the non-top taxonomies into "Other"
            other_taxa_df = (
                subj_df[~subj_df["taxonomy"].isin(top_taxa)]
                .groupby(["subject", "day", "plot_group"], as_index=False)
                .agg({"pct": "sum"})
                .assign(taxonomy="Other")
            )
            result_df = pd.concat([top_taxa_df, other_taxa_df], ignore_index=True)
            taxa = list(top_taxa) + ["Other"]


            # Draw the stacked barplot
            ax = fig.add_subplot(grid_rows, grid_cols, grid_pos)
            title = group_name.get(group, '?') + r": $\bf{"+str(subject)+"}$"
            show_xlabel = grid_pos > grid_cols*(grid_rows-1)
            show_ylabel = grid_pos % grid_cols == 1
            if level == "p":
                stacked_barplot(ax, result_df, taxa, phyla_colors, title, show_xlabel=show_xlabel, show_ylabel=show_ylabel)
            else:
                stacked_barplot(ax, result_df, taxa, twenty_colors, title, border=False, show_xlabel=show_xlabel, show_ylabel=show_ylabel, show_legend=True)

    # Draw the shared legend for phylum plot
    if level == "p":
        legend_handles = [mpatches.Patch(color=color, label=phylum) for phylum, color in phyla_colors.items()]
        fig.legend(
            handles=legend_handles,
            title="Phylum",
            loc='center left',
            bbox_to_anchor=(0.91, 0.835)
        )

    plt.savefig(out_filename, bbox_inches="tight")


if __name__ == "__main__":
    plot_means("p", "../figures/mean_phylum_composition_supplementary_figure_S5.pdf")
    plot_individual("s", "../figures/individual_species_composition_supplementary_figure_S7.pdf")
    plot_individual("p", "../figures/individual_phylum_composition_supplementary_figure_S6.pdf")


