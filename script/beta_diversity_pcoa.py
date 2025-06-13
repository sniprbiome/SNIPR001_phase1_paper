#!/usr/bin/env python
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import skbio
from config import *

def put_text(string, ax):
    y1, y2 = ax.get_ylim()
    height = y2 - y1
    x1, x2 = ax.get_xlim()
    width = x2 - x1
    ax.text(x1 + width*0.06, y2 - height*0.12, string)

def main():

    # Read full species abundance matrix
    species_data_all = pd.read_table("../data/species_abundance_matrix.tsv", index_col=0)

    # Read metagenomics sample metadata
    metagenomics_metadata_df = pd.read_table("../data/metagenomics_sequencing_metadata.tsv")#, index_col=0)

    sns.set(font_scale=1.1)
    sns.set_style("white")
    fig = plt.figure(figsize=(22, 6), dpi=250)

    all_permanova = []
    for i, cohort in enumerate(["cohort1", "cohort2", "cohort3"]):
        cohort_colors["Active"] = cohort_colors[cohort]

        # Select relevant samples
        metadata = metagenomics_metadata_df.query(f"cohort == @cohort and Day < 14").copy()

        # Set the status of the sample
        metadata.loc[:, "Status"] = np.where((metadata["treatment_group"] != "Placebo") & (metadata["Day"] > 0), "Active", "Placebo/Pre")

        # Subsample the abundance data to the relevant samples
        sample_ids = list(metadata["unblinded_sample_id"])
        species_data = species_data_all.loc[sample_ids]

        # Calculate bray-curtis dissimilarity matrix for selected samples
        bray_curtis_distance_matrix = skbio.diversity.beta_diversity("braycurtis", species_data.to_numpy(), sample_ids)

        # PERMANOVA test
        metadata = metadata.set_index("unblinded_sample_id")
        permanova = skbio.stats.distance.permanova(bray_curtis_distance_matrix, metadata["Status"])
        permanova["cohort"] = cohort
        all_permanova.append(permanova)

        # Principal coordinates plot (2d)
        bray_curtis_pcoa = skbio.stats.ordination.pcoa(bray_curtis_distance_matrix)
        two_dim = bray_curtis_pcoa.samples[['PC1', 'PC2']]
        two_dim = two_dim.merge(metadata[["Status","unblinded_subject_id"]], right_on="unblinded_sample_id", left_index=True).sort_values("Status")
        ax = fig.add_subplot(1, 3, i+1)
        sns.scatterplot(two_dim, x='PC1', y='PC2', hue="Status", palette=cohort_colors, edgecolor="black", ax=ax, s=50)
        put_text(f"F={permanova['test statistic']:.2f}\np={permanova['p-value']:.2f}", ax)

        plt.title(group_name[cohort])

    plt.savefig("../figures/beta_diversity_pcoa_supplementary_figure_S4.pdf", bbox_inches="tight")

    print(pd.DataFrame(all_permanova))

if __name__ == "__main__":
    main()
