#!/usr/bin/env python
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from Bio import Phylo
from Bio.Phylo import BaseTree
from matplotlib.lines import Line2D

phylogroup_colors = {
    'A': "#1f77b4",
    'B1': "#ff7f0e",
    'B2': "#2ca02c",
    'D': "#d62728",
    'E': "#9467bd",
    'F': "#8c564b",
    'G': "#990022",
    'U': "#7f7f7f",
}


def plot_tree(tree, phylo_df):

    # find and remove the outgroup node
    outgroup_terminal = None
    for terminal in tree.get_terminals():
        if terminal.name == "E_albertii":
            tree.prune(terminal)
            break

    strain_to_phylogroup = dict(zip(phylo_df['strain_id'], phylo_df['Phylogroup']))

    # Setup the plot
    fig = plt.figure(figsize=(8,15))
    matplotlib.rc('font', size=7)
    ax_tree = fig.add_axes([0.1, 0.1, 0.6, 0.8])
    ax_boxes = fig.add_axes([0.675, 0.1, 0.05, 0.8])
    ax_legend = fig.add_axes([0.72, 0.6, 0.1, 0.3])

    # plot the actual tree
    Phylo.draw(tree, axes=ax_tree, do_show=False, show_confidence=False,
               branch_labels=None, label_func=lambda x: x.name if x.name else '')

    ax_tree.set_title('')
    ax_tree.spines[['top','right','left']].set_visible(False)
    ax_tree.get_yaxis().set_visible(False)


    # find the positions of all the terminals, to allow placing the annotation boxes in the right place
    leaf_names = [leaf.name for leaf in tree.get_terminals()]
    y_positions = {}
    for i, leaf_name in enumerate(leaf_names):
        for text_obj in ax_tree.texts:
            if text_obj.get_text() == ' '+leaf_name:
                y_positions[leaf_name] = text_obj.get_position()[1]
                break

    # draw annotations boxes for phylogroups
    ax_boxes.set_xlim(0, 1)
    ax_boxes.set_ylim(ax_tree.get_ylim())
    box_height = (ax_tree.get_ylim()[1] - ax_tree.get_ylim()[0]) / len(leaf_names)
    for leaf_name in leaf_names:
        ax_boxes.add_patch(
            patches.Rectangle(
                (0, y_positions[leaf_name] - box_height/2),
                0.5, box_height,
                linewidth=0,
                facecolor=phylogroup_colors[strain_to_phylogroup[leaf_name]]
            )
        )
    ax_boxes.spines[['top','right', 'bottom', 'left']].set_visible(False)
    ax_boxes.axis('off')


    # plot the legend
    phylogroup_colors['Unknown'] = phylogroup_colors['U']
    phylogroup_colors.pop("U")
    legend_elements = [Line2D([0], [0], marker='s', color='w',
                             markerfacecolor=phylogroup_colors[pg], markersize=10,
                             label=pg) for pg in sorted(phylogroup_colors)]

    ax_legend.legend(handles=legend_elements, loc='upper left', frameon=True, title='Phylogroups', title_fontsize=12, fontsize=10)
    ax_legend.set_xlim(0, 1)
    ax_legend.set_ylim(0, 1)
    ax_legend.axis('off')

    # Save it
    plt.savefig("../figures/annotated_tree_phylogroup_figure_S19.pdf", bbox_inches="tight")


if __name__ == "__main__":
    tree = Phylo.read("../data/strains.tree", 'newick')
    phylo_df = pd.read_table("../data/strains_phylogroups.tsv")
    plot_tree(tree, phylo_df)
