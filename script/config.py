def rgb(r, g, b):
    return (r/255, g/255, b/255)

cohort_colors = {
    "Placebo/Pre": rgb(235,144,4),
    "Placebo": rgb(235,144,4),
    "cohort1": rgb(0,62,104),
    "cohort2": rgb(4,122,127),
    "cohort3": rgb(160,226,185),
}


group_name = {
    'Placebo': 'Placebo',
    'cohort1': "SNIPR001 $10^{8}$",
    'cohort2': "SNIPR001 $10^{10}$",
    'cohort3': "SNIPR001 $10^{12}$",
}

phyla_colors = {
    'Firmicutes': '#1982c4',
    'Bacteroidetes': '#ff924c',
    'Actinobacteria': '#8ac926',
    'Proteobacteria': '#ff595e',
    'Verrucomicrobia': '#ffca3a',
    'Synergistetes': '#6a4c93',
    'Other': '#cccccc'
}

twenty_colors = [
    '#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#42d4f4',
    '#f032e6', '#bfef45', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8',
    '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#a9a9a9', '#000000'
]

settings = {
    's': {'data_file': '../data/species_abundance_data.tsv','hspace': 0.15, 'wspace': 1,    'font_scale': 0.8},
    'p': {'data_file': '../data/phyla_abundance_data.tsv',  'hspace': 0.3,  'wspace': 0.07, 'font_scale': 2},
}
