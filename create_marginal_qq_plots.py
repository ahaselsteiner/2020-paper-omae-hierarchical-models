import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from virocon import (get_OMAE2020_Hs_Tz, get_OMAE2020_V_Hs, 
    GlobalHierarchicalModel, plot_marginal_quantiles, read_ec_benchmark_dataset)

matplotlib.rcParams["svg.fonttype"] = "none" # Do avoid outputting font as path.
matplotlib.rcParams.update({"font.size": 6})
matplotlib.rcParams["font.sans-serif"] = "Arial"
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["mathtext.fontset"] = "custom"
matplotlib.rcParams["mathtext.rm"] = "Arial"
matplotlib.rcParams["mathtext.it"] = "Arial:italic"
matplotlib.rcParams["mathtext.bf"] = "Arial:bold"

# Read dataset A, B  or C.
DATASET_CHARS = ["A", "B", "C"]

fig, axes = plt.subplots(2, 3, figsize=(7, 5), sharey="row")
for i, dataset_char in enumerate(DATASET_CHARS):
    file_path = "datasets/" + dataset_char + ".txt"
    sample = read_ec_benchmark_dataset(file_path)


    # Define the structure of the joint distribution model and fit it to the data.
    dist_descriptions, fit_descriptions, semantics = get_OMAE2020_Hs_Tz()
    model = GlobalHierarchicalModel(dist_descriptions)
    model.fit(sample)

    two_axes = [axes[0, i], axes[1, i]]
    plot_marginal_quantiles(model, sample, semantics, two_axes)

    for j, ax in enumerate(two_axes):
        if j == 0:
            ax.set_title("Dataset $" + dataset_char + "$")
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        if i > 0:
            ax.set_ylabel("")
        


fig.tight_layout()
fig.savefig("figs/marginals-qq-datasets-abc.svg", dpi=300)

# Read dataset D, E, F
DATASET_CHARS = ["D", "E", "F"]

fig, axes = plt.subplots(2, 3, figsize=(7, 5), sharey="row")
for i, dataset_char in enumerate(DATASET_CHARS):
    file_path = "datasets/" + dataset_char + ".txt"
    sample = read_ec_benchmark_dataset(file_path)


    # Define the structure of the joint distribution model and fit it to the data.
    dist_descriptions, fit_descriptions, semantics = get_OMAE2020_V_Hs()
    model = GlobalHierarchicalModel(dist_descriptions)
    model.fit(sample)

    two_axes = [axes[0, i], axes[1, i]]
    plot_marginal_quantiles(model, sample, semantics, two_axes)

    for j, ax in enumerate(two_axes):
        if j == 0:
            ax.set_title("Dataset $" + dataset_char + "$")
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        if i > 0:
            ax.set_ylabel("")
        


fig.tight_layout()
fig.savefig("figs/marginals-qq-datasets-def.svg", dpi=300)

plt.show()
