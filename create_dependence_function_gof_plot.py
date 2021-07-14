import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from virocon import (get_OMAE2020_Hs_Tz, get_OMAE2020_V_Hs, 
    GlobalHierarchicalModel, plot_dependence_functions, read_ec_benchmark_dataset)

matplotlib.rcParams["svg.fonttype"] = "none" # Do avoid outputting font as path.
matplotlib.rcParams.update({"font.size": 6, "axes.labelsize": 8,})
matplotlib.rcParams["font.sans-serif"] = "Arial"
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["mathtext.fontset"] = "custom"
matplotlib.rcParams["mathtext.rm"] = "Arial"
matplotlib.rcParams["mathtext.it"] = "Arial:italic"
matplotlib.rcParams["mathtext.bf"] = "Arial:bold"

# Read dataset A, B  or C.
DATASET_CHARS = ["A", "B", "C"]

fig, axes = plt.subplots(2, 3, figsize=(8, 5), sharey="row")
for i, dataset_char in enumerate(DATASET_CHARS):
    file_path = "datasets/" + dataset_char + ".txt"
    sample = read_ec_benchmark_dataset(file_path)


    # Define the structure of the joint distribution model and fit it to the data.
    dist_descriptions, fit_descriptions, semantics = get_OMAE2020_Hs_Tz()
    model = GlobalHierarchicalModel(dist_descriptions)
    model.fit(sample)

    two_axes = [axes[0, i], axes[1, i]]
    par_rename = {"mu": "$\mu$", "sigma": "$\sigma$"}
    plot_dependence_functions(model, semantics, par_rename, two_axes)

    for j, ax in enumerate(two_axes):
        if j == 0:
            ax.set_title("Dataset $" + dataset_char + "$")
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.set_xlabel("")
        if i > 0:
            ax.set_ylabel("")
            
    
    # Thanks to: https://stackoverflow.com/questions/16150819/common-xlabel-ylabel-for-matplotlib-subplots
    fig.add_subplot(111, frameon=False)
    # Hide tick and tick label of the big axis.
    plt.tick_params(labelcolor="none", which="both", top=False, bottom=False, left=False, right=False)
    plt.xlabel("Significant wave height, $h_s$ (m)")
fig.tight_layout()
fig.savefig("figs/dependence-functions-abc.svg", dpi=300)

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
    #par_rename = {"alpha": "$\alpha$", "beta": "$\beta$"}
    plot_dependence_functions(model, semantics, par_rename, two_axes)

    for j, ax in enumerate(two_axes):
        if j == 0:
            ax.set_title("Dataset $" + dataset_char + "$")
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.set_xlabel("")
        if i > 0:
            ax.set_ylabel("")
            
    # Thanks to: https://stackoverflow.com/questions/16150819/common-xlabel-ylabel-for-matplotlib-subplots
    fig.add_subplot(111, frameon=False)
    # Hide tick and tick label of the big axis.
    plt.tick_params(labelcolor="none", which="both", top=False, bottom=False, left=False, right=False)
    plt.xlabel("Wind speed, $v$ (m s$^{-1}$)")
fig.tight_layout()
fig.savefig("figs/dependence-functions-def.svg", dpi=300)

plt.show()
