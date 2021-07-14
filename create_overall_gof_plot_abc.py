import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from viroconcom.fitting import Fit
from viroconcom.contours import HighestDensityContour as HDC
from plot import plot_contour, PlottedSample, plot_marginal_fit, plot_dependence_functions
from contour_statistics import points_outside, sort_points_to_form_continous_line
from read_write import read_dataset


from datetime import datetime
import matplotlib
from palettable.colorbrewer.qualitative import Set2_3 as my_colors
from palettable.colorbrewer.sequential import BuPu_5 as my_seq_colors
from scipy.stats import gaussian_kde
from virocon import get_OMAE2020_Hs_Tz, GlobalHierarchicalModel

# Read dataset A, B  or C.
DATASET_CHARS = ["A", "B", "C"]
for dataset_char in DATASET_CHARS:
    file_path = "datasets/" + dataset_char + ".txt"
    hs, tz, label_hs, label_tz= read_dataset(file_path)
    label_hs = "Significant wave height (m)"
    label_tz = "Zero-up-crossing period (s)"
    every_x_point = 1
    hs = hs[0::every_x_point]
    tz = tz[0::every_x_point]

    # Create figure, start with scatter plot.
    matplotlib.rcParams["svg.fonttype"] = "none" # Do avoid outputting font as path.
    matplotlib.rcParams.update({"font.size": 6})
    matplotlib.rcParams["font.sans-serif"] = "Arial"
    matplotlib.rcParams["font.family"] = "sans-serif"
    matplotlib.rcParams["mathtext.fontset"] = "custom"
    matplotlib.rcParams["mathtext.rm"] = "Arial"
    matplotlib.rcParams["mathtext.it"] = "Arial:italic"
    matplotlib.rcParams["mathtext.bf"] = "Arial:bold"
    fig1, ax = plt.subplots(figsize=(2.7, 2.7))
    color = my_colors.mpl_colors[2]
    levels = np.array([0.001, 0.05])
    h_step = 0.1
    t_step = 0.2
    hgrid, tgrid = np.mgrid[0:12:h_step, 0:20:t_step]
    positions = np.vstack([hgrid.ravel(), tgrid.ravel()])
    CS_empirical = []
    CS = []
    ms = 10
    h_scatter = ax.scatter(
        tz, hs, color=[0.5, 0.5, 0.5], s=ms, label="Dataset $" + dataset_char + "$", rasterized=True
    )


    # Kernel density estimate.
    kernel = gaussian_kde(np.vstack((hs, tz)))
    Z = np.reshape(kernel(positions).T, hgrid.shape)
    CS_empirical.append(
        ax.contour(
            tgrid,
            hgrid,
            Z,
            levels=levels[-2:],
            linestyles="-",
            linewidths=1,
            colors="k",
            zorder=2,
        )
    )
    CS_empirical[-1].collections[0].set_label("KDE, constant density")


    # Define the structure of the joint distribution model and fit it to the data.
    dist_descriptions, fit_descriptions, semantics = get_OMAE2020_Hs_Tz()
    model = GlobalHierarchicalModel(dist_descriptions)
    data = np.array([hs, tz])
    data = data.T
    model.fit(data)


    f = np.empty_like(hgrid)
    for i in range(hgrid.shape[0]):
        for j in range(hgrid.shape[1]):
            f[i, j] = model.pdf([hgrid[i, j], tgrid[i, j]])

    CS.append(
        ax.contour(
            tgrid,
            hgrid,
            f,
            levels=levels,
            zorder=2,
            colors="b",
            linestyles="--",
            linewidths=1,
        )
    )
    CS[-1].collections[0].set_label("Model, constant density")

    ax.legend(
        handles=[h_scatter, CS_empirical[0].collections[0], CS[0].collections[0]],
        loc=2,
    )

    ax.set_xlabel("Zero-up-crossing period (s)")
    ax.set_ylabel("Significant wave height (m)")
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    max_hs = 15
    ax.set_xlim((0, 18))
    ax.set_ylim((0, max_hs))
    ax.set_title("")
    fig1.tight_layout()
    fig1.savefig("constant-density-dataset-" + dataset_char + ".svg", dpi=300)
    
    fig2, ax = plt.subplots(figsize=(2.7, 2.7))
    h_scatter = ax.scatter(
        tz, hs, color=[0.5, 0.5, 0.5], s=ms, label="Dataset $" + dataset_char + "$", rasterized=True
    )
    hs_step = 0.2
    interval_centers = np.arange(0.5 * hs_step, max(hs), hs_step)
    tz_medians_empirical = np.zeros(np.shape(interval_centers))
    tz_medians_empirical[:] = np.nan
    tz_medians_model = np.zeros(np.shape(interval_centers))
    tz_medians_model[:] = np.nan
    dist = model.distributions[1]
    dep_func = dist.conditional_parameters["mu"]
    
    for i, center in enumerate(interval_centers):
        lower = center - 0.5 * hs_step
        upper = center + 0.5 * hs_step
        mask = (hs > lower) & (hs < upper)
        tz_medians_empirical[i] = np.median(tz[mask])
    h_empirical = ax.plot(tz_medians_empirical, interval_centers, "-k", label="Median $T_z|H_s$ of dataset")
    hs_values = np.arange(0, max_hs, hs_step)
    a = dep_func.parameters["a"]
    b = dep_func.parameters["b"]
    tz_medians_model = a + b * np.sqrt(np.divide(hs_values, 9.81))
    h_model = ax.plot(tz_medians_model, hs_values, "--b", label="Median $T_z|H_s$ of model")
    ax.set_xlabel("Zero-up-crossing period (s)")
    ax.set_ylabel("Significant wave height (m)")
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.set_xlim((0, 18))
    ax.set_ylim((0, max_hs))
    ax.legend(handles=[h_scatter, h_empirical[0], h_model[0]], loc=2)
    fig2.tight_layout()
    fig2.savefig("median-comparison-dataset-" + dataset_char + ".svg", dpi=300)
plt.show()
