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
from virocon import get_OMAE2020_V_Hs, GlobalHierarchicalModel

# Read dataset D, E or F.
DATASET_CHARS = ["D", "E", "F"]
for dataset_char in DATASET_CHARS:
    file_path = "datasets/" + dataset_char + ".txt"
    v, hs, label_v, label_hs = read_dataset(file_path)
    label_v = "Wind speed (m s$^{-1}$)"
    label_hs = "Significant wave height (m)"
    every_x_point = 1
    hs = hs[0::every_x_point]
    v = v[0::every_x_point]

    # Create figure, start with scatter plot.
    matplotlib.rcParams["svg.fonttype"] = "none" # Do avoid outputting font as path.
    matplotlib.rcParams.update({"font.size": 6})
    matplotlib.rcParams["font.family"] = "sans-serif"
    matplotlib.rcParams["font.sans-serif"] = "Arial"
    matplotlib.rcParams["mathtext.fontset"] = "custom"
    matplotlib.rcParams["mathtext.rm"] = "Arial"
    matplotlib.rcParams["mathtext.it"] = "Arial:italic"
    matplotlib.rcParams["mathtext.bf"] = "Arial:bold"
    fig1, ax = plt.subplots(figsize=(2.7, 2.7))
    color = my_colors.mpl_colors[2]
    levels = np.array([0.0001, 0.01])
    h_step = 0.1
    v_step = 0.2
    hgrid, vgrid = np.mgrid[0:12:h_step, 0:30:v_step]
    positions = np.vstack([hgrid.ravel(), vgrid.ravel()])
    CS_empirical = []
    CS = []
    ms = 10
    h_scatter = ax.scatter(
        v, hs, color=[0.5, 0.5, 0.5], s=ms, label="Dataset $" + dataset_char + "$", rasterized=True
    )


    # Kernel density estimate.
    kernel = gaussian_kde(np.vstack((hs, v)))
    Z = np.reshape(kernel(positions).T, hgrid.shape)
    CS_empirical.append(
        ax.contour(
            vgrid,
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
    dist_descriptions, fit_descriptions, semantics = get_OMAE2020_V_Hs()
    model = GlobalHierarchicalModel(dist_descriptions)
    data = np.array([v, hs])
    data = data.T
    model.fit(data)


    f = np.empty_like(hgrid)
    for i in range(hgrid.shape[0]):
        for j in range(hgrid.shape[1]):
            f[i, j] = model.pdf([vgrid[i, j], hgrid[i, j]])

    CS.append(
        ax.contour(
            vgrid,
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

    ax.set_xlabel("Wind speed (m s$^{-1}$)")
    ax.set_ylabel("Significant wave height (m)")
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    max_hs = 18
    max_v = 35
    ax.set_xlim((0, max_v))
    ax.set_ylim((0, max_hs))
    ax.set_title("")
    fig1.tight_layout()
    fig1.savefig("constant-density-dataset-" + dataset_char + ".svg", dpi=300)

    fig2, ax = plt.subplots(figsize=(2.7, 2.7))
    h_scatter = ax.scatter(
        v, hs, color=[0.5, 0.5, 0.5], s=ms, label="Dataset $" + dataset_char + "$", rasterized=True
    )
    v_step = 1
    interval_centers = np.arange(0.5 * v_step, max(v), v_step)
    hs_medians_empirical = np.empty(np.shape(interval_centers))
    hs_medians_empirical[:] = np.nan
    hs_medians_model = np.empty(np.shape(interval_centers))
    hs_medians_model[:] = np.nan
    dist = model.distributions[1]
    dep_func = dist.conditional_parameters["alpha"]

    for i, center in enumerate(interval_centers):
        lower = center - 0.5 * v_step
        upper = center + 0.5 * v_step
        mask = (v > lower) & (v < upper)
        hs_medians_empirical[i] = np.median(hs[mask])
    h_empirical = ax.plot(interval_centers, hs_medians_empirical, "-k", label="Median $H_s|V$ of dataset")
    v_values = np.arange(0, max_v, v_step)
    a = dep_func.parameters["a"]
    b = dep_func.parameters["b"]
    c = dep_func.parameters["c"]
    hs_medians_model = a + b * np.power(v_values, c)
    h_model = ax.plot(v_values, hs_medians_model, "--b", label="Median $H_s|V$ of model")
    ax.set_xlabel("Wind speed (m s$^{-1}$)")
    ax.set_ylabel("Significant wave height (m)")
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.set_xlim((0, max_v))
    ax.set_ylim((0, max_hs))
    ax.legend(handles=[h_scatter, h_empirical[0], h_model[0]], loc=2)
    fig2.tight_layout()
    fig2.savefig("median-comparison-dataset-" + dataset_char + ".svg", dpi=300)
plt.show()
