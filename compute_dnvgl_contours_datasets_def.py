import numpy as np
import matplotlib.pyplot as plt

from viroconcom.fitting import Fit
from viroconcom.contours import HighestDensityContour as HDC
from plot import plot_contour, PlottedSample, plot_marginal_fit, plot_dependence_functions
from contour_statistics import points_outside, sort_points_to_form_continous_line
from read_write import read_dataset, determine_file_name_e1, write_contour, read_contour

# Read dataset D, E  or F.
DATASET_CHAR = 'D'
file_path = 'datasets/' + DATASET_CHAR + '.txt'
sample_v, sample_hs, label_v, label_hs = read_dataset(file_path)
label_v = 'wind speed (m s$^{-1}$)'

# Define the structure of the probabilistic model that will be fitted to the
# dataset. We will use the model that is recommended in DNV-RP-C205 (2010) on
# page 38 and that is called 'conditonal modeling approach' (CMA).
dist_description_hs = {'name': 'Weibull_3p',
                      'dependency': (None, None, None),
                      'width_of_intervals': 0.5}
dist_description_v = {'name': 'Weibull_2p',
                      'dependency': (0,  None, 0), #Shape, Location, Scale
                      'functions': ('power3', None, 'power3'), #Shape, Location, Scale
                      'min_datapoints_for_fit': 10
                      }

# Fit the model to the data.
fit = Fit((sample_hs, sample_v), (dist_description_hs, dist_description_v))

dist0 = fit.mul_var_dist.distributions[0]


fig = plt.figure(figsize=(12.5, 3.5), dpi=150)
ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)
plot_marginal_fit(sample_hs, dist0, fig=fig, ax=ax1, label='$h_s$ (m)',
                  dataset_char=DATASET_CHAR)
plot_dependence_functions(fit=fit, fig=fig, ax1=ax2, ax2=ax3, unconditonal_variable_label=label_hs)
fig.suptitle('Dataset ' + DATASET_CHAR)
fig.subplots_adjust(wspace=0.25, bottom=0.15)

# Compute highest density contours with return periods of 0.01, 1 and 50 years.
ts = 1 # Sea state duration in hours.
limits = [(0.10239, 20), (0, 45)] # Limits of the computational domain.
deltas = [0.2, 0.2] # Dimensions of the grid cells.

return_period_lowest = 0.01
hdc_contour_lowest = HDC(fit.mul_var_dist, return_period_lowest, ts, limits, deltas)
return_period_1 = 1
hdc_contour_1 = HDC(fit.mul_var_dist, return_period_1, ts, limits, deltas)
return_period_50 = 50
hdc_contour_50 = HDC(fit.mul_var_dist, return_period_50, ts, limits, deltas)

c = sort_points_to_form_continous_line(hdc_contour_lowest.coordinates[0][0],
                                       hdc_contour_lowest.coordinates[0][1],
                                       do_search_for_optimal_start=True)
contour_v_lowest = c[1]
contour_hs_lowest = c[0]

contour_v_1 = hdc_contour_1.coordinates[0][1]
contour_hs_1 = hdc_contour_1.coordinates[0][0]
c = sort_points_to_form_continous_line(hdc_contour_1.coordinates[0][0],
                                       hdc_contour_1.coordinates[0][1],
                                       do_search_for_optimal_start=True)
contour_v_1 = c[1]
contour_hs_1 = c[0]

contour_v_50 = hdc_contour_50.coordinates[0][1]
contour_hs_50 = hdc_contour_50.coordinates[0][0]
c = sort_points_to_form_continous_line(hdc_contour_50.coordinates[0][0],
                                       hdc_contour_50.coordinates[0][1],
                                       do_search_for_optimal_start=True)
contour_v_50 = c[1]
contour_hs_50 = c[0]

fig = plt.figure(figsize=(5, 5), dpi=150)
ax = fig.add_subplot(111)

# Plot the two lower contours.
plot_contour(x=contour_v_lowest,
             y=contour_hs_lowest,
             ax=ax,
             contour_label=str(return_period_lowest) + '-yr contour',
             line_style='r-')
plot_contour(x=contour_v_1,
             y=contour_hs_1,
             ax=ax,
             contour_label=str(return_period_1) + '-yr contour',
             line_style='r-')

# Plot the 50-year contour and the sample.
plotted_sample = PlottedSample(x=np.asarray(sample_v),
                               y=np.asarray(sample_hs),
                               ax=ax,
                               return_period=return_period_50)

plot_contour(x=contour_v_50,
             y=contour_hs_50,
             ax=ax,
             contour_label=str(return_period_50) + '-yr contour',
             x_label=label_v,
             y_label=label_hs,
             line_style='b-',
             plotted_sample=plotted_sample,
             x_lim=(0, 35))

plt.title('Dataset ' + DATASET_CHAR)
plt.show()
