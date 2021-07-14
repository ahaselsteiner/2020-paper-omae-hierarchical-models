import numpy as np
import matplotlib.pyplot as plt

from viroconcom.fitting import Fit
from viroconcom.contours import HighestDensityContour as HDC
from plot import plot_contour, PlottedSample, plot_marginal_fit, plot_dependence_functions
from contour_statistics import points_outside, sort_points_to_form_continous_line
from read_write import read_dataset, determine_file_name_e1, write_contour, read_contour

# Read dataset D, E  or F.
DATASET_CHAR = 'F'
file_path = 'datasets/' + DATASET_CHAR + '.txt'
sample_v, sample_hs, label_v, label_hs = read_dataset(file_path)
label_v = 'Wind speed (m s$^{-1}$)'
label_hs = 'Significant wave height (m)'

# Define the structure of the probabilistic model that will be fitted to the
# dataset.
dist_description_v = {'name': 'Weibull_Exp',
                      'dependency': (None, None, None, None),
                      'width_of_intervals': 2}
dist_description_hs = {'name': 'Weibull_Exp',
                       'fixed_parameters': (None, None, None, 5),
                       # shape, location, scale, shape2
                       'dependency': (0, None, 0, None),
                       # shape, location, scale, shape2
                       'functions': ('logistics4', None, 'alpha3', None),
                       # shape, location, scale, shape2
                       'min_datapoints_for_fit': 50,
                       'do_use_weights_for_dependence_function': True}


# Fit the model to the data.
fit = Fit((sample_v, sample_hs),
          (dist_description_v, dist_description_hs))

dist0 = fit.mul_var_dist.distributions[0]

fig = plt.figure(figsize=(12.5, 3.5), dpi=150)
ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)
plot_marginal_fit(sample_v, dist0, fig=fig, ax=ax1, label='$v$ (m s$^{-1}$)',
                  dataset_char=DATASET_CHAR)
plot_dependence_functions(fit=fit, fig=fig, ax1=ax2, ax2=ax3, unconditonal_variable_label=label_v)
fig.suptitle('Dataset ' + DATASET_CHAR)
fig.subplots_adjust(wspace=0.25, bottom=0.15)

# Compute highest density contours with return periods of 0.01, 1 and 50 years.
ts = 1 # Sea state duration in hours.
limits = [(0, 45), (0, 20)] # Limits of the computational domain.
deltas = [0.05, 0.05] # Dimensions of the grid cells.

return_period_lowest = 0.01
hdc_contour_lowest = HDC(fit.mul_var_dist, return_period_lowest, ts, limits, deltas)
return_period_1 = 1
hdc_contour_1 = HDC(fit.mul_var_dist, return_period_1, ts, limits, deltas)
return_period_50 = 50
hdc_contour_50 = HDC(fit.mul_var_dist, return_period_50, ts, limits, deltas)

c = sort_points_to_form_continous_line(hdc_contour_lowest.coordinates[0][0],
                                       hdc_contour_lowest.coordinates[0][1],
                                       do_search_for_optimal_start=True)
contour_v_lowest = c[0]
contour_hs_lowest = c[1]

contour_v_1 = hdc_contour_1.coordinates[0][0]
contour_hs_1 = hdc_contour_1.coordinates[0][1]
c = sort_points_to_form_continous_line(hdc_contour_1.coordinates[0][0],
                                       hdc_contour_1.coordinates[0][1],
                                       do_search_for_optimal_start=True)
contour_v_1 = c[0]
contour_hs_1 = c[1]

contour_v_50 = hdc_contour_50.coordinates[0][0]
contour_hs_50 = hdc_contour_50.coordinates[0][1]
c = sort_points_to_form_continous_line(hdc_contour_50.coordinates[0][0],
                                       hdc_contour_50.coordinates[0][1],
                                       do_search_for_optimal_start=True)
contour_v_50 = c[0]
contour_hs_50 = c[1]

fig = plt.figure(figsize=(5, 5), dpi=150)
ax = fig.add_subplot(111)

# Plot the two lower contours.
plot_contour(x=contour_v_lowest,
             y=contour_hs_lowest,
             ax=ax,
             line_style='r-')
plot_contour(x=contour_v_1,
             y=contour_hs_1,
             ax=ax,
             line_style='r-')

# Compute the median hs conditonal on v.
x = np.linspace(0, 35, 100)
d1 = fit.mul_var_dist.distributions[1]
a = d1.scale.a
b = d1.scale.b
c = d1.scale.c
y = a + b * np.power(x, c)

# Plot the 50-year contour and the sample.
plotted_sample = PlottedSample(x=np.asarray(sample_v),
                               y=np.asarray(sample_hs),
                               ax=ax)
plot_contour(x=contour_v_50,
             y=contour_hs_50,
             ax=ax,
             contour_label='Constant density',
             x_label=label_v,
             y_label=label_hs,
             line_style='r-',
             plotted_sample=plotted_sample,
             x_lim=(0, 35),
             upper_ylim=20,
             median_x=x,
             median_y=y,
             median_label='Median of $H_s | V$',
             median_style='b--')

plt.title('Dataset ' + DATASET_CHAR)
plt.show()
