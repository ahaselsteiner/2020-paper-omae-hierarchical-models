import numpy as np
import matplotlib.pyplot as plt
from scipy import stats


def plot_sample(plotted_sample, ax=None, do_plot_rasterized=True):
    """
    Plots the sample of metocean data.

    Parameters
    ----------
    plotted_sample : PlottedSample,
        The sample that should be plotted and its meta information.
    """
    if ax is None:
        ax = plotted_sample.ax
    ps = plotted_sample
    x = ps.x
    y = ps.y
    if ps.x_inside is not None and ps.y_inside is not None:
        if ps.return_period:
            inside_label = 'inside contour'
            outside_label = 'outside contour'
        else:
            inside_label = 'inside contour'
            outside_label = 'outside contour'
        ax.scatter(ps.x_inside, ps.y_inside, s=11, alpha=0.5, c='k',
                      marker='o', label=inside_label, rasterized=do_plot_rasterized)
        ax.scatter(ps.x_outside, ps.y_outside, s=9, alpha=0.5, c='r',
                      marker='D', label=outside_label, rasterized=do_plot_rasterized)
    else:
        if ps.label:
            ax.scatter(x, y, s=40, alpha=0.5, c='k', marker='.',
                          label=ps.label, rasterized=do_plot_rasterized)
        else:
            ax.scatter(x, y, s=40, alpha=0.5, c='k', marker='.',
                       label='Observation', rasterized=do_plot_rasterized)
    x_extremes = np.empty((4,1,))
    x_extremes[:] = np.nan
    y_extremes = np.empty((4,1,))
    y_extremes[:] = np.nan
    if ps.do_plot_extreme[0]:
        index_min_x = np.argmin(x)
        x_extremes[0] = x[index_min_x]
        y_extremes[0] = y[index_min_x]
    if ps.do_plot_extreme[1]:
        index_max_x = np.argmax(x)
        x_extremes[1] = x[index_max_x]
        y_extremes[1] = y[index_max_x]
    if ps.do_plot_extreme[2]:
        index_min_y = np.argmin(y)
        x_extremes[2] = x[index_min_y]
        y_extremes[2] = y[index_min_y]
    if ps.do_plot_extreme[3]:
        index_max_y = np.argmax(y)
        x_extremes[3] = x[index_max_y]
        y_extremes[3] = y[index_max_y]
    if any(d == True for d in ps.do_plot_extreme):
        ax.scatter(x_extremes, y_extremes,  s=40, c='g', marker='*',
                   label='observed extreme')


def plot_marginal_fit(sample, rv, fig, ax=None, label=None, color_sample='k',
                      marker_sample='x', marker_size_sample=3, color_fit='b',
                      dataset_char='?',legend_fontsize=8):
    #Plot Q-Q plot.
    if ax is None:
        ax = fig.add_subplot(111)
    plt.sca(ax)
    stats.probplot(sample, dist=rv, plot=ax)
    ax.get_lines()[0].set_markerfacecolor(color_sample)
    ax.get_lines()[0].set_markeredgecolor(color_sample)
    ax.get_lines()[0].set_marker(marker_sample)
    ax.get_lines()[0].set_markersize(marker_size_sample)
    ax.get_lines()[1].set_color(color_fit)
    ax.title.set_text('')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    if max(sample) < 12:
        plt.xlim((0, 12))
        plt.ylim((0, 15.5))
    else:
        plt.xlim((0, 35))
        plt.ylim((0, 35))
    if rv.name == 'ExponentiatedWeibull':
        dist_description = 'Exponentiated Weibull\n' \
                           '($\\alpha$=' + str('%.3g' % rv.scale(0)) + ', ' \
                           '$\\beta$=' + str('%.3g' % rv.shape(0)) + ', ' \
                           '$\\delta$=' + str('%.3g' % rv.shape2(0)) +')'
    elif rv.name == 'Weibull':
        dist_description = 'Weibull, ' \
                           '$\\alpha$=' + str('%.3g' % rv.scale(0)) + ', ' \
                           '$\\beta$=' + str('%.3g' % rv.shape(0))
    plt.legend(['Dataset '+ dataset_char, dist_description], loc='upper left',
               frameon=False, prop={'size': legend_fontsize})
    plt.xlabel('theoretical quantiles, ' + str(label).lower())
    plt.ylabel('ordered values, ' + str(label).lower())


def plot_dependence_functions(
        fit, fig, ax1=None, ax2=None, unconditonal_variable_label=None,
        marker_discrete='o',
        markersize_discrete=5, markerfacecolor_discrete='lightgray',
        markeredgecolor_discrete='k', style_dependence_function='b-',
        legend_fontsize=8):
    """

    Parameters
    ----------
    fit : Fit
    fig : Figure
    unconditonal_variable_label : str
    style_discrete_values : str
        Style of the discrete values for the parameters that were derived from
        marginal distribution fitting.
    style_dependence_function : str
        Style of the fitted dependence function.
    """
    if ax1 is None:
        ax1 = fig.add_subplot(121)
    if ax2 is None:
        ax2 = fig.add_subplot(122)

    plt.sca(ax1)
    scale_at = fit.multiple_fit_inspection_data[1].scale_at
    x1 = np.linspace(0, max(scale_at)*1.1, 100)
    if fit.mul_var_dist.distributions[1].scale.func_name == 'power3':
        dp_function = '$' + str('%.3g' % fit.mul_var_dist.distributions[1].scale.a) + \
                      '+' + str('%.3g' % fit.mul_var_dist.distributions[1].scale.b) + \
                      '\cdot h_s^{' + str('%.3g' % fit.mul_var_dist.distributions[1].scale.c) + '}$'
    elif fit.mul_var_dist.distributions[1].scale.func_name == 'lnsquare2':
        dp_function = '$\ln(' + str('%.3g' % fit.mul_var_dist.distributions[1].scale.a) + \
                      '+' + str('%.3g' % fit.mul_var_dist.distributions[1].scale.b) + \
                      '\sqrt{h_s / g})$'
    elif fit.mul_var_dist.distributions[1].scale.func_name == 'alpha3':
        dp_function = '$(' + str('%.3g' % fit.mul_var_dist.distributions[1].scale.a) + \
                      '+' + str('%.3g' % fit.mul_var_dist.distributions[1].scale.b) + \
                      '\cdot v^{' + str('%.3g' % fit.mul_var_dist.distributions[1].scale.c) + \
                      '}) / 2.0445^{(1 / \\beta_{hs})}$'
    else:
        dp_function = str(fit.mul_var_dist.distributions[1].scale)

    if fit.mul_var_dist.distributions[1].name == 'Lognormal':
        plt.plot(scale_at, np.log(fit.multiple_fit_inspection_data[1].scale_value),
                 marker_discrete,
                 markersize=markersize_discrete,
                 markerfacecolor=markerfacecolor_discrete,
                 markeredgecolor=markeredgecolor_discrete,
                 label='from marginal distribution')
        plt.plot(x1, np.log(fit.mul_var_dist.distributions[1].scale(x1)),
                 style_dependence_function, label=dp_function)
        ylabel = '$μ_{tz}$'
        plt.xlim((0, 6))
        plt.ylim((0.9, 2.15))
    if fit.mul_var_dist.distributions[1].name == 'Weibull' or \
                    fit.mul_var_dist.distributions[1].name == 'ExponentiatedWeibull':
        plt.plot(scale_at, fit.multiple_fit_inspection_data[1].scale_value,
                 marker_discrete,
                 markersize=markersize_discrete,
                 markerfacecolor=markerfacecolor_discrete,
                 markeredgecolor=markeredgecolor_discrete,
                 label='from marginal distribution')
        plt.plot(x1, fit.mul_var_dist.distributions[1].scale(x1),
                 style_dependence_function, label=dp_function)
        ylabel = '$α_{hs}$'
        plt.xlim((0, 30))
        plt.ylim((0, 10))
    plt.xlabel(unconditonal_variable_label)
    plt.legend(frameon=False, prop={'size': legend_fontsize})
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    plt.ylabel(ylabel)

    plt.sca(ax2)
    shape_at = fit.multiple_fit_inspection_data[1].shape_at
    x1 = np.linspace(0, max(shape_at)*1.1, 100)
    plt.plot(shape_at, fit.multiple_fit_inspection_data[1].shape_value,
             marker_discrete,
             markersize=markersize_discrete,
             markerfacecolor=markerfacecolor_discrete,
             markeredgecolor=markeredgecolor_discrete,)
    plt.plot(x1, fit.mul_var_dist.distributions[1].shape(x1),
             style_dependence_function)
    plt.xlabel(unconditonal_variable_label)
    if fit.mul_var_dist.distributions[1].name == 'Lognormal':
        plt.xlim((0, 6))
        plt.ylim((0.065, 0.33))
        ylabel = '$σ_{tz}$'
        if fit.mul_var_dist.distributions[1].shape.func_name == 'exp3':
            dp_function = '$' + str('%.3g' % fit.mul_var_dist.distributions[1].shape.a) + \
                          '+' + str('%.3g' % fit.mul_var_dist.distributions[1].shape.b) + \
                          '\exp (' + str('%.3g' % fit.mul_var_dist.distributions[1].shape.c) + \
                          'h_s)$'
        elif fit.mul_var_dist.distributions[1].shape.func_name == 'powerdecrease3':
            dp_function = '$' + str('%.4f' % fit.mul_var_dist.distributions[1].shape.a) + \
                          '+ 1 / (h_s + ' + str('%.3g' % fit.mul_var_dist.distributions[1].shape.b) + \
                          ')^{' + str('%.3g' % fit.mul_var_dist.distributions[1].shape.c) + \
                          '}$'
        elif fit.mul_var_dist.distributions[1].shape.func_name == 'asymdecrease3':
            dp_function = '$' + str('%.4f' % fit.mul_var_dist.distributions[1].shape.a) + \
                          ' + ' + str('%.3g' % fit.mul_var_dist.distributions[1].shape.b) + \
                          ' / (1 + ' + str('%.3g' % fit.mul_var_dist.distributions[1].shape.c) + \
                          ' h_s )$'
    if fit.mul_var_dist.distributions[1].name == 'Weibull'  or \
                    fit.mul_var_dist.distributions[1].name == 'ExponentiatedWeibull':
        ylabel = '$β_{h_s}$'
        plt.xlim((0, 30))
        plt.ylim((0, 3.5))
        if fit.mul_var_dist.distributions[1].shape.func_name == 'power3':
            dp_function = '$' + str('%.3g' % fit.mul_var_dist.distributions[1].shape.a) + \
                          '+' + str('%.3g' % fit.mul_var_dist.distributions[1].shape.b) + \
                          '\cdot h_s^{' + str('%.3g' % fit.mul_var_dist.distributions[1].shape.c) + '}$'
        elif fit.mul_var_dist.distributions[1].shape.func_name == 'logistics4':
            # logistics4 uses np.abs(c), to display it nicer, abs(c) is shown.
            absOfC = np.abs(fit.mul_var_dist.distributions[1].shape.c)
            dp_function = '$' + str('%.3g' % fit.mul_var_dist.distributions[1].shape.a) + \
                          '+' + str('%.3g' % fit.mul_var_dist.distributions[1].shape.b) + \
                          '/ [1 + e^{-' + str('%.3g' % absOfC) + \
                          '(v - ' + str('%.3g' % fit.mul_var_dist.distributions[1].shape.d) + \
                          ')}]$'
        else:
            dp_function = str(fit.mul_var_dist.distributions[1].shape)
    plt.legend(['from marginal distribution', dp_function], frameon=False, prop={'size': legend_fontsize})
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    plt.ylabel(ylabel)


def plot_contour(x, y, ax, contour_label=None, x_label=None, y_label=None,
                 line_style='b-', alpha=1, plotted_sample=None, x_lim = None,
                 upper_ylim=None, median_x=None, median_y=None, median_style='r-',
                 median_label='median of x2|x1'):
    """
    Plots the environmental contour.

    The method expects the coordinates to be ordered by angle.

    Parameters
    ----------
    x : ndarray of doubles
        The contour's coordinates in the x-direction.
    y : ndarray of doubles
        The contour's coordiantes in the y-direction.
    ax : Axes
        Axes of the figure where the contour should be plotted.
    contour_label : str, optional
        The environmental contour's label that will be used in the legend.
    x_label : str, optional
        Label for the x-axis.
    y_label : str, optional
        Label for the y-axis.
    line_style : str, optional
        Matplotlib line style.
    alpha : float, optional
        Alpha value (transparency) for the contour's line.
    plotted_sample : PlottedSample,
        The sample that should be plotted and its meta information.
    """
    # For generating a closed contour: add the first coordinate at the end.
    xplot = x.tolist()
    xplot.append(x[0])
    yplot = y.tolist()
    yplot.append(y[0])

    # Plot the contour and, if provided, also the sample.
    if plotted_sample:
        plot_sample(plotted_sample, ax=ax)
    if contour_label:
        ax.plot(xplot, yplot, line_style, alpha=alpha, label=contour_label)
    else:
        ax.plot(xplot, yplot, line_style, alpha=alpha)
    if median_x is not None:
        ax.plot(median_x, median_y, median_style, label=median_label)

    # Format the figure.
    if contour_label:
        plt.legend(loc='upper left', frameon=False)
    if x_label:
        plt.xlabel(x_label)
    if y_label:
        plt.ylabel(y_label)
    if x_lim:
        plt.xlim(x_lim)
    y_lim_factor = 1.2
    if plotted_sample and upper_ylim is None:
        # If there is not enough space for the legend in the upper left corner:
        # make space for it.
        max_index = np.where(plotted_sample.y == max(plotted_sample.y))
        if plotted_sample.x[max_index] < 0.6 * max(max(x), max(plotted_sample.x)):
            y_lim_factor = 1.35

        upper_ylim = max(max(y), max(plotted_sample.y)) * y_lim_factor
    elif upper_ylim is None:
        upper_ylim = max(y) * y_lim_factor
    plt.ylim((0, upper_ylim))


    # Remove axis on the right and on the top (Matlab 'box off').
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')


class PlottedSample():
    """
    Class that holds a plotted sample and its meta information.

    Attributes
    ----------
    x : ndarray of floats
        The sample's first environmental variable.
    y : ndarray of floats
        The sample's second environmental variable.
    ax : Axes
        Axes of the figure where the scatter plot should be drawn.
    label : str
        Label that will be used in the legend for the sample.
    x_inside : ndarray of floats
        Values in the first dimension of the points inside the contour.
    y_inside : ndarray of floats
        Values in the second dimension of the points inside the contour.
    x_outside : ndarray of floats
        Values in the first dimension of the points outside the contour.
    y_outside : ndarray of floats
        Values in the second dimension of the points outside the contour.
    return_period : int, optional
        Return period in years. Is used in legend for describing the inside and
        outside datapoints.
    do_plot_extreme : ndarray of booleans, optional
        Specifies which extremes should be plotted.
        The order is [min(x), max(x), min(y), max(y)].
    """

    def __init__(self, x, y, ax, label=None, x_inside=None, y_inside=None,
                 x_outside=None, y_outside=None, return_period=None,
                 do_plot_extreme=[False, False, False, False]):
        """

        Parameters
        ----------
        x : ndarray of floats
            The sample's first environmental variable.
        y : ndarray of floats
            The sample's second environmental variable.
        ax : Axes
            Axes of the figure where the scatter plot should be drawn.
        label : str
            Label that will be used in the legend for the sample.
        x_inside : ndarray of floats
            Values in the first dimension of the points inside the contour.
        y_inside : ndarray of floats
            Values in the second dimension of the points inside the contour.
        x_outside : ndarray of floats
            Values in the first dimension of the points outside the contour.
        y_outside : ndarray of floats
            Values in the second dimension of the points outside the contour.
        return_period : int, optional
            Return period in years. Is used in legend for describing the inside and
            outside datapoints.
        do_plot_extreme : ndarray of booleans, optional
            Specifies which extremes should be plotted.
            The order is [min(x), max(x), min(y), max(y)].
        """
        self.x = x
        self.y = y
        self.ax = ax
        self.label = label
        self.x_inside = x_inside
        self.y_inside = y_inside
        self.x_outside = x_outside
        self.y_outside = y_outside
        self.return_period = return_period
        self.do_plot_extreme = do_plot_extreme


def plot_confidence_interval(x_median, y_median, x_bottom, y_bottom, x_upper,
                             y_upper, ax, x_label=None, y_label=None,
                             contour_labels=[None, None, None], plotted_sample=None):
    """
    Plots the confidence interval (median, bottom, upper) in a standardized
    appearance.

    Parameters
    ----------
    x_median : ndarray of doubles
        The 50-percentile contour's coordinates in the x-direction.
    y_median : ndarray of doubles
        The 50-percentile contour's coordinates in the y-direction.
    x_bottom : ndarray of doubles
        The bottom percentile contour's coordinates in the x-direction.
    y_bottom : ndarray of doubles
        The bottom percentile contour's coordinates in the y-direction.
    x_upper : ndarray of doubles
        The upper percentile contour's coordinates in the x-direction.
    y_upper : ndarray of doubles
        The upper percentile contour's coordinates in the y-direction.
    ax : Axes
        Axes of the figure where the contour should be plotted.
    x_label : str, optional
        Label for the x-axis.
    y_label : str, optional
        Label for the y-axis.
    contour_labels : list of str, optional
        Label for th environmental contours that will be used in the legend.
    plotted_sample : PlottedSample, optional
        If provided, this sample is drawn together with the contours.
    """
    plot_contour(x=x_median,
                 y=y_median,
                 ax=ax,
                 x_label=x_label,
                 y_label=y_label,
                 line_style='b-',
                 contour_label=contour_labels[0],
                 plotted_sample=plotted_sample)
    plot_contour(x=x_bottom,
                 y=y_bottom,
                 contour_label=contour_labels[1],
                 line_style='r--',
                 ax=ax)
    plot_contour(x=x_upper,
                 y=y_upper,
                 contour_label=contour_labels[2],
                 line_style='r--',
                 ax=ax)


def plot_wave_breaking_limit(ax, bottom_tz=0, upper_tz=20, steps=100):
    """
    Plots the wave breaking limit on a given axes.

    Assumes that x = zero-up-crossing period and y = sig. wave height.

    Parameters
    ----------
    ax : Axes
    bottom_tz : float
        Bottom limit for the curve.
    upper_tz : float
        Upper limit for the curve.
    steps : int
        Number of points that are plotted on the curve.
    """
    tz_lim = np.linspace(bottom_tz, upper_tz, steps)
    hs_lim = hs_from_limiting_sig_wave_steepness(tz_lim)
    ax.plot(tz_lim, hs_lim, linestyle='-.', color=[0.5, 0.5, 0.5])


def hs_from_limiting_sig_wave_steepness(tz):
    """
    Calculates highest Hs value for a given Tz based on wave steepness.

    The calculaion uses the 'limiting significant wave steepness' described in
    DNV GL's DNV-GL-RP-C205:2017 (section 3.5.3.5. and 3.5.4)

    Parameters
    ----------
    tz : ndarray of doubles
        Zero-up-crossing period in seconds.

    Returns
    -------
    hs : ndarray of doubles
        Significant wave height in meters.
    """
    TZLOW_STEEPNESS_VALUE = 0.1 # is 1/10 in DNVG-GL-RP-C205:2017, 3.5.4
    TZHIGH_STEEPNESS_VALUE = 0.0666666  # is 1/15 in DNVG-GL-RP-C205:2017, 3.5.4

    G = 9.81

    ss = np.empty(tz.size)
    ss[:] = np.nan


    for i in range(tz.size):
        if tz[i] <= 6:
            ss[i] = 0.1
        elif tz[i] < 12:
            ss[i] = \
                TZLOW_STEEPNESS_VALUE\
                + (TZHIGH_STEEPNESS_VALUE - TZLOW_STEEPNESS_VALUE) / 6.0 * (tz[i] - 6)
        else:
            ss[i] = TZHIGH_STEEPNESS_VALUE

    hs = ss * tz**2 * G / (2 * np.pi)

    return hs
