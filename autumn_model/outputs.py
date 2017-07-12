
import tool_kit
from matplotlib import pyplot
import platform
import os


def find_smallest_factors_of_integer(n):
    """
    Quick method to iterate through integers to find the smallest whole number fractions.
    Written only to be called by find_subplot_numbers.

    Args:
        n: Integer to be factorised
    Returns:
        answer: The two smallest factors of the integer
    """

    answer = [1e3, 1e3]
    for i in range(1, n + 1):
        if n % i == 0 and i+(n/i) < sum(answer):
            answer = [i, n/i]
    return answer


def get_nice_font_size(subplot_grid):
    """
    Simple function to return a reasonable font size as appropriate to the number of rows of subplots in the figure.
    """

    return 2. + 8. / subplot_grid[0]


def find_reasonable_year_ticks(start_time, end_time):
    """
    Function to find a reasonable spacing between years for x-ticks.

    Args:
        start_time: Float for left x-limit in years
        end_time: Float for right x-limit in years
    Returns:
        times: The times for the x ticks
    """

    duration = end_time - start_time
    if duration > 1e3:
        spacing = 1e2
    elif duration > 75.:
        spacing = 25.
    elif duration > 25.:
        spacing = 10.
    elif duration > 15.:
        spacing = 5.
    elif duration > 5.:
        spacing = 2.
    else:
        spacing = 1.

    times = []
    working_time = start_time
    while working_time < end_time:
        times.append(working_time)
        working_time += spacing
    return times


def find_standard_output_styles(labels, lightening_factor=1.):
    """
    Function to find some standardised colours for the outputs we'll typically
    be reporting on - i.e. incidence, prevalence, mortality and notifications.
    Incidence is black/grey, prevalence green, mortality red and notifications blue.

    Args:
        labels: List containing strings for the outputs that colours are needed for.
        lightening_factor: Float between zero and one that specifies how much lighter to make
            the colours - with 0. being no additional lightening (black or dark green/red/blue)
            and 1. being completely lightened to reach white.
    Returns:
        colour: Colour for plotting
        indices: List of strings to be used to find the data in the data object
        yaxis_label: Unit of measurement for outcome
        title: Title for plot (so far usually a subplot)
        patch_colour: Colour half way between colour and white
    """

    colour = []
    indices = []
    yaxis_label = []
    title = []
    patch_colour = []

    if 'incidence' in labels:
        colour += [(lightening_factor, lightening_factor, lightening_factor)]
        indices += ['e_inc_100k']
        yaxis_label += ['Per 100,000 per year']
        title += ['Incidence']
    if 'mortality' in labels:
        colour += [(1., lightening_factor, lightening_factor)]
        indices += ['e_mort_exc_tbhiv_100k']
        yaxis_label += ['Per 100,000 per year']
        title += ['Mortality']
    if 'prevalence' in labels:
        colour += [(lightening_factor, 0.5 + 0.5 * lightening_factor, lightening_factor)]
        indices += ['e_prev_100k']
        yaxis_label += ['Per 100,000']
        title += ['Prevalence']
    if 'notifications' in labels:
        colour += [(lightening_factor, lightening_factor, 0.5 + 0.5 * lightening_factor)]
        yaxis_label += ['']
        title += ['Notifications']
    if 'perc_incidence' in labels:
        colour += [(lightening_factor, lightening_factor, lightening_factor)]
        yaxis_label += ['Percentage']
        title += ['Proportion of incidence']

    # create a colour half-way between the line colour and white for patches
    for i in range(len(colour)):
        patch_colour += [[]]
        for j in range(len(colour[i])):
            patch_colour[i] += [1. - (1. - colour[i][j]) / 2.]

    return colour, indices, yaxis_label, title, patch_colour


def find_subplot_numbers(n):

    """
    Method to find a good number of rows and columns for subplots of figure.

    Args:
        n: Total number of subplots.
    Returns:
        answer: List of two elements, being the rows and columns of the subplots.

    """

    # Find a nice number of subplots for a panel plot
    answer = find_smallest_factors_of_integer(n)
    i = 0
    while i < 10:
        if abs(answer[0] - answer[1]) > 3:
            n = n + 1
            answer = find_smallest_factors_of_integer(n)
        i = i + 1

    return answer


def scale_axes(vals, max_val, y_sig_figs):

    """
    General function to scale a set of axes and produce text that can be added to the axis label. Written here as a
    separate function from the tidy_axis method below because it can then be applied to both x- and y-axes.

    Args:
        vals: List of the current y-ticks
        max_val: The maximum value of this list
        y_sig_figs: The number of significant figures for the ticks
    Returns:
        labels: List of the modified tick labels
        axis_modifier: The text to be added to the axis
    """

    y_number_format = '%.' + str(y_sig_figs) + 'f'
    if max_val < 5e-9:
        labels = [y_number_format % (v * 1e12) for v in vals]
        axis_modifier = 'Trillionth '
    elif max_val < 5e-6:
        labels = [y_number_format % (v * 1e9) for v in vals]
        axis_modifier = 'Billionth '
    elif max_val < 5e-3:
        labels = [y_number_format % (v * 1e6) for v in vals]
        axis_modifier = 'Millionth '
    elif max_val < 5:
        labels = [y_number_format % (v * 1e3) for v in vals]
        axis_modifier = 'Thousandth '
    elif max_val < 5e3:
        labels = [y_number_format % v for v in vals]
        axis_modifier = ''
    elif max_val < 5e6:
        labels = [y_number_format % (v / 1e3) for v in vals]
        axis_modifier = 'Thousand '
    elif max_val < 5e9:
        labels = [y_number_format % (v / 1e6) for v in vals]
        axis_modifier = 'Million '
    else:
        labels = [y_number_format % (v / 1e9) for v in vals]
        axis_modifier = 'Billion '

    return labels, axis_modifier


class Project:

    def __init__(self, model_runner, epi_outputs_to_analyse=[]):
        """
        Initialises an object of class Project, that will contain all the information (data + outputs) for writing a
        report for a country.

        Args:
            models: dictionary such as: models = {'baseline': model, 'scenario_1': model_1,  ...}
        """

        self.epi_outputs_to_analyse = epi_outputs_to_analyse

        self.model_runner = model_runner
        # self.inputs = self.model_runner.inputs
        self.country = self.model_runner.country
        self.name = 'test_' + self.country
        self.out_dir_project = os.path.join('projects', self.name)
        if not os.path.isdir(self.out_dir_project):
            os.makedirs(self.out_dir_project)
        self.figure_number = 1
        self.output_colours = {}
        self.program_colours = {}
        self.suptitle_size = 13
        self.grid = False
        # self.plot_rejected_runs = False

        self.scenarios = self.model_runner.scenarios_to_run
        self.gtb_available_outputs = ['incidence', 'prevalence']
        # self.level_conversion_dict = {'lower_limit': '_lo', 'upper_limit': '_hi', 'point_estimate': ''}
        #
        # comes up so often that we need to find this index, that best to do once in instantiation
        self.start_time_index \
            = tool_kit.find_first_list_element_at_least_value(self.model_runner.epi_outputs[0]['times'], 1990)

    #################################
    # General methods for use below #
    #################################

    def set_and_update_figure(self):
        """
        If called at the start of each plotting function, will create a figure that is numbered according to
        self.figure_number, which is then updated at each call. This stops figures plotting on the same axis
        and saves you having to worry about how many figures are being opened.
        """

        fig = pyplot.figure(self.figure_number)
        self.figure_number += 1
        return fig

    def make_single_axis(self, fig):
        """
        Create axes for a figure with a single plot with a reasonable amount of space around.

        Returns:
            ax: The axes that can be plotted on
        """

        ax = fig.add_axes([0.1, 0.1, 0.6, 0.75])
        return ax

    def make_default_line_styles(self, n, return_all=True):
        """
        Produces a standard set of line styles that isn't adapted to the data being plotted.

        Args:
            n: The number of line-styles
            return_all: Whether to return all of the styles up to n or just the last one
        Returns:
            line_styles: A list of standard line-styles, or if return_all is False,
                then the single item (for methods that are iterating through plots.
        """

        # iterate through a standard set of line styles
        for i in range(n):
            line_styles = []
            for line in ["-", ":", "-.", "--"]:
                for colour in "krbgmcy":
                    line_styles.append(line + colour)

        if return_all:
            return line_styles
        else:
            return line_styles[n - 1]

    def tidy_axis(self, ax, subplot_grid, title='', start_time=0., legend=False, x_label='', y_label='',
                  x_axis_type='time', y_axis_type='scaled', x_sig_figs=0, y_sig_figs=0):
        """
        Method to make cosmetic changes to a set of plot axes.
        """

        # add the sub-plot title with slightly larger titles than the rest of the text on the panel
        if title: ax.set_title(title, fontsize=get_nice_font_size(subplot_grid) + 2.)

        # add a legend if needed
        if legend == 'for_single':
            ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., frameon=False, prop={'size': 7})
        elif legend:
            ax.legend(fontsize=get_nice_font_size(subplot_grid), frameon=False)

        # sort x-axis
        if x_axis_type == 'time':
            ax.set_xlim((start_time, 2035))
            ax.set_xticks(find_reasonable_year_ticks(start_time, 2035))
        elif x_axis_type == 'scaled':
            vals = list(ax.get_xticks())
            max_val = max([abs(v) for v in vals])
            labels, axis_modifier = scale_axes(vals, max_val, x_sig_figs)
            ax.set_xticklabels(labels)
            ax.set_xlabel(x_label + axis_modifier, fontsize=get_nice_font_size(subplot_grid), labelpad=1)
        elif x_axis_type == 'proportion':
            ax.set_xlabel(x_label, fontsize=get_nice_font_size(subplot_grid), labelpad=1)
            ax.set_xlim((0., 1.))
        elif x_axis_type == 'individual_years':
            ax.set_xlim((start_time, self.inputs.model_constants['plot_end_time']))
            ax.set_xticks(range(int(start_time), int(self.inputs.model_constants['plot_end_time']), 1))
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_rotation(45)
        else:
            ax.set_xlabel(x_label, fontsize=get_nice_font_size(subplot_grid), labelpad=1)

        # sort y-axis
        ax.set_ylim(bottom=0.)
        vals = list(ax.get_yticks())
        max_val = max([abs(v) for v in vals])
        if y_axis_type == 'time':
            ax.set_ylim((start_time, self.inputs.model_constants['plot_end_time']))
            ax.set_yticks(find_reasonable_year_ticks(start_time, self.inputs.model_constants['plot_end_time']))
        elif y_axis_type == 'scaled':
            labels, axis_modifier = scale_axes(vals, max_val, y_sig_figs)
            ax.set_yticklabels(labels)
            ax.set_ylabel(axis_modifier + y_label, fontsize=get_nice_font_size(subplot_grid), labelpad=1)
        elif y_axis_type == 'proportion':
            ax.set_ylim(top=1.)
            ax.set_ylabel(y_label, fontsize=get_nice_font_size(subplot_grid), labelpad=1)
        else:
            ax.set_ylim(top=max_val * 1.1)
            ax.set_ylabel(y_label, fontsize=get_nice_font_size(subplot_grid), labelpad=1)

        # set size of font for x-ticks and add a grid if requested
        for axis_to_change in [ax.xaxis, ax.yaxis]:
            for tick in axis_to_change.get_major_ticks():
                tick.label.set_fontsize(get_nice_font_size(subplot_grid))
            axis_to_change.grid(self.grid)

    def save_figure(self, fig, last_part_of_name_for_figure):
        """
        Simple method to standardise names for output figure files.

        Args:
            last_part_of_name_for_figure: The part of the figure name that is variable and input from the
                plotting method.
        """

        png = os.path.join(self.out_dir_project, self.country + last_part_of_name_for_figure + '.png')
        fig.savefig(png, dpi=300)

    def save_opti_figure(self, fig, last_part_of_name_for_figure):

        """
        Same as previous method, when applied to optimisation outputs.
        """

        png = os.path.join(self.model_runner.opti_outputs_dir, self.country + last_part_of_name_for_figure + '.png')
        fig.savefig(png, dpi=300)

    #####################################################
    ### Methods for outputting to Office applications ###
    #####################################################

    def master_outputs_runner(self):
        """
        Method to work through all the fundamental output methods, which then call all the specific output
        methods for plotting and writing as required.
        """

        # master plotting method
        self.run_plotting()
        self.open_output_directory()

    ########################
    ### Plotting methods ###
    ########################

    def run_plotting(self):
        """
        Master plotting method to call all the methods that produce specific plots.
        """

        # find some general output colours
        output_colours = self.make_default_line_styles(5, True)
        for scenario in self.scenarios:
            self.output_colours[scenario] = output_colours[scenario]
            self.program_colours[scenario] = output_colours[scenario]

        # plot main outputs
        self.plot_epi_outputs(self.epi_outputs_to_analyse, ci_plot=None)

    def plot_epi_outputs(self, outputs, ci_plot=None):
        """
        Produces the plot for the main outputs, loops over multiple scenarios.

        Args:
            outputs: A list of the outputs to be plotted.
            ci_plot: Whether to plot uncertainty intervals around the estimates generated from uncertainty runs.
        """

        # standard preliminaries
        start_time = 2000
        colour, indices, yaxis_label, title, patch_colour = find_standard_output_styles(outputs, lightening_factor=0.3)
        subplot_grid = find_subplot_numbers(len(outputs) - 1)
        fig = self.set_and_update_figure()

        # loop through indicators
        for o, output in enumerate(outputs[1:]):

            # preliminaries
            ax = fig.add_subplot(subplot_grid[0], subplot_grid[1], o + 1)

            # plot scenarios without uncertainty
            if ci_plot is None:

                # plot model estimates
                for scenario in self.scenarios[::-1]:  # reversing ensures black baseline plotted over top
                    ax.plot(self.model_runner.epi_outputs[scenario]['times'][self.start_time_index:],
                            self.model_runner.epi_outputs[scenario][output][self.start_time_index:],
                            color=self.output_colours[scenario][1], linestyle=self.output_colours[scenario][0],
                            linewidth=1.5, label='Scenario ' + str(scenario))

            # # plot with uncertainty confidence intervals
            # elif ci_plot and self.gui_inputs['output_uncertainty']:
            #     for scenario in self.scenarios[::-1]:
            #         scenario_name = tool_kit.find_scenario_string_from_number(scenario)
            #         start_index = self.find_start_index(scenario)
            #
            #         # median
            #         ax.plot(
            #             self.model_runner.epi_outputs_uncertainty['uncertainty_' + scenario_name][
            #                 'times'][start_index:],
            #             self.model_runner.epi_outputs_uncertainty_centiles['uncertainty_' + scenario_name][
            #                 output][self.model_runner.percentiles.index(50), :][start_index:],
            #             color=self.output_colours[scenario][1], linestyle=self.output_colours[scenario][0],
            #             linewidth=1.5, label=tool_kit.capitalise_and_remove_underscore(scenario_name))
            #
            #         # upper and lower confidence bounds
            #         for centile in [2.5, 97.5]:
            #             ax.plot(
            #                 self.model_runner.epi_outputs_uncertainty['uncertainty_' + scenario_name]['times'][
            #                     start_index:],
            #                 self.model_runner.epi_outputs_uncertainty_centiles['uncertainty_' + scenario_name][output][
            #                 self.model_runner.percentiles.index(centile), :][start_index:],
            #                 color=self.output_colours[scenario][1], linestyle='--', linewidth=.5, label=None)
            #         end_filename = '_ci'
            #
            # # plot progressive model run outputs for baseline scenario
            # elif self.gui_inputs['output_uncertainty']:
            #     for run in range(len(self.model_runner.epi_outputs_uncertainty['uncertainty_baseline'][output])):
            #         if run not in self.model_runner.accepted_indices and self.plot_rejected_runs:
            #             ax.plot(self.model_runner.epi_outputs_uncertainty[
            #                         'uncertainty_baseline']['times'][self.start_time_index:],
            #                     self.model_runner.epi_outputs_uncertainty[
            #                         'uncertainty_baseline'][output][run, self.start_time_index:],
            #                     linewidth=.2, color='y', label=tool_kit.capitalise_and_remove_underscore('baseline'))
            #         else:
            #             ax.plot(self.model_runner.epi_outputs_uncertainty[
            #                         'uncertainty_baseline']['times'][self.start_time_index:],
            #                     self.model_runner.epi_outputs_uncertainty[
            #                         'uncertainty_baseline'][output][run, self.start_time_index:],
            #                     linewidth=1.2,
            #                     color=str(1. - float(run) / float(len(
            #                         self.model_runner.epi_outputs_uncertainty['uncertainty_baseline'][output]))),
            #                     label=tool_kit.capitalise_and_remove_underscore('baseline'))
            #         end_filename = '_progress'

            self.tidy_axis(ax, subplot_grid, title=title[o], start_time=start_time,
                           legend=(o == len(outputs) - 1 and len(self.scenarios) > 1),
                           y_axis_type='raw', y_label=yaxis_label[o])

        # add main title and save
        fig.suptitle(self.country + ' model outputs', fontsize=self.suptitle_size)
        self.save_figure(fig, '_gtb')

    def plot_likelihoods(self):
        """
        Method to plot likelihoods over runs, differentiating accepted and rejected runs to illustrate progression.
        """

        # plotting prelims
        fig = self.set_and_update_figure()
        ax = fig.add_subplot(1, 1, 1)

        # find accepted likelihoods
        accepted_log_likelihoods = [self.model_runner.loglikelihoods[i] for i in self.model_runner.accepted_indices]

        # plot the rejected values
        for i in self.model_runner.rejected_indices:

            # find the index of the last accepted index before the rejected one we're currently interested in
            last_acceptance_before = [j for j in self.model_runner.accepted_indices if j < i][-1]

            # plot from the previous acceptance to the current rejection
            ax.plot([last_acceptance_before, i],
                    [self.model_runner.loglikelihoods[last_acceptance_before],
                     self.model_runner.loglikelihoods[i]], marker='o', linestyle='--', color='.5')

        # plot the accepted values
        ax.plot(self.model_runner.accepted_indices, accepted_log_likelihoods, marker='o', color='k')

        # finishing up
        fig.suptitle('Progression of likelihood', fontsize=self.suptitle_size)
        ax.set_xlabel('All runs', fontsize=get_nice_font_size([1, 1]), labelpad=1)
        ax.set_ylabel('Likelihood', fontsize=get_nice_font_size([1, 1]), labelpad=1)
        self.save_figure(fig, '_likelihoods')

    ############################
    ### Miscellaneous method ###
    ############################

    def open_output_directory(self):
        """
        Opens the directory into which all the outputs have been placed.
        """

        operating_system = platform.system()
        if 'Windows' in operating_system:
            os.system('start ' + ' ' + self.out_dir_project)
        elif 'Darwin' in operating_system:
            os.system('open ' + ' ' + self.out_dir_project)


