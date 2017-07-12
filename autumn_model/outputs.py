
import tool_kit
from matplotlib import pyplot
import platform
import os


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


def find_standard_output_styles(labels):
    """
    Function to find some standardised colours for the outputs we'll typically
    be reporting on - i.e. incidence, prevalence, mortality and notifications.
    Incidence is black/grey, prevalence green, mortality red and notifications blue.

    Args:
        labels: List containing strings for the outputs that colours are needed for.
    Returns:
        yaxis_label: Unit of measurement for outcome
        title: Title for plot (so far usually a subplot)
    """

    yaxis_label = {}
    title = {}
    if 'incidence' in labels:
        yaxis_label['incidence'] = 'Per 100,000 per year'
        title['incidence'] = 'Incidence'
    if 'prevalence' in labels:
        yaxis_label['prevalence'] = 'Per 100,000'
        title['prevalence'] = 'Prevalence'
    return yaxis_label, title


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

    def __init__(self, model_runner):
        """
        Initialises an object of class Project, that will contain all the information (data + outputs) for writing a
        report for a country.

        Args:
            model_runner: The entire model runner object that was responsible for coordinating all model runs
        """

        self.model_runner = model_runner
        self.mode = self.model_runner.mode
        self.epi_outputs_to_analyse = self.model_runner.epi_outputs_to_analyse
        self.country = self.model_runner.country
        self.name = 'test_' + self.country
        self.out_dir_project = os.path.join('projects', self.name)
        if not os.path.isdir(self.out_dir_project):
            os.makedirs(self.out_dir_project)
        self.output_colours = {}
        self.program_colours = {}
        self.suptitle_size = 13
        self.grid = False
        self.scenarios = self.model_runner.scenarios_to_run

        if self.mode == 'manual':
            self.start_time_index \
                = tool_kit.find_first_list_element_at_least_value(self.model_runner.epi_outputs[0]['times'], 1990)
        elif self.mode == 'uncertainty':
            self.start_time_index \
                = tool_kit.find_first_list_element_at_least_value(self.model_runner.epi_outputs['uncertainty']['times'],
                                                                  1990)

    def make_default_line_styles(self, n):
        """
        Produces a standard set of line styles that isn't adapted to the data being plotted.

        Args:
            n: The number of line-styles
            return_all: Whether to return all of the styles up to n or just the last one
        Returns:
            line_styles: A list of standard line-styles, or if return_all is False,
                then the single item (for methods that are iterating through plots.
        """

        for i in range(n):
            line_styles = []
            for line in ["-", ":", "-.", "--"]:
                for colour in "krbgmcy":
                    line_styles.append(line + colour)
        return line_styles

    def tidy_axis(self, ax, title='', start_time=0., legend=False, x_label='', y_label='',
                  x_axis_type='time', y_axis_type='scaled', x_sig_figs=0, y_sig_figs=0):
        """
        Method to make cosmetic changes to a set of plot axes.
        """

        # add the sub-plot title with slightly larger titles than the rest of the text on the panel
        if title: ax.set_title(title, fontsize=12.)

        # add a legend if needed
        if legend == 'for_single':
            ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., frameon=False, prop={'size': 7})
        elif legend:
            ax.legend(fontsize=10., frameon=False)

        # sort x-axis
        if x_axis_type == 'time':
            ax.set_xlim((start_time, 2035))
            ax.set_xticks(find_reasonable_year_ticks(start_time, 2035))
        elif x_axis_type == 'scaled':
            vals = list(ax.get_xticks())
            max_val = max([abs(v) for v in vals])
            labels, axis_modifier = scale_axes(vals, max_val, x_sig_figs)
            ax.set_xticklabels(labels)
            ax.set_xlabel(x_label + axis_modifier, fontsize=10., labelpad=1)
        elif x_axis_type == 'proportion':
            ax.set_xlabel(x_label, fontsize=10., labelpad=1)
            ax.set_xlim((0., 1.))
        elif x_axis_type == 'individual_years':
            ax.set_xlim((start_time, self.inputs.model_constants['plot_end_time']))
            ax.set_xticks(range(int(start_time), int(self.inputs.model_constants['plot_end_time']), 1))
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_rotation(45)
        else:
            ax.set_xlabel(x_label, fontsize=10., labelpad=1)

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
            ax.set_ylabel(axis_modifier + y_label, fontsize=10., labelpad=1)
        elif y_axis_type == 'proportion':
            ax.set_ylim(top=1.)
            ax.set_ylabel(y_label, fontsize=10., labelpad=1)
        else:
            ax.set_ylim(top=max_val * 1.1)
            ax.set_ylabel(y_label, fontsize=10., labelpad=1)

        # set size of font for x-ticks and add a grid if requested
        for axis_to_change in [ax.xaxis, ax.yaxis]:
            for tick in axis_to_change.get_major_ticks():
                tick.label.set_fontsize(10.)
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

    def master_outputs_runner(self):
        """
        Method to work through all the fundamental output methods, which here call only the methods relevant to plotting
        as an example.
        """

        # master plotting method
        self.run_plotting()
        self.open_output_directory()

    def run_plotting(self):
        """
        Master plotting method to call all the methods that produce specific plots.
        """

        # find some general output colours
        output_colours = self.make_default_line_styles(5)
        for scenario in self.scenarios:
            self.output_colours[scenario] = output_colours[scenario]
            self.program_colours[scenario] = output_colours[scenario]

        # plot main outputs
        self.plot_epi_outputs(self.epi_outputs_to_analyse)

    def plot_epi_outputs(self, outputs):
        """
        Produces the plot for the main outputs, loops over multiple scenarios.

        Args:
            outputs: A list of the outputs to be plotted.
            ci_plot: Whether to plot uncertainty intervals around the estimates generated from uncertainty runs.
        """

        # standard preliminaries
        start_time = 2000
        yaxis_label, title = find_standard_output_styles(outputs)
        fig = pyplot.figure(1)

        # loop through indicators (ignoring population)
        for o, output in enumerate(outputs[1:]):

            # preliminaries
            ax = fig.add_subplot(1, len(outputs) - 1, o + 1)

            # plotting (manual or uncertainty)
            if self.mode == 'manual':
                for scenario in self.scenarios[::-1]:  # reversing ensures black baseline plotted over top
                    ax.plot(self.model_runner.epi_outputs[scenario]['times'][self.start_time_index:],
                            self.model_runner.epi_outputs[scenario][output][self.start_time_index:],
                            color=self.output_colours[scenario][1], linestyle=self.output_colours[scenario][0],
                            linewidth=1.5, label='Scenario ' + str(scenario))
            elif self.mode == 'uncertainty':
                for centile in range(3):
                    linewidths = {0: .5, 1: 1.5, 2: .5}
                    ax.plot(self.model_runner.epi_outputs_uncertainty['times'][self.start_time_index:],
                            self.model_runner.epi_outputs_uncertainty_centiles[output][centile][self.start_time_index:],
                            linewidth=linewidths[centile], color='k')

            self.tidy_axis(ax, title=title[output], start_time=start_time,
                           legend=(o == len(outputs) - 2 and len(self.scenarios) > 1),
                           y_axis_type='raw', y_label=yaxis_label[output])

        # add main title and save
        fig.suptitle(self.country + ' model outputs', fontsize=self.suptitle_size)
        self.save_figure(fig, '_outputs')

    def open_output_directory(self):
        """
        Opens the directory into which all the outputs have been placed.
        """

        operating_system = platform.system()
        if 'Windows' in operating_system:
            os.system('start ' + ' ' + self.out_dir_project)
        elif 'Darwin' in operating_system:
            os.system('open ' + ' ' + self.out_dir_project)


