
from model_runner import ModelRunner
from outputs import Project

"""
This script runs the AuTuMN_framework platform, setting all the input parameters, selections and output options required
by the remainder of the platform. As all the user interaction can occur through this script, graphical user interfaces
can be constructed out of the code contained here.
"""

# user inputs __________________________________________________________________________________________________________

# basic model features__________________________________________________________________________________________________
mode = 'manual'  # must be either 'manual' or 'uncertainty'
country = 'India'  # must accord with the country string used in the Global TB Report
epi_outputs_to_analyse = ['population', 'incidence', 'prevalence']  # epidemiological outputs to be assessed
plot_start_time = 2005  # left border of x-axes on output plots
fixed_parameters = {'demo_rate_birth': 20. / 1e3,
                    'demo_rate_death': 1. / 65,
                    'tb_n_contact': 20.,
                    'tb_rate_earlyprogress': .1 / .5,
                    'tb_rate_lateprogress': .1 / 20.,
                    'tb_rate_stabilise': .9 / .5,
                    'tb_rate_recover': .6 / 3.,
                    'tb_rate_death': .4 / 3.,
                    'program_prop_completion_infect': .98,
                    'program_prop_default_infect': .01,
                    'program_prop_death_infect': .01,
                    'program_prop_completion_noninfect': .9,
                    'program_prop_default_noninfect': .05,
                    'program_prop_death_noninfect': .05,
                    'int_vaccine_efficacy': .5,
                    'time_early_treatment': 1. / 52.,
                    'time_treatment': .5}  # fixed value input parameters for the model (some needed further processing)

# user-defined time-variant parameter dictionary with keys interventions and values dictionaries for coverage levels
# coverage dictionaries have keys years and values for parameter values
time_variant_parameters = {'prop_vaccination': {1921: 0., 1980: .8, 2015: .85}}

# scenario input________________________________________________________________________________________________________
# scenario implementation, list of dictionaries with first element None for the baseline scenario and then standardised
# keys as shown:
scenario_implementation = [None,
                           {'intervention': 'prop_vaccination', 'year': 2020, 'coverage': .99},
                           {'intervention': 'program_prop_detect', 'year': 2017, 'coverage': .9}]
scenarios_to_run = range(len(scenario_implementation))

# uncertainty inputs____________________________________________________________________________________________________
# param_ranges_unc is a dictionary of uncertainty parameters, with standardised keys. To add more parameters,
# add another list element with the same set of keys as the following:
param_ranges_unc = [{'name': 'tb_n_contact',
                     'start': 25.,
                     'lower_bound': 0.,
                     'upper_bound': 50.,
                     'search_width': 5.,
                     'distribution': 'beta'}]
uncertainty_accepted_runs = 50  # how many accepted runs needed before uncertainty analysis finishes (including burn-in)
burn_in = 5  # number of runs to discard (both accepted and rejected)
integration_times = {'start': 1850,
                     'finish': 2035,
                     'step': .05}  # dictionary for constructing integration times with fixed keys expected by runner
target_incidence = {'indicator': 'incidence',
                    'estimate': 150.,
                    'sd': 30.,
                    'year': 2016}  # dictionary for output comparison with fixed keys expected by model runner

# code to start the model running (not for user interaction)____________________________________________________________
model_runner = ModelRunner(country, fixed_parameters, time_variant_parameters, mode, scenarios_to_run, param_ranges_unc,
                           epi_outputs_to_analyse, scenario_implementation, uncertainty_accepted_runs, burn_in,
                           integration_times, target_incidence)
model_runner.master_runner()
project = Project(model_runner, plot_start_time)
project.master_outputs_runner()
