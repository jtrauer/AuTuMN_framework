
from model_runner import ModelRunner
from outputs import Project

# user inputs (which could be moved to a GUI as required)
mode = 'uncertainty'  # in AuTuMN_framework, must be either 'manual' or 'uncertainty'
country = 'India'  # must accord with the country string used in the Global TB Report
scenarios_to_run = [0, 1, 2]  # scenarios to be run
epi_outputs_to_analyse = ['population', 'incidence', 'prevalence']  # epidemiological outputs to be assessed
input_parameters = {'demo_rate_birth': 20. / 1e3,
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

# dictionary of uncertainty parameters, with standardised keys
# to add more parameters, add another list element with the same set of keys
param_ranges_unc = [{'name': 'tb_n_contact',
                     'start': 25.,
                     'lower_bound': 0.,
                     'upper_bound': 50.,
                     'search_width': 5.,
                     'distribution': 'beta'}]
uncertainty_accepted_runs = 4  # how many accepted runs needed before uncertainty analysis finishes (including burn-in)
burn_in = 2  # number of runs to discard (both accepted and rejected)
integration_times = [1850, 2035, .05]  # must contain three-element list of 0: start time, 1: end time, 2: time step
target_incidence = {'indicator': 'incidence', 'estimate': 150., 'sd': 30., 'year': 2016}

# code to start the model running
model_runner = ModelRunner(country,
                           input_parameters,
                           mode,
                           scenarios_to_run,
                           param_ranges_unc,
                           epi_outputs_to_analyse,
                           uncertainty_accepted_runs,
                           burn_in,
                           integration_times,
                           target_incidence)
model_runner.master_runner()
project = Project(model_runner, mode, epi_outputs_to_analyse)
project.master_outputs_runner()
