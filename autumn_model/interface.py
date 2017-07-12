
from model_runner import ModelRunner
from outputs import Project

# user inputs (which could be moved to a GUI as required)

# mode = 'manual'
mode = 'uncertainty'

country = 'India'

scenarios_to_run = [0, 1, 2]

input_parameters = {
    'demo_rate_birth': 20. / 1e3,
    'demo_rate_death': 1. / 65,
    'tb_n_contact': 60.,
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
    'time_treatment': .5
}

# dictionary of uncertainty parameters, with standard keys
param_ranges_unc = [{'name': 'tb_n_contact',
                     'start': 25.,
                     'lower_bound': 0.,
                     'upper_bound': 50.,
                     'search_width': 5.,
                     'distribution': 'beta'}]

# runnning
model_runner = ModelRunner(country=country,
                           fixed_parameters=input_parameters,
                           mode=mode,
                           scenarios_to_run=scenarios_to_run,
                           param_ranges_unc=param_ranges_unc)
model_runner.master_runner()
project = Project(model_runner)
project.master_outputs_runner()
