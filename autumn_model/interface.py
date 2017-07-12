
from model_runner import ModelRunner
from outputs import Project


mode = 'manual'
# mode = 'uncertainty'

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
model_runner \
    = ModelRunner(country='India', fixed_parameters=input_parameters, mode=mode, scenarios_to_run=[0, 1, 2])
model_runner.master_runner()
project = Project(model_runner)
project.master_outputs_runner()
