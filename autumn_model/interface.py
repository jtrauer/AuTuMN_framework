
from model_runner import ModelRunner
from outputs import Project


mode = 'manual'
# mode = 'uncertainty'



time_early_treatment = 1. / 52.
time_late_treatment = .5 - time_early_treatment
fixed_parameters = {
    'demo_rate_birth': 20. / 1e3,
    'demo_rate_death': 1. / 65,
    'tb_n_contact': 60.,
    'tb_rate_earlyprogress': .1 / .5,
    'tb_rate_lateprogress': .1 / 20.,
    'tb_rate_stabilise': .9 / .5,
    'tb_rate_recover': .6 / 3.,
    'tb_rate_death': .4 / 3.,
    'program_rate_completion_infect': .98 / time_early_treatment,
    'program_rate_default_infect': .01 / time_early_treatment,
    'program_rate_death_infect': .01 / time_early_treatment,
    'program_rate_completion_noninfect': .9 / time_late_treatment,
    'program_rate_default_noninfect': .05 / time_late_treatment,
    'program_rate_death_noninfect': .05 / time_late_treatment,
    'int_vaccine_efficacy': .5
}
model_runner \
    = ModelRunner(country='India', fixed_parameters=fixed_parameters, mode=mode, scenarios_to_run=[0, 1, 2])
model_runner.master_runner()
project = Project(model_runner)
project.master_outputs_runner()
