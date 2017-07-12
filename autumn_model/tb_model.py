import os
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir))
sys.path.append('C:/Users/James/Desktop/popdynamics/popdynamics/')
from basepop import BaseModel


########################################################
### Static functions for graphing and managing files ###
########################################################


def make_plots(model, out_dir):
    """
    Make some basic graphs of the scaling time-variant parameters and the basic model outputs.

    Args:
        model: Instance of the model object to be interrogated
        out_dir: The directory to put the graphs in
    """

    import pylab

    # scaling case detection rate
    pylab.clf()
    scaleup_fn = model.scaleup_fns['program_rate_detect']
    y_vals = map(scaleup_fn, model.times)
    pylab.plot(model.times, y_vals)
    pylab.title('scaleup test curve')
    pylab.savefig(os.path.join(out_dir, 'scaleup.png'))

    # main epidemiological outputs
    pylab.clf()
    for var_key in ['mortality', 'incidence', 'prevalence']:
        soln = model.get_var_soln(var_key)
        pylab.plot(model.times, soln, label=var_key)
    pylab.legend()
    pylab.savefig(os.path.join(out_dir, 'fraction.png'))


#################################
### Define model object class ###
#################################

class SimpleTbModel(BaseModel):
    """
    Initial TB model by James Trauer.
    Nested inheritance from BaseModel, which applies to any infectious disease generally.
    """

    def __init__(self, fixed_parameters, inputs, scenario=0):
        """
        Inputs:
            interventions: List of interventions to be simulated in the run of the model
        """

        BaseModel.__init__(self)

        self.inputs = inputs
        self.scenario = scenario
        self.riskgroups = ['']
        self.riskgroup_proportions = {'': 1.}
        # self.riskgroups = ['_nocomorb', '_hiv']
        # self.riskgroup_proportions = {'_nocomorb': .9, '_hiv': .1}

        # define all compartments, initialise as empty and then populate
        model_compartments \
            = ['susceptible', 'susceptible_vaccinated', 'latent_early', 'latent_late', 'active', 'treatment_infect',
               'treatment_noninfect']
        for each_compartment in model_compartments:
            for riskgroup in self.riskgroups:
                self.set_compartment(each_compartment + riskgroup, 0.)
        self.set_compartment('susceptible' + self.riskgroups[0], 1e6)
        self.set_compartment('active' + self.riskgroups[0], 1.)

        # parameter setting
        for parameter, value in fixed_parameters.items():
            self.set_param(parameter, value)
        for parameter in self.inputs.scaleup_fns[self.scenario]:
            self.set_scaleup_fn(parameter, self.inputs.scaleup_fns[self.scenario][parameter])

    def calculate_vars(self):
        """
        Calculate values that change with time over the course of model integration.
        """

        # demographic
        self.vars['population'] = sum(self.compartments.values())
        for riskgroup in self.riskgroups:
            self.vars['rate_birth' + riskgroup] \
                = self.params['demo_rate_birth'] * self.vars['population'] * self.riskgroup_proportions[riskgroup]

        # infection
        self.vars['infectious_population'] = 0.
        for label in self.labels:
            if 'active' in label or '_infect' in label:
                self.vars['infectious_population'] += self.compartments[label]
        self.vars['rate_force'] = \
            self.params['tb_n_contact'] * self.vars['infectious_population'] / self.vars['population']
        self.vars['rate_force_vaccinated'] = self.vars['rate_force'] * self.params['int_vaccine_efficacy']

        # vaccination
        for riskgroup in self.riskgroups:
            self.vars['rate_birth_vaccinated' + riskgroup] \
                = self.vars['rate_birth' + riskgroup] * self.vars['prop_vaccination']
            self.vars['rate_birth' + riskgroup] -= self.vars['rate_birth_vaccinated' + riskgroup]

        # detection
        self.vars['program_rate_detect'] \
            = self.vars['program_prop_detect'] / (1. - self.vars['program_prop_detect']) \
              * (self.params['tb_rate_death'] + self.params['tb_rate_recover'])

    def set_flows(self):
        """
        Set inter-compartmental flows, whether time-variant or constant
        """

        for riskgroup in self.riskgroups:

            # demographic
            self.set_var_entry_rate_flow(
                'susceptible' + riskgroup,
                'rate_birth' + riskgroup)
            self.set_background_death_rate(
                'demo_rate_death')
            self.set_var_entry_rate_flow(
                'susceptible_vaccinated' + riskgroup,
                'rate_birth_vaccinated' + riskgroup)

            # infection
            self.set_var_transfer_rate_flow(
                'susceptible' + riskgroup,
                'latent_early' + riskgroup,
                'rate_force')
            self.set_var_transfer_rate_flow(
                'susceptible_vaccinated' + riskgroup,
                'latent_early' + riskgroup,
                'rate_force_vaccinated')

            # natural history of infection and disease
            self.set_fixed_transfer_rate_flow(
                'latent_early' + riskgroup,
                'active' + riskgroup,
                'tb_rate_earlyprogress')
            self.set_fixed_transfer_rate_flow(
                'latent_early' + riskgroup,
                'latent_late' + riskgroup,
                'tb_rate_stabilise')
            self.set_fixed_transfer_rate_flow(
                'latent_late' + riskgroup,
                'active' + riskgroup,
                'tb_rate_lateprogress')
            self.set_fixed_transfer_rate_flow(
                'active' + riskgroup,
                'latent_late' + riskgroup,
                'tb_rate_recover')
            self.set_infection_death_rate_flow(
                'active' + riskgroup,
                'tb_rate_death')

            # programmatic
            self.set_var_transfer_rate_flow(
                'active' + riskgroup,
                'treatment_infect' + riskgroup,
                'program_rate_detect')
            self.set_fixed_transfer_rate_flow(
                'treatment_infect' + riskgroup,
                'treatment_noninfect' + riskgroup,
                'program_rate_completion_infect')
            self.set_fixed_transfer_rate_flow(
                'treatment_infect' + riskgroup,
                'active' + riskgroup,
                'program_rate_default_infect')
            self.set_fixed_transfer_rate_flow(
                'treatment_noninfect' + riskgroup,
                'susceptible' + riskgroup,
                'program_rate_completion_noninfect')
            self.set_fixed_transfer_rate_flow(
                'treatment_noninfect' + riskgroup,
                'active' + riskgroup,
                'program_rate_default_noninfect')
            self.set_infection_death_rate_flow(
                'treatment_infect' + riskgroup,
                'program_rate_death_infect')
            self.set_infection_death_rate_flow(
                'treatment_noninfect' + riskgroup,
                'program_rate_death_noninfect')

    def calculate_diagnostic_vars(self):
        """
        Calculate output variables from model quantities.
        """

        # main epidemiological indicators
        self.vars['prevalence'] = self.vars['infectious_population'] / self.vars['population'] * 1e5
        self.vars['mortality'] = self.vars['rate_infection_death'] / self.vars['population'] * 1e5
        rate_incidence = 0.
        for from_label, to_label, rate in self.fixed_transfer_rate_flows:
            val = self.compartments[from_label] * rate
            if 'latent' in from_label and 'active' in to_label:
                rate_incidence += val
        self.vars['incidence'] = rate_incidence / self.vars['population'] * 1e5

        # proportion latently infected
        self.vars['latent'] = 0.
        for label in self.labels:
            if 'latent' in label:
                self.vars['latent'] += (self.compartments[label] / self.vars['population'] * 1e5)

