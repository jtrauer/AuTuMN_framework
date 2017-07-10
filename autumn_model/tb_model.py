import os
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir))
sys.path.append('C:/Users/James/Desktop/popdynamics/popdynamics/')
from basepop import BaseModel

"""
This file presents a transmission dynamic model for TB based on those most commonly used by the AuTuMN team
See http://www.tb-modelling.com/home/index.php for more information on AuTuMN.

Similar to seir_model.py, it uses methods from the BaseModel class in the basepop.py file from this module
(one directory above) to create the model objects for two similar TB transmission dynamic models.

It advances upon the seir_model.py examples by incorporating time-variant functions used to represent changes in
the programmatic response to TB over time, so that their epidemiological impact on the endemic burden can be
illustrated.

The first section of this code (to line 96) presents static functions for use in the master execution below.

The second section of the code (from line 99 to line 241) presents the creation of the two TB model classes, first with
varying case detection rate and second also adding the intervention of BCG vaccination.

The last section of the code (from line 244 to the end of the script) presents the execution of the example models and
calls the functions to graph their outputs.

Suggestion to get started:
- Adjust some parameters within parameter values set in model object instantiation and observe how model outputs change.
- Try adding another intervention in a similar way to that presented for BCG vaccination (e.g. treatment of LTBI,
 active case finding).
"""


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

    def __init__(self, fixed_parameters, scenario=0, time_variant_parameters={}):
        """
        Inputs:
            interventions: List of interventions to be simulated in the run of the model
        """

        BaseModel.__init__(self)

        self.scenario = scenario
        self.time_variant_parameters = time_variant_parameters

        # define all compartments, initialise as empty and then populate
        model_compartments \
            = ['susceptible', 'susceptible_vaccinated', 'latent_early', 'latent_late', 'active', 'treatment_infect',
               'treatment_noninfect']
        for each_compartment in model_compartments:
            self.set_compartment(each_compartment, 0.)
        self.set_compartment('susceptible', 1e6)
        self.set_compartment('active', 1.)

        # parameter setting
        for parameter, value in fixed_parameters.items():
            self.set_param(parameter, value)
        for parameter in self.time_variant_parameters[scenario]:
            self.set_scaleup_fn(parameter, self.time_variant_parameters[scenario][parameter])

    def calculate_vars(self):
        """
        Calculate values that change with time over the course of model integration.
        """

        # demographic
        self.vars['population'] = sum(self.compartments.values())
        self.vars['rate_birth'] = self.params['demo_rate_birth'] * self.vars['population']

        # infection
        self.vars['infectious_population'] = 0.
        for label in self.labels:
            if 'active' in label or '_infect' in label:
                self.vars['infectious_population'] += self.compartments[label]
        self.vars['rate_force'] = \
            self.params['tb_n_contact'] * self.vars['infectious_population'] / self.vars['population']

        self.vars['rate_birth_vaccinated'] = self.vars['rate_birth'] * self.vars['prop_vaccination']
        self.vars['rate_birth'] -= self.vars['rate_birth_vaccinated']
        self.vars['rate_force_vaccinated'] = self.vars['rate_force'] * self.params['int_vaccine_efficacy']

    def set_flows(self):
        """
        Set inter-compartmental flows, whether time-variant or constant
        """

        # demographic
        self.set_var_entry_rate_flow('susceptible', 'rate_birth')
        self.set_background_death_rate('demo_rate_death')

        # infection
        self.set_var_transfer_rate_flow('susceptible', 'latent_early', 'rate_force')

        # natural history of infection and disease
        self.set_fixed_transfer_rate_flow('latent_early', 'active', 'tb_rate_earlyprogress')
        self.set_fixed_transfer_rate_flow('latent_early', 'latent_late', 'tb_rate_stabilise')
        self.set_fixed_transfer_rate_flow('latent_late', 'active', 'tb_rate_lateprogress')
        self.set_fixed_transfer_rate_flow('active', 'latent_late', 'tb_rate_recover')
        self.set_infection_death_rate_flow('active', 'tb_rate_death')

        # programmatic
        self.set_var_transfer_rate_flow('active', 'treatment_infect', 'rate_program_detect')
        self.set_fixed_transfer_rate_flow('treatment_infect', 'treatment_noninfect', 'program_rate_completion_infect')
        self.set_fixed_transfer_rate_flow('treatment_infect', 'active', 'program_rate_default_infect')
        self.set_fixed_transfer_rate_flow('treatment_noninfect', 'susceptible', 'program_rate_completion_noninfect')
        self.set_fixed_transfer_rate_flow('treatment_noninfect', 'active', 'program_rate_default_noninfect')
        self.set_infection_death_rate_flow('treatment_infect', 'program_rate_death_infect')
        self.set_infection_death_rate_flow('treatment_noninfect', 'program_rate_death_noninfect')

        self.set_var_entry_rate_flow('susceptible_vaccinated', 'rate_birth_vaccinated')
        self.set_var_transfer_rate_flow('susceptible_vaccinated', 'latent_early', 'rate_force_vaccinated')

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

