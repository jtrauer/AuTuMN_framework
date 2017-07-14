
from basepop import BaseModel

"""
A simple TB model for use in this framework platform. It is built upon the BaseModel class inherited from the basepop
repository (which is a dependency). This version is considerably simpler than most country-level AuTuMN levels, having
one stratification example (by risk group). The example only provides the code structure for elaborating the model and
the risk groups do not have any function in this version. Further stratifications can be added by including loops for
age groups, organ involvement (e.g. smear-positive, smear-negative, extrapulmonary), health system (e.g. public and
private), etc. Note that not all such stratifications would be incorporated in exactly the same way - for example, the
loop for organ involvement would only apply to compartments representing patients with active disease.
"""


class TbModel(BaseModel):
    """
    Initial TB model. Nested inheritance from BaseModel, which applies to any infectious disease generally.
    """

    def __init__(self, fixed_parameters, inputs, scenario=0, additional_riskgroups={}):
        """
        Inputs:
            fixed_parameters: Fixed constant model parameters (including those that can be over-ridden in uncertainty)
            inputs: Other user and spreadsheet inputs
        """

        BaseModel.__init__(self)

        self.inputs = inputs
        self.scenario = scenario
        self.riskgroups = additional_riskgroups.keys()
        self.riskgroups.append('')
        self.riskgroup_proportions = {'': 1. - sum(additional_riskgroups.values())}
        self.riskgroup_proportions.update(additional_riskgroups)

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
        Calculate values that change with time over the course of model integration, which require information from the
        model for their calculation.
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
        Set inter-compartmental flows (whether time-variant or constant).
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
            if 'latent' in from_label and 'active' in to_label: rate_incidence += val
        self.vars['incidence'] = rate_incidence / self.vars['population'] * 1e5

        # proportion latently infected
        self.vars['latent'] = 0.
        for label in self.labels:
            if 'latent' in label: self.vars['latent'] += (self.compartments[label] / self.vars['population'] * 1e5)

