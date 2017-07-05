
# import tool_kit
# import model
# import os
# import data_processing
# import numpy
# import datetime
# from scipy.stats import norm, beta
# from Tkinter import *
# from scipy.optimize import minimize
# from random import uniform
# import matplotlib.pyplot as plt
# from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
# import outputs
# import autumn.economics
# import itertools
# import time
# import eventlet
# from flask_socketio import emit
# from numpy import isfinite
import tb_model
from basepop import BaseModel, make_sigmoidal_curve



def generate_candidates(n_candidates, param_ranges_unc):

    """
    Function for generating candidate parameters.
    """

    # Dictionary for storing candidates
    param_candidates = {}
    for param_dict in param_ranges_unc:

        # Find bounds of parameter
        bound_low, bound_high = param_dict['bounds'][0], param_dict['bounds'][1]

        # Draw from distribution
        if param_dict['distribution'] == 'beta':
            x = bound_low + numpy.random.beta(2., 2., n_candidates) * (bound_high - bound_low)
        elif param_dict['distribution'] == 'uniform':
            x = numpy.random.uniform(bound_low, bound_high, n_candidates)
        else:
            x = 0.5*(param_dict['bounds'][0] + param_dict['bounds'][1])
            print "Unsupported statistical distribution specified to generate a candidate. Returned the midpoint of the range."

        # Return values
        param_candidates[param_dict['key']] = x
    return param_candidates


def elementwise_list_addition(increment, list_to_increment):

    """
    Simple method to element-wise increment a list by the values in another list of the same length.
    """

    assert len(increment) == len(list_to_increment), 'Attempted to add two lists of different lengths'
    return [sum(x) for x in zip(list_to_increment, increment)]


def elementwise_list_division(numerator, denominator):

    """
    Simple method to element-wise increment a list by the values in another list of the same length.
    """

    assert len(numerator) == len(denominator), 'Attempted to divide two lists of different lengths'
    return [n / d for n, d in zip(numerator, denominator)]


def find_integer_dict_from_float_dict(float_dict):

    # Method may be redundant with optimal code

    integer_dict = {}
    times = float_dict.keys()
    times.sort()
    start = numpy.floor(times[0])
    finish = numpy.floor(times[-1])
    float_years = numpy.linspace(start, finish, finish - start + 1.)
    for year in float_years:
        key = [t for t in times if t >= year][0]
        integer_dict[int(key)] = float_dict[key]
    return integer_dict


def extract_integer_dicts(models_to_analyse={}, dict_to_extract_from={}):

    # Method may be redundant with optimal code

    integer_dict = {}
    for scenario in models_to_analyse:
        integer_dict[scenario] = {}
        for output in dict_to_extract_from[scenario]:
            integer_dict[scenario][output] \
                = find_integer_dict_from_float_dict(dict_to_extract_from[scenario][output])
    return integer_dict


def get_output_dicts_from_lists(models_to_analyse={}, output_dict_of_lists={}):

    """
    Convert output lists to dictionaries. Also may ultimately be unnecessary.
    """

    output_dictionary = {}
    for scenario in models_to_analyse:
        output_dictionary[scenario] = {}
        for output in output_dict_of_lists[scenario]:
            if output != 'times':
                output_dictionary[scenario][output] \
                    = dict(zip(output_dict_of_lists[scenario]['times'], output_dict_of_lists[scenario][output]))
    return output_dictionary


def find_uncertainty_output_weights(list, method, relative_weights=[1., 2.]):

    """
    Creates a set of "weights" to determine the proportion of the log-likelihood to be contributed by the years
    considered in the calibration.

    Args:
        list: A list of the years that the weights are to be applied to.
        method: Choice of method.
        relative_weights: Relative size of the starting and ending weights if method is 1.
    """

    # Linearly scaling weights summing to one
    if method == 1:
        weights = []
        if len(list) == 1:
            return [1.]
        else:
            weights = numpy.linspace(relative_weights[0], relative_weights[1], num=len(list))
            return [i / sum(weights) for i in weights]

    # Equally distributed weights summing to one
    elif method == 2:
        return [1. / float(len(list))] * len(list)

    # All weights equal to one
    elif method == 3:
        return [1.] * len(list)


def is_parameter_value_valid(parameter):
    """
    Determine whether a number (typically a parameter value) is finite and positive.
    """

    return isfinite(parameter) and parameter > 0.


class ModelRunner:

    def __init__(self, country, fixed_parameters, time_variant_parameters={}, mode='manual', scenarios_to_run=0):
        """
        Instantiation method for model runner - currently including many attributes that should be set externally, e.g.
        in the GUI(s).

        Args:
            country: String for country being simulated
            fixed_parameters: Dictionary of parameter set used to run manual calibration
            mode: Whether scenario or uncertainty being run, set to either 'manual' or 'uncertainty'
        """

        self.country = country
        self.fixed_parameters = fixed_parameters
        self.time_variant_parameters = time_variant_parameters
        self.mode = mode
        self.scenarios_to_run = scenarios_to_run

        # loading of inputs
        # self.inputs = data_processing.Inputs(gui_inputs, runtime_outputs, js_gui=js_gui)
        # self.inputs.read_and_load_data()

        # preparing for basic runs
        self.model_dict = {}
        # self.interventions_to_cost = self.inputs.interventions_to_cost

        # Uncertainty-related attributes
        # self.is_last_run_success = False
        # self.loglikelihoods = []
        # self.outputs_unc = [{'key': 'incidence',
        #                      'posterior_width': None,
        #                      'width_multiplier': 2.  # Width of normal posterior relative to range of parameter values allowed
        #                      }]
        # self.all_parameters_tried = {}
        # self.whether_accepted_list = []
        # self.accepted_indices = []
        # self.rejected_indices = []
        # self.solns_for_extraction = ['compartment_soln', 'fraction_soln']
        # self.arrays_for_extraction = ['flow_array', 'fraction_array', 'soln_array', 'var_array', 'costs']
        # self.acceptance_dict = {}
        # self.rejection_dict = {}
        # self.uncertainty_percentiles = {}
        # self.percentiles = [2.5, 50., 97.5]
        # self.accepted_no_burn_in_indices = []
        # self.random_start = False  # whether to start from a random point, as opposed to the manually calibrated value
        #
        # # Optimisation attributes
        # self.optimisation = False  # Leave True even if loading optimisation results
        # self.indicator_to_minimise = 'incidence'  # Currently must be 'incidence' or 'mortality'
        # self.annual_envelope = [112.5e6]  # Size of funding envelope in scenarios to be run
        # self.save_opti = True
        # self.load_opti = False  # Optimisation will not be run if on
        # self.total_funding = None  # Funding for entire period
        # self.f_tol = {'incidence': 0.5,
        #               'mortality': 0.05}  # Stopping condition for optimisation algorithm (differs by indicator)
        # self.year_end_opti = 2035.  # Model is run until that date during optimisation
        # self.acceptable_combinations = []  # List of intervention combinations that can be considered with funding
        # self.opti_results = {}  # Store all the results that we need for optimisation
        # self.optimised_combinations = []
        # self.optimal_allocation = {}
        # self.interventions_considered_for_opti \
        #     = ['engage_lowquality', 'xpert', 'cxrxpertacf_prison', 'cxrxpertacf_urbanpoor', 'ipt_age0to5', 'intensive_screening']  # Interventions that must appear in optimal plan
        # self.interventions_forced_for_opti \
        #     = ['engage_lowquality', 'ipt_age0to5', 'intensive_screening']
        #
        # # Output-related attributes
        # self.epi_outputs_to_analyse = ['population', 'incidence', 'true_incidence', 'prevalence', 'true_prevalence',
        #                                'mortality', 'true_mortality', 'notifications']
        # self.epi_outputs = {}
        # self.epi_outputs_dict = {}
        # self.epi_outputs_integer_dict = {}
        # self.epi_outputs_uncertainty = {}
        # self.epi_outputs_uncertainty_centiles = None
        # self.cost_outputs = {}
        # self.cost_outputs_dict = {}
        # self.cost_outputs_integer_dict = {}
        # self.cost_outputs_uncertainty = {}
        # self.cost_outputs_uncertainty_centiles = None
        # self.additional_cost_types = ['inflated', 'discounted', 'discounted_inflated']
        # self.cost_types = self.additional_cost_types + ['raw']

    ###############################################
    ### Master methods to run all other methods ###
    ###############################################

    def master_runner(self):
        """
        Calls methods to run model with each of the three fundamental approaches.
        """

        if self.mode == 'manual':
            self.run_manual_calibration()
        elif self.mode == 'uncertainty':
            self.run_uncertainty()

    def run_manual_calibration(self):
        """
        Runs the scenarios a single time, starting from baseline with parameter values as specified in spreadsheets.
        """

        for scenario in self.scenarios_to_run:

            # Name and initialise model
            # scenario_name = 'manual_' + tool_kit.find_scenario_string_from_number(scenario)

            self.model_dict[scenario] \
                = tb_model.SimpleTbModel(fixed_parameters,
                                         scenario,
                                         time_variant_parameters=self.time_variant_parameters)

            # sort out times for scenario runs
            # if scenario is not None:
            #     scenario_name = 'manual_' + tool_kit.find_scenario_string_from_number(scenario)
            #     self.prepare_new_model_from_baseline('manual', scenario_name)

            # describe model and integrate
            print('Running scenario ' + str(scenario) + ' conditions for ' + self.country +
                  ' using single parameter set')
            self.model_dict[scenario].make_times(1900, 2050, .05)
            self.model_dict[scenario].integrate(method='explicit')

            # model interpretation for each scenario
            # self.epi_outputs[scenario_name] \
            #     = self.find_epi_outputs(scenario_name,
            #                             outputs_to_analyse=self.epi_outputs_to_analyse,
            #                             stratifications=[self.model_dict[scenario_name].agegroups,
            #                                              self.model_dict[scenario_name].riskgroups])
            # if len(self.model_dict[scenario_name].interventions_to_cost) > 0:
            #     self.find_cost_outputs(scenario_name)

    # def prepare_new_model_from_baseline(self, run_type, scenario_name):
    #
    #     """
    #     Method to set the start time of a model and load the compartment values from the baseline run.
    #
    #     Args:
    #         run_type: The type of run for the model object to be set
    #         scenario_name: Either the scenario name or optimisation if during optimisation run
    #     """
    #
    #     scenario_start_time_index = \
    #         self.model_dict[run_type + '_baseline'].find_time_index(self.inputs.model_constants['recent_time'])
    #     start_time = self.model_dict[run_type + '_baseline'].times[scenario_start_time_index]
    #     self.model_dict[scenario_name].start_time = start_time
    #     self.model_dict[scenario_name].next_time_point = start_time
    #     self.model_dict[scenario_name].loaded_compartments = \
    #         self.model_dict[run_type + '_baseline'].load_state(scenario_start_time_index)
    #
    # ####################################
    # ### Model interpretation methods ###
    # ####################################
    #
    # def find_epi_outputs(self, scenario, outputs_to_analyse, stratifications=[]):
    #     """
    #     Method to extract all requested epidemiological outputs from the models. Intended ultimately to be flexible\
    #     enough for use for analysis of scenarios, uncertainty and optimisation.
    #     """
    #
    #     epi_outputs = {'times': self.model_dict[scenario].times}
    #
    #     # Unstratified outputs______________________________________________________________________________________
    #     # Initialise lists
    #     for output in outputs_to_analyse:
    #         epi_outputs[output] = [0.] * len(epi_outputs['times'])
    #         for strain in self.model_dict[scenario].strains:
    #             epi_outputs[output + strain] = [0.] * len(epi_outputs['times'])
    #
    #     # Population
    #     if 'population' in outputs_to_analyse:
    #         for compartment in self.model_dict[scenario].compartments:
    #             epi_outputs['population'] \
    #                 = elementwise_list_addition(self.model_dict[scenario].get_compartment_soln(compartment),
    #                                             epi_outputs['population'])
    #     # Replace zeroes with small numbers for division
    #     total_denominator = tool_kit.prepare_denominator(epi_outputs['population'])
    #
    #     # To allow calculation by strain and the total output
    #     strains = self.model_dict[scenario].strains + ['']
    #
    #     # Incidence
    #     if 'incidence' in outputs_to_analyse:
    #         for strain in strains:
    #             # Variable flows
    #             for from_label, to_label, rate in self.model_dict[scenario].var_transfer_rate_flows:
    #                 if 'latent' in from_label and 'active' in to_label and strain in to_label:
    #                     incidence_increment = self.model_dict[scenario].get_compartment_soln(from_label) \
    #                                              * self.model_dict[scenario].get_var_soln(rate) \
    #                                              / total_denominator \
    #                                              * 1e5
    #                     epi_outputs['true_incidence' + strain] \
    #                         = elementwise_list_addition(incidence_increment,
    #                                                     epi_outputs['true_incidence' + strain])
    #                     # Reduce paediatric contribution
    #                     if '_age' in from_label and tool_kit.is_upper_age_limit_at_or_below(from_label, 15.):
    #                         incidence_increment *= self.inputs.model_constants['program_prop_child_reporting']
    #                     epi_outputs['incidence' + strain] \
    #                         = elementwise_list_addition(incidence_increment,
    #                                                     epi_outputs['incidence' + strain])
    #             # Fixed flows
    #             for from_label, to_label, rate in self.model_dict[scenario].fixed_transfer_rate_flows:
    #                 if 'latent' in from_label and 'active' in to_label and strain in to_label:
    #                     incidence_increment = self.model_dict[scenario].get_compartment_soln(from_label) \
    #                                              * rate / total_denominator * 1e5
    #                     epi_outputs['true_incidence' + strain] \
    #                         = elementwise_list_addition(incidence_increment,
    #                                                     epi_outputs['true_incidence' + strain])
    #                     # Reduce paedatric contribution
    #                     if '_age' in from_label and tool_kit.is_upper_age_limit_at_or_below(from_label, 15.):
    #                         incidence_increment *= self.inputs.model_constants['program_prop_child_reporting']
    #                     epi_outputs['incidence' + strain] \
    #                         = elementwise_list_addition(incidence_increment,
    #                                                     epi_outputs['incidence' + strain])
    #         # Find percentage incidence by strain
    #         if len(self.model_dict[scenario].strains) > 1:
    #             for strain in self.model_dict[scenario].strains:
    #                 epi_outputs['perc_incidence' + strain] \
    #                     = [i / j * 1e2 for i, j in zip(epi_outputs['incidence' + strain],
    #                                                    tool_kit.prepare_denominator(epi_outputs['incidence']))]
    #
    #     # Notifications
    #     if 'notifications' in outputs_to_analyse:
    #         for strain in strains:
    #             for from_label, to_label, rate in self.model_dict[scenario].var_transfer_rate_flows:
    #                 if 'active' in from_label and 'detect' in to_label and strain in to_label:
    #                     notifications_increment \
    #                         = self.model_dict[scenario].get_compartment_soln(from_label) \
    #                           * self.model_dict[scenario].get_var_soln(rate)
    #                     if '_age' in from_label and tool_kit.is_upper_age_limit_at_or_below(from_label, 15.):
    #                         notifications_increment *= self.inputs.model_constants['program_prop_child_reporting']
    #                     epi_outputs['notifications' + strain] \
    #                         = elementwise_list_addition(notifications_increment, epi_outputs['notifications' + strain])
    #
    #     # Prevalence
    #     if 'prevalence' in outputs_to_analyse:
    #         for strain in strains:
    #             for label in self.model_dict[scenario].labels:
    #                 if 'susceptible' not in label and 'latent' not in label and strain in label:
    #                     prevalence_increment = self.model_dict[scenario].get_compartment_soln(label) \
    #                                            / total_denominator \
    #                                            * 1e5
    #                     epi_outputs['true_prevalence' + strain] \
    #                         = elementwise_list_addition(prevalence_increment,
    #                                                     epi_outputs['true_prevalence' + strain])
    #                     # Reduce paediatric contribution
    #                     if '_age' in label and tool_kit.is_upper_age_limit_at_or_below(label, 15.):
    #                         prevalence_increment *= self.inputs.model_constants['program_prop_child_reporting']
    #                     epi_outputs['prevalence' + strain] \
    #                         = elementwise_list_addition(prevalence_increment,
    #                                                     epi_outputs['prevalence' + strain])
    #
    #     # Stratified outputs________________________________________________________________________________________
    #     # Currently not bothering to do this for each strain
    #     for stratification in stratifications:
    #         if len(stratification) > 1:
    #             for stratum in stratification:
    #
    #                 # Initialise lists
    #                 for output in outputs_to_analyse:
    #                     epi_outputs[output + stratum] = [0.] * len(epi_outputs['times'])
    #
    #                 # Population
    #                 if 'population' in outputs_to_analyse:
    #                     for compartment in self.model_dict[scenario].compartments:
    #                         if stratum in compartment:
    #                             epi_outputs['population' + stratum] \
    #                                 = elementwise_list_addition(self.model_dict[scenario].get_compartment_soln(compartment),
    #                                                             epi_outputs['population' + stratum])
    #
    #                 # The population denominator to be used with zeros replaced with small numbers
    #                 stratum_denominator \
    #                     = tool_kit.prepare_denominator(epi_outputs['population' + stratum])
    #
    #                 # Incidence
    #                 if 'incidence' in outputs_to_analyse:
    #                     # Variable flows
    #                     for from_label, to_label, rate in self.model_dict[scenario].var_transfer_rate_flows:
    #                         if 'latent' in from_label and 'active' in to_label and stratum in from_label:
    #                             incidence_increment = self.model_dict[scenario].get_compartment_soln(from_label) \
    #                                                   * self.model_dict[scenario].get_var_soln(rate) \
    #                                                   / stratum_denominator \
    #                                                   * 1e5
    #                             epi_outputs['true_incidence' + stratum] \
    #                                 = elementwise_list_addition(incidence_increment,
    #                                                             epi_outputs['true_incidence' + stratum])
    #                             # Reduce paediatric contribution
    #                             if '_age' in from_label \
    #                                     and tool_kit.is_upper_age_limit_at_or_below(from_label, 15.):
    #                                 incidence_increment *= self.inputs.model_constants[
    #                                     'program_prop_child_reporting']
    #                             epi_outputs['incidence' + stratum] \
    #                                 = elementwise_list_addition(incidence_increment,
    #                                                             epi_outputs['incidence' + stratum])
    #
    #                     # Fixed flows
    #                     for from_label, to_label, rate in self.model_dict[scenario].fixed_transfer_rate_flows:
    #                         if 'latent' in from_label and 'active' in to_label and stratum in from_label:
    #                             incidence_increment = self.model_dict[scenario].get_compartment_soln(from_label) \
    #                                                   * rate \
    #                                                   / stratum_denominator \
    #                                                   * 1e5
    #                             epi_outputs['true_incidence' + stratum] \
    #                                 = elementwise_list_addition(incidence_increment,
    #                                                             epi_outputs['true_incidence' + stratum])
    #                             # Reduce paediatric contribution
    #                             if '_age' in from_label \
    #                                     and tool_kit.is_upper_age_limit_at_or_below(from_label, 15.):
    #                                 incidence_increment \
    #                                     *= self.inputs.model_constants['program_prop_child_reporting']
    #                             epi_outputs['incidence' + stratum] \
    #                                 = elementwise_list_addition(incidence_increment,
    #                                                             epi_outputs['incidence' + stratum])
    #                 # Prevalence
    #                 if 'prevalence' in outputs_to_analyse:
    #                     for label in self.model_dict[scenario].labels:
    #                         if 'susceptible' not in label and 'latent' not in label and stratum in label:
    #                             prevalence_increment = self.model_dict[scenario].get_compartment_soln(label) \
    #                                                    / stratum_denominator \
    #                                                    * 1e5
    #                             epi_outputs['true_prevalence' + stratum] \
    #                                 = elementwise_list_addition(prevalence_increment,
    #                                                             epi_outputs['true_prevalence' + stratum])
    #                         # Reduce paediatric contribution
    #                         if '_age' in label and tool_kit.is_upper_age_limit_at_or_below(label, 15.):
    #                             prevalence_increment *= self.inputs.model_constants['program_prop_child_reporting']
    #                         epi_outputs['prevalence' + stratum] \
    #                             = elementwise_list_addition(prevalence_increment, epi_outputs['prevalence' + stratum])
    #
    #     return epi_outputs
    #
    # def find_population_fractions(self, stratifications=[]):
    #
    #     """
    #     Find the proportion of the population in various stratifications.
    #     The stratifications must apply to the entire population, so not to be used for strains, etc.
    #     """
    #
    #     for scenario in self.model_dict:
    #         for stratification in stratifications:
    #             if len(stratification) > 1:
    #                 for stratum in stratification:
    #                     self.epi_outputs[scenario]['fraction' + stratum] \
    #                         = elementwise_list_division(self.epi_outputs[scenario]['population' + stratum],
    #                                                     self.epi_outputs[scenario]['population'])
    #
    # def find_cost_outputs(self, scenario_name):
    #
    #     """
    #     Master method to call methods to find and update costs below.
    #
    #     Args:
    #         scenario_name: String for the name of the model being costed.
    #     """
    #
    #     self.cost_outputs[scenario_name] = self.find_raw_cost_outputs(scenario_name)
    #     self.cost_outputs[scenario_name]['raw_cost_all_programs'] = self.find_costs_all_programs(scenario_name)
    #     self.cost_outputs[scenario_name].update(self.find_adjusted_costs(scenario_name))
    #
    # def find_raw_cost_outputs(self, scenario_name):
    #
    #     """
    #     Add cost dictionaries to cost_outputs attribute.
    #     """
    #
    #     cost_outputs = {'times': self.model_dict[scenario_name].cost_times}
    #     for i, intervention \
    #             in enumerate(self.interventions_to_cost[tool_kit.find_scenario_number_from_string(scenario_name)]):
    #         cost_outputs['raw_cost_' + intervention] = self.model_dict[scenario_name].costs[:, i]
    #     return cost_outputs
    #
    # def find_costs_all_programs(self, scenario_name):
    #
    #     """
    #     Sum costs across all programs and populate to cost_outputs dictionary for each scenario.
    #     """
    #
    #     costs_all_programs \
    #         = [0.] * len(self.cost_outputs[scenario_name]['raw_cost_' + self.interventions_to_cost[
    #         tool_kit.find_scenario_number_from_string(scenario_name)][0]])
    #     for i in self.interventions_to_cost[tool_kit.find_scenario_number_from_string(scenario_name)]:
    #         costs_all_programs \
    #             = elementwise_list_addition(self.cost_outputs[scenario_name]['raw_cost_' + i], costs_all_programs)
    #     return costs_all_programs
    #
    # def find_adjusted_costs(self, scenario_name):
    #
    #     """
    #     Find costs adjusted for inflation and discounting.
    #
    #     Args:
    #         scenario_name: Scenario being costed
    #     """
    #
    #     # Get some preliminary parameters
    #     year_current = self.inputs.model_constants['recent_time']
    #     current_cpi = self.inputs.scaleup_fns[None]['econ_cpi'](year_current)
    #     discount_rate = self.inputs.model_constants['econ_discount_rate']
    #
    #     # Loop over interventions to be costed and cost types to calculate costs
    #     cost_outputs = {}
    #     for intervention in self.interventions_to_cost[tool_kit.find_scenario_number_from_string(scenario_name)] \
    #             + ['all_programs']:
    #         for cost_type in self.additional_cost_types:
    #             cost_outputs[cost_type + '_cost_' + intervention] = []
    #             for t, time in enumerate(self.cost_outputs[scenario_name]['times']):
    #                 cost_outputs[cost_type + '_cost_' + intervention].append(
    #                     autumn.economics.get_adjusted_cost(self.cost_outputs[scenario_name]['raw_cost_'
    #                                                                                         + intervention][t],
    #                                                        cost_type,
    #                                                        current_cpi,
    #                                                        self.inputs.scaleup_fns[None]['econ_cpi'](time),
    #                                                        discount_rate,
    #                                                        max(0., (time - year_current))))
    #
    #     return cost_outputs
    #
    # def find_uncertainty_centiles(self, full_uncertainty_outputs):
    #
    #     """
    #     Find percentiles from uncertainty dictionaries.
    #
    #     Updates:
    #         self.percentiles: Adds all the required percentiles to this dictionary.
    #     """
    #
    #     uncertainty_centiles = {}
    #     self.accepted_no_burn_in_indices = [i for i in self.accepted_indices if i >= self.gui_inputs['burn_in_runs']]
    #
    #     # Loop through scenarios and outputs
    #     for scenario in full_uncertainty_outputs:
    #         uncertainty_centiles[scenario] = {}
    #         for output in full_uncertainty_outputs[scenario]:
    #             if output != 'times':
    #
    #                 # To deal with the fact that we are currently saving all baseline runs but only accepted scenarios:
    #                 if scenario == 'uncertainty_baseline':
    #                     matrix_to_analyse = full_uncertainty_outputs[scenario][output][
    #                                         self.accepted_no_burn_in_indices, :]
    #                 else:
    #                     matrix_to_analyse = full_uncertainty_outputs[scenario][output]
    #
    #                 # Find the actual centiles
    #                 uncertainty_centiles[scenario][output] \
    #                     = numpy.percentile(matrix_to_analyse, self.percentiles, axis=0)
    #
    #     # Return result to make usable in other situations
    #     return uncertainty_centiles
    #
    # ###########################
    # ### Uncertainty methods ###
    # ###########################
    #
    # def run_uncertainty(self):
    #
    #     """
    #     Main method to run all the uncertainty processes.
    #     """
    #
    #     self.add_comment_to_gui_window('Uncertainty analysis commenced')
    #
    #     # If not doing an adaptive search, only need to start with a single parameter set
    #     if self.gui_inputs['adaptive_uncertainty']:
    #         n_candidates = 1
    #     else:
    #         n_candidates = self.gui_inputs['uncertainty_runs'] * 10
    #
    #     # Decide whether to start analysis from a random point or the manual values of the parameters
    #     if not self.gui_inputs['adaptive_uncertainty'] or self.random_start:
    #         param_candidates = generate_candidates(n_candidates, self.inputs.param_ranges_unc)
    #     else:
    #         param_candidates = {}
    #         for param in self.inputs.param_ranges_unc:
    #             param_candidates[param['key']] = [self.inputs.model_constants[param['key']]]
    #
    #     # Find weights for outputs that are being calibrated to
    #     normal_char = self.get_fitting_data()
    #     years_to_compare = range(1990, 2015)
    #     weights = find_uncertainty_output_weights(years_to_compare, 1, [1., 2.])
    #     self.add_comment_to_gui_window('"Weights": \n' + str(weights))
    #
    #     # Prepare for uncertainty loop
    #     n_accepted = 0
    #     prev_log_likelihood = -1e10
    #     for param in self.inputs.param_ranges_unc:
    #         self.all_parameters_tried[param['key']] = []
    #         self.acceptance_dict[param['key']] = {}
    #         self.rejection_dict[param['key']] = {}
    #         self.rejection_dict[param['key']][n_accepted] = []
    #
    #     # Instantiate uncertainty model objects
    #     for scenario in self.gui_inputs['scenarios_to_run']:
    #         scenario_name = 'uncertainty_' + tool_kit.find_scenario_string_from_number(scenario)
    #         self.model_dict[scenario_name] = model.ConsolidatedModel(scenario, self.inputs, self.gui_inputs)
    #
    #     # Until a sufficient number of parameters are accepted
    #     run = 0
    #     while n_accepted < self.gui_inputs['uncertainty_runs']:
    #
    #         # Set timer
    #         start_timer_run = datetime.datetime.now()
    #
    #         # If we are using existing parameters
    #         if run == 0 or not self.gui_inputs['adaptive_uncertainty']:
    #             new_param_list = []
    #             for param in self.inputs.param_ranges_unc:
    #                 new_param_list.append(param_candidates[param['key']][run])
    #                 params = new_param_list
    #
    #         # If we need to get a new parameter set from the old accepted set
    #         else:
    #             new_param_list = self.update_params(params)
    #
    #         # Run baseline integration (includes parameter checking, parameter setting and recording success/failure)
    #         self.run_with_params(new_param_list, model_object='uncertainty_baseline')
    #
    #         # Now storing regardless of acceptance, provided run was completed successfully
    #         if self.is_last_run_success:
    #
    #             # Get outputs for calibration and store results
    #             self.store_uncertainty('uncertainty_baseline', epi_outputs_to_analyse=self.epi_outputs_to_analyse)
    #             integer_dictionary \
    #                 = extract_integer_dicts(['uncertainty_baseline'],
    #                                         get_output_dicts_from_lists(models_to_analyse=['uncertainty_baseline'],
    #                                                                     output_dict_of_lists=self.epi_outputs))
    #
    #             # Calculate prior
    #             prior_log_likelihood = 0.
    #             for p, param in enumerate(self.inputs.param_ranges_unc):
    #                 param_val = new_param_list[p]
    #                 self.all_parameters_tried[param['key']].append(new_param_list[p])
    #
    #                 # Calculate the density of param_val
    #                 bound_low, bound_high = param['bounds'][0], param['bounds'][1]
    #
    #                 # Normalise value and find log of PDF from appropriate distribution
    #                 if param['distribution'] == 'beta':
    #                     prior_log_likelihood \
    #                         += beta.logpdf((param_val - bound_low) / (bound_high - bound_low), 2., 2.)
    #                 elif param['distribution'] == 'uniform':
    #                     prior_log_likelihood += numpy.log(1. / (bound_high - bound_low))
    #
    #             # Calculate posterior
    #             posterior_log_likelihood = 0.
    #             for output_dict in self.outputs_unc:
    #
    #                 # The GTB values for the output of interest
    #                 working_output_dictionary = normal_char[output_dict['key']]
    #                 for y, year in enumerate(years_to_compare):
    #                     if year in working_output_dictionary.keys():
    #                         model_result_for_output = integer_dictionary['uncertainty_baseline']['incidence'][year]
    #                         mu, sd = working_output_dictionary[year][0], working_output_dictionary[year][1]
    #                         posterior_log_likelihood += norm.logpdf(model_result_for_output, mu, sd) * weights[y]
    #
    #             # Determine acceptance
    #             log_likelihood = prior_log_likelihood + posterior_log_likelihood
    #             accepted = numpy.random.binomial(n=1, p=min(1., numpy.exp(log_likelihood - prev_log_likelihood)))
    #
    #             # Explain progression of likelihood
    #             self.add_comment_to_gui_window('Previous log likelihood:\n' + str(prev_log_likelihood)
    #                                            + '\nLog likelihood this run:\n' + str(log_likelihood)
    #                                            + '\nAcceptance probability:\n'
    #                                            + str(min(1., numpy.exp(log_likelihood - prev_log_likelihood)))
    #                                            + '\nWhether accepted:\n' + str(bool(accepted)) + '\n________________\n')
    #             self.loglikelihoods.append(log_likelihood)
    #
    #             # Record information for all runs
    #             if bool(accepted):
    #                 self.whether_accepted_list.append(True)
    #                 self.accepted_indices += [run]
    #                 n_accepted += 1
    #                 for p, param in enumerate(self.inputs.param_ranges_unc):
    #                     self.acceptance_dict[param['key']][n_accepted] = new_param_list[p]
    #                     self.rejection_dict[param['key']][n_accepted] = []
    #
    #                 # Update likelihood and parameter set for next run
    #                 prev_log_likelihood = log_likelihood
    #                 params = new_param_list
    #
    #                 # Run scenarios other than baseline and store uncertainty (only if accepted)
    #                 for scenario in self.gui_inputs['scenarios_to_run']:
    #                     if scenario is not None:
    #                         scenario_name = 'uncertainty_' + tool_kit.find_scenario_string_from_number(scenario)
    #                         self.prepare_new_model_from_baseline('uncertainty', scenario_name)
    #                     scenario_name = 'uncertainty_' + tool_kit.find_scenario_string_from_number(scenario)
    #                     self.run_with_params(new_param_list, model_object=scenario_name)
    #                     self.store_uncertainty(scenario_name, epi_outputs_to_analyse=self.epi_outputs_to_analyse)
    #
    #             else:
    #                 self.whether_accepted_list.append(False)
    #                 self.rejected_indices += [run]
    #                 for p, param in enumerate(self.inputs.param_ranges_unc):
    #                     self.rejection_dict[param['key']][n_accepted].append(new_param_list[p])
    #
    #             # Plot parameter progression and report on progress
    #             self.plot_progressive_parameters()
    #             self.add_comment_to_gui_window(str(n_accepted) + ' accepted / ' + str(run) +
    #                                            ' candidates. Running time: '
    #                                            + str(datetime.datetime.now() - start_timer_run))
    #             run += 1
    #
    #         # Generate more candidates if required -
    #         if not self.gui_inputs['adaptive_uncertainty'] and run >= len(param_candidates.keys()):
    #             param_candidates = generate_candidates(n_candidates, self.inputs.param_ranges_unc)
    #             run = 0
    #
    # def set_model_with_params(self, param_dict, model_object='baseline'):
    #
    #     """
    #     Populates baseline model with params from uncertainty calculations.
    #
    #     Args:
    #         param_dict: Dictionary of the parameters to be set within the model (keys parameter name strings and values
    #             parameter values).
    #     """
    #
    #     for key in param_dict:
    #         if key in self.model_dict[model_object].params:
    #             self.model_dict[model_object].set_parameter(key, param_dict[key])
    #         else:
    #             raise ValueError("%s not in model_object params" % key)
    #
    # def convert_param_list_to_dict(self, params):
    #
    #     """
    #     Extract parameters from list into dictionary that can be used for setting in the model through the
    #     set_model_with_params method.
    #
    #     Args:
    #         params: The parameter names for extraction.
    #     Returns:
    #         param_dict: The dictionary returned in appropriate format.
    #     """
    #
    #     param_dict = {}
    #     for names, vals in zip(self.inputs.param_ranges_unc, params):
    #         param_dict[names['key']] = vals
    #     return param_dict
    #
    # def get_fitting_data(self):
    #
    #     """
    #     Define the characteristics (mean and standard deviation) of the normal distribution for model outputs
    #     (incidence, mortality).
    #
    #     Returns:
    #         normal_char: Dictionary with keys outputs and values dictionaries. Sub-dictionaries have keys years
    #             and values lists, with first element of list means and second standard deviations.
    #     """
    #
    #     # Dictionary storing the characteristics of the normal distributions
    #     normal_char = {}
    #     for output_dict in self.inputs.outputs_unc:
    #         normal_char[output_dict['key']] = {}
    #
    #         # Mortality
    #         if output_dict['key'] == 'mortality':
    #             sd = output_dict['posterior_width'] / (2. * 1.96)
    #             for year in self.inputs.data_to_fit[output_dict['key']].keys():
    #                 mu = self.inputs.data_to_fit[output_dict['key']][year]
    #                 normal_char[output_dict['key']][year] = [mu, sd]
    #
    #         # Incidence
    #         elif output_dict['key'] == 'incidence':
    #             for year in self.inputs.data_to_fit[output_dict['key']].keys():
    #                 low = self.inputs.data_to_fit['incidence_low'][year]
    #                 high = self.inputs.data_to_fit['incidence_high'][year]
    #                 sd = output_dict['width_multiplier'] * (high - low) / (2. * 1.96)
    #                 mu = (high + low) / 2.
    #                 normal_char[output_dict['key']][year] = [mu, sd]
    #
    #     return normal_char
    #
    # def update_params(self, old_params):
    #
    #     """
    #     Update all the parameter values being used in the uncertainty analysis.
    #
    #     Args:
    #         old_params:
    #
    #     Returns:
    #         new_params: The new parameters to be used in the next model run.
    #
    #     """
    #
    #     new_params = []
    #
    #     # Iterate through the parameters being used
    #     for p, param_dict in enumerate(self.inputs.param_ranges_unc):
    #         bounds = param_dict['bounds']
    #         sd = self.gui_inputs['search_width'] * (bounds[1] - bounds[0]) / (2.0 * 1.96)
    #         random = -100.
    #
    #         # Search for new parameters
    #         while random < bounds[0] or random > bounds[1]:
    #             random = norm.rvs(loc=old_params[p], scale=sd, size=1)
    #
    #         # Add them to the dictionary
    #         new_params.append(random[0])
    #
    #     return new_params
    #
    # def run_with_params(self, params, model_object='uncertainty_baseline'):
    #
    #     """
    #     Integrate the model with the proposed parameter set.
    #
    #     Args:
    #         params: The parameters to be set in the model.
    #     """
    #
    #     # Check whether parameter values are acceptable
    #     for p, param in enumerate(params):
    #
    #         # Whether the parameter value is valid
    #         if not is_parameter_value_valid(param):
    #             print 'Warning: parameter%d=%f is invalid for model' % (p, param)
    #             self.is_last_run_success = False
    #             return
    #         bounds = self.inputs.param_ranges_unc[p]['bounds']
    #
    #         # Whether the parameter value is within acceptable ranges
    #         if (param < bounds[0]) or (param > bounds[1]):
    #             # print 'Warning: parameter%d=%f is outside of the allowed bounds' % (p, param)
    #             self.is_last_run_success = False
    #             return
    #
    #     param_dict = self.convert_param_list_to_dict(params)
    #
    #     # Set parameters and run
    #     self.set_model_with_params(param_dict, model_object)
    #     self.is_last_run_success = True
    #     try:
    #         self.model_dict[model_object].integrate()
    #     except:
    #         print "Warning: parameters=%s failed with model" % params
    #         self.is_last_run_success = False
    #
    # def store_uncertainty(self, scenario_name, epi_outputs_to_analyse):
    #
    #     """
    #     Add model results from one uncertainty run to the appropriate outputs dictionary, vertically stacking
    #     results on to the previous matrix.
    #
    #     Args:
    #         scenario_name: The scenario being run.
    #         epi_outputs_to_analyse: The epidemiological outputs of interest.
    #     Updates:
    #         self.epi_outputs_uncertainty
    #         self.cost_outputs_uncertainty
    #     """
    #
    #     # Get outputs
    #     self.epi_outputs[scenario_name] \
    #         = self.find_epi_outputs(scenario_name, outputs_to_analyse=self.epi_outputs_to_analyse)
    #     if len(self.model_dict[scenario_name].interventions_to_cost) > 0:
    #         self.find_cost_outputs(scenario_name)
    #
    #     # Initialise dictionaries if needed
    #     if scenario_name not in self.epi_outputs_uncertainty:
    #         self.epi_outputs_uncertainty[scenario_name] = {'times': self.epi_outputs[scenario_name]['times']}
    #         self.cost_outputs_uncertainty[scenario_name] = {'times': self.cost_outputs[scenario_name]['times']}
    #         for output in epi_outputs_to_analyse:
    #             self.epi_outputs_uncertainty[scenario_name][output] \
    #                 = numpy.empty(shape=[0, len(self.epi_outputs[scenario_name]['times'])])
    #         for output in self.cost_outputs[scenario_name]:
    #             self.cost_outputs_uncertainty[scenario_name][output] \
    #                 = numpy.empty(shape=[0, len(self.cost_outputs[scenario_name]['times'])])
    #
    #     # Add uncertainty data to dictionaries
    #     for output in epi_outputs_to_analyse:
    #         self.epi_outputs_uncertainty[scenario_name][output] \
    #             = numpy.vstack([self.epi_outputs_uncertainty[scenario_name][output],
    #                             self.epi_outputs[scenario_name][output]])
    #     for output in self.cost_outputs[scenario_name]:
    #         self.cost_outputs_uncertainty[scenario_name][output] \
    #             = numpy.vstack([self.cost_outputs_uncertainty[scenario_name][output],
    #                             self.cost_outputs[scenario_name][output]])
    #


if __name__ == '__main__':
    time_treatment = .5
    fixed_parameters = {
        'demo_rate_birth': 20. / 1e3,
        'demo_rate_death': 1. / 65,
        'tb_n_contact': 40.,
        'tb_rate_earlyprogress': .1 / .5,
        'tb_rate_lateprogress': .1 / 100.,
        'tb_rate_stabilise': .9 / .5,
        'tb_rate_recover': .6 / 3.,
        'tb_rate_death': .4 / 3.,
        'program_rate_completion_infect': .9 / time_treatment,
        'program_rate_default_infect': .05 / time_treatment,
        'program_rate_death_infect': .05 / time_treatment,
        'program_rate_completion_noninfect': .9 / time_treatment,
        'program_rate_default_noninfect': .05 / time_treatment,
        'program_rate_death_noninfect': .05 / time_treatment,
        'int_vaccine_efficacy': .5
    }
    scenarios_to_run = [0, 1]
    time_variant_parameters = {}
    for scenario in scenarios_to_run:
        time_variant_parameters[scenario] = {}
        curve1 = make_sigmoidal_curve(y_high=2, y_low=0, x_start=1950, x_inflect=1970, multiplier=4)
        curve2 = make_sigmoidal_curve(y_high=4, y_low=2, x_start=1995, x_inflect=2003, multiplier=3)
        time_variant_parameters[scenario]['rate_program_detect'] \
            = lambda x: curve1(x) if x < 1990 else curve2(x)
        time_variant_parameters[scenario]['prop_vaccination'] \
            = make_sigmoidal_curve(y_high=.95, y_low=0, x_start=1921, x_inflect=2006, multiplier=3)

    model_runner = ModelRunner(country='India',
                               fixed_parameters=fixed_parameters,
                               time_variant_parameters=time_variant_parameters,
                               mode='manual',
                               scenarios_to_run=scenarios_to_run)
    model_runner.master_runner()


