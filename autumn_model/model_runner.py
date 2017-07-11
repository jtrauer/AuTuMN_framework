
import tb_model
import data_processing
from basepop import BaseModel, make_sigmoidal_curve
from outputs import Project

def prepare_denominator(list_to_prepare):
    """
    Method to safely divide a list of numbers while ignoring zero denominators.

    Args:
        list_to_prepare: The list to be used as a denominator
    Returns:
        The list with zeros replaced with small numbers
    """

    return [list_to_prepare[i] if list_to_prepare[i] > 0. else 1e-10 for i in range(len(list_to_prepare))]


def elementwise_list_addition(increment, list_to_increment):
    """
    Simple method to element-wise increment a list by the values in another list of the same length.
    """

    assert len(increment) == len(list_to_increment), 'Attempted to add two lists of different lengths'
    return [sum(x) for x in zip(list_to_increment, increment)]


class ModelRunner:

    def __init__(self, country, fixed_parameters, mode='manual', scenarios_to_run=[0]):
        """
        Instantiation method for model runner.

        Args:
            country: String for country being simulated
            fixed_parameters: Dictionary of parameter set used to run manual calibration
            mode: Whether scenario or uncertainty being run, set to either 'manual' or 'uncertainty'
        """

        self.country = country
        self.fixed_parameters = fixed_parameters
        # self.time_variant_parameters = time_variant_parameters
        self.mode = mode
        self.scenarios_to_run = scenarios_to_run

        self.inputs = data_processing.Inputs(self.country, self.scenarios_to_run)
        self.inputs.read_and_load_data()

        # loading of inputs
        # self.inputs.read_and_load_data()

        # preparing for basic runs
        self.model_dict = {}

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
        # output-related attributes
        self.epi_outputs_to_analyse = ['population', 'incidence', 'prevalence']
        self.epi_outputs = {}
        # self.epi_outputs_uncertainty = {}

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

            # name and initialise model
            self.model_dict[scenario] \
                = tb_model.SimpleTbModel(fixed_parameters, self.inputs, scenario)

            # describe model and integrate
            print('Running scenario ' + str(scenario) + ' conditions for ' + self.country +
                  ' using single parameter set')
            self.model_dict[scenario].make_times(1850, 2035, .05)
            self.model_dict[scenario].integrate(method='explicit')

            # model interpretation for each scenario
            self.epi_outputs[scenario] \
                = self.find_epi_outputs(scenario)

    ####################################
    ### Model interpretation methods ###
    ####################################

    def find_epi_outputs(self, scenario):
        """
        Method to extract requested epidemiological outputs from the models.
        """

        # first create a list of model times as de facto keys for the series of lists created below
        epi_outputs = {'times': self.model_dict[scenario].times}

        # initialise lists
        for output in self.epi_outputs_to_analyse:
            epi_outputs[output] = [0.] * len(epi_outputs['times'])

        # population
        if 'population' in self.epi_outputs_to_analyse:
            for compartment in self.model_dict[scenario].compartments:
                epi_outputs['population'] \
                    = elementwise_list_addition(self.model_dict[scenario].get_compartment_soln(compartment),
                                                epi_outputs['population'])

        # replace zeroes with small numbers for division
        total_denominator = prepare_denominator(epi_outputs['population'])

        # incidence
        if 'incidence' in self.epi_outputs_to_analyse:

            # fixed flows
            for from_label, to_label, rate in self.model_dict[scenario].fixed_transfer_rate_flows:
                if 'latent' in from_label and 'active' in to_label:
                    incidence_increment \
                        = self.model_dict[scenario].get_compartment_soln(from_label) * rate / total_denominator * 1e5
                    epi_outputs['incidence'] = elementwise_list_addition(incidence_increment, epi_outputs['incidence'])

            # variable flows (note that there are currently none to which this is applicable, but could be in theory)
            for from_label, to_label, rate in self.model_dict[scenario].var_transfer_rate_flows:
                if 'latent' in from_label and 'active' in to_label:
                    incidence_increment \
                        = self.model_dict[scenario].get_compartment_soln(from_label) \
                          * self.model_dict[scenario].get_var_soln(rate) / total_denominator * 1e5
                    epi_outputs['incidence'] = elementwise_list_addition(incidence_increment, epi_outputs['incidence'])

        # prevalence
        if 'prevalence' in self.epi_outputs_to_analyse:
            for label in self.model_dict[scenario].labels:
                if 'susceptible' not in label and 'latent' not in label:
                    prevalence_increment = self.model_dict[scenario].get_compartment_soln(label) \
                                           / total_denominator * 1e5
                    epi_outputs['prevalence'] \
                        = elementwise_list_addition(prevalence_increment, epi_outputs['prevalence'])

        return epi_outputs

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
        = ModelRunner(country='India', fixed_parameters=fixed_parameters, mode='manual', scenarios_to_run=[0, 1, 2])
    model_runner.master_runner()
    project = Project(model_runner)
    project.master_outputs_runner()


