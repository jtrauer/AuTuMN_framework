
import tb_model
import data_processing
from basepop import BaseModel, make_sigmoidal_curve
from outputs import Project
from scipy.stats import norm, beta
import numpy


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
        self.epi_outputs_uncertainty = []

        # parameter name, mean of prior, sd of prior, lower bound, upper bound
        self.param_ranges_unc = [['tb_n_contact', 30., 15., 10., 50., 5., 'beta']]

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

    ###########################
    ### Uncertainty methods ###
    ###########################

    def run_uncertainty(self):
        """
        Main method to run all the uncertainty processes.
        """

        print('Uncertainty analysis commenced')

        # if not doing an adaptive search, only need to start with a single parameter set
        n_candidates = 1

        # prepare for uncertainty loop
        n_accepted = 0
        prev_log_likelihood = -1e10

        # instantiate uncertainty model object
        self.model_dict['uncertainty'] = tb_model.SimpleTbModel(fixed_parameters, self.inputs, 0)

        # until a sufficient number of parameters are accepted
        run = 0
        uncertainty_runs = 20
        while n_accepted < uncertainty_runs:

            # if we are using the first parameter set
            if run == 0:
                new_param_list = []
                for i in range(len(self.param_ranges_unc)):
                    new_param_list.append(self.param_ranges_unc[i][1])
                    params = new_param_list

            # if we need to get a new parameter set from the old accepted set
            else:
                new_param_list = self.update_params(params)

            # run baseline integration (includes parameter checking, parameter setting and recording success/failure)
            self.run_with_params(new_param_list)

            # store regardless of acceptance, provided run was completed successfully
            if self.is_last_run_success:

                # get outputs for calibration and store results
                self.store_uncertainty()

                # calculate prior
                prior_log_likelihood = 0.
                for i in range(len(self.param_ranges_unc)):
                    param_val = new_param_list[i]

                    # calculate the density of param_val
                    bound_low, bound_high = self.param_ranges_unc[i][3], self.param_ranges_unc[i][4]

                    # normalise value and find log of PDF from appropriate distribution
                    if self.param_ranges_unc[i][6] == 'beta':
                        prior_log_likelihood += beta.logpdf((param_val - bound_low) / (bound_high - bound_low), 2., 2.)
                    elif self.param_ranges_unc[i][6] == 'uniform':
                        prior_log_likelihood += numpy.log(1. / (bound_high - bound_low))

                # calculate posterior
                posterior_log_likelihood = 0.
                for output_dict in self.outputs_unc:

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

    def set_model_with_params(self, param_dict):
        """
        Populates baseline model with params from uncertainty calculations.

        Args:
            param_dict: Dictionary of the parameters to be set within the model (keys parameter name strings and values
                parameter values).
        """

        for key in param_dict:
            self.model_dict['uncertainty'].set_param(key, param_dict[key])

    def convert_param_list_to_dict(self, params):
        """
        Extract parameters from list into dictionary that can be used for setting in the model through the
        set_model_with_params method.

        Args:
            params: The parameter names for extraction.
        Returns:
            param_dict: The dictionary returned in appropriate format.
        """

        param_dict = {}
        for i in range(len(self.param_ranges_unc)):
            param_dict[self.param_ranges_unc[i][0]] = params[i]
        return param_dict

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

    def update_params(self, old_params):
        """
        Update all the parameter values being used in the uncertainty analysis.

        Args:
            old_params:

        Returns:
            new_params: The new parameters to be used in the next model run.
        """

        new_params = []

        # iterate through the parameters being used
        for i in range(len(self.param_ranges_unc)):
            search_width = self.param_ranges_unc[i][5]
            random = -100.

            # search for new parameters
            while random < self.param_ranges_unc[i][3] or random > self.param_ranges_unc[i][4]:
                random = norm.rvs(loc=old_params[i], scale=search_width, size=1)

            # add them to the list
            new_params.append(float(random[0]))

        return new_params

    def run_with_params(self, params):
        """
        Integrate the model with the proposed parameter set.

        Args:
            params: The parameters to be set in the model.
        """

        param_dict = self.convert_param_list_to_dict(params)

        # set parameters and run
        self.set_model_with_params(param_dict)
        self.is_last_run_success = True
        try:
            self.model_dict['uncertainty'].make_times(1850, 2035, .05)
            self.model_dict['uncertainty'].integrate()
        except:
            print "Warning: parameters=%s failed with model" % params
            self.is_last_run_success = False

    def store_uncertainty(self):
        """
        Add model results from one uncertainty run to the appropriate outputs dictionary, vertically stacking
        results on to the previous matrix.

        Updates:
            self.epi_outputs_uncertainty
        """

        # get outputs
        self.epi_outputs['uncertainty'] = self.find_epi_outputs('uncertainty')

        # initialise dictionaries if needed
        if not self.epi_outputs_uncertainty:
            self.epi_outputs_uncertainty = {'times': self.epi_outputs['uncertainty']['times']}
            for output in self.epi_outputs_to_analyse:
                self.epi_outputs_uncertainty[output] \
                    = numpy.empty(shape=[0, len(self.epi_outputs['uncertainty']['times'])])

        # add uncertainty data to dictionaries
        for output in self.epi_outputs_to_analyse:
            self.epi_outputs_uncertainty[output] \
                = numpy.vstack([self.epi_outputs_uncertainty[output], self.epi_outputs['uncertainty'][output]])


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
        = ModelRunner(country='India', fixed_parameters=fixed_parameters, mode='uncertainty', scenarios_to_run=[0, 1, 2])
    model_runner.master_runner()
    project = Project(model_runner)
    project.master_outputs_runner()


