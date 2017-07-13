
import tb_model
import data_processing
from basepop import BaseModel, make_sigmoidal_curve
from outputs import Project
from scipy.stats import norm, beta
import numpy

"""
This module provides an object-oriented structure for running model objects. Manual calibration and uncertainty are
provided as examples (being largely user-coded, with few dependencies for the main running processes).
"""


def prepare_denominator(list_to_prepare):
    """
    Method to safely divide a list of numbers while ignoring zero denominators by adding a very small value to any
    zeroes.

    Args:
        list_to_prepare: The list to be used as a denominator
    Returns:
        The list with zeros replaced with small numbers
    """

    return [list_to_prepare[i] if list_to_prepare[i] > 0. else 1e-10 for i in range(len(list_to_prepare))]


def elementwise_list_addition(increment, list_to_increment):
    """
    Simple method to element-wise increment a list by the values in another list of the same length
    (as is needed in output generation).

    Args:
        increment: A list of values to be added to the previous list
        list_to_increment: The original list to be incremented
    Returns:
        The resulting list with elements being the sum of the elements of the two lists
    """

    assert len(increment) == len(list_to_increment), 'Attempted to increment lists with two lists of different lengths'
    return [sum(x) for x in zip(list_to_increment, increment)]


class ModelRunner:
    """
    Object to coordinate running of all the functions of the platform, such as manual scenario running and uncertainty.
    Stores model objects as attributes to itself.
    """

    def __init__(self, country, fixed_parameters, time_variant_parameters, mode='manual',
                 param_ranges_unc=[], epi_outputs_to_analyse=[], scenario_implementation=[],
                 uncertainty_accepted_runs=50, burn_in=5,
                 integration_times={'start': 1900, 'finish': 2035, 'step': .05},
                 target={'indicator': 'incidence', 'estimate': 150., 'sd': 30., 'year': 2016}):
        """
        Instantiation method for model runner.

        Args:
            country: String for country being simulated
            fixed_parameters: Dictionary of parameter set used to run manual calibration
            mode: Whether scenario or uncertainty being run, set to either 'manual' or 'uncertainty'
            scenarios_to_run: List of values for scenarios to be simulated (with 0 being baseline)
            param_ranges_unc: List of dictionaries for the uncertainty parameters to be considered
            epi_outputs_to_analyse: List of strings for the epidemiological outcomes of interest (e.g. incidence)
            scenario_implementation: List of dictionaries for scenarios to be implemented as described in interface
            uncertainty_accepted_runs: How many accepted uncertainty runs are required
            burn_in: How many runs to be discarded as a burn in
            target: Dictionary describing values for model to be calibrated against
        """

        # convert arguments to attributes
        self.country = country
        self.fixed_parameters = fixed_parameters
        self.time_variant_parameters = time_variant_parameters
        self.mode = mode
        self.param_ranges_unc = param_ranges_unc
        self.epi_outputs_to_analyse = epi_outputs_to_analyse
        self.scenario_implementation = scenario_implementation
        self.scenarios_to_run = range(len(scenario_implementation))
        self.uncertainty_accepted_runs = uncertainty_accepted_runs
        self.burn_in = burn_in
        self.integration_times = integration_times
        self.target = target

        # inputs obtained from spreadsheet reading and data processing
        self.inputs = data_processing.Inputs(self.country, self.fixed_parameters,
                                             self.time_variant_parameters, self.scenario_implementation)
        self.inputs.read_and_load_data()

        # dictionary for storing models
        self.model_dict = {}

        # output-related attributes
        self.epi_outputs = {}
        self.epi_outputs_uncertainty = []
        self.epi_outputs_uncertainty_centiles = {}

        # uncertainty-related attributes
        self.accepted_indices = []
        self.is_last_run_success = True

    ###############################################
    ### Master methods to run all other methods ###
    ###############################################

    def master_runner(self):
        """
        Calls methods to run model with either of the two fundamental approaches presented here.
        """

        if self.mode == 'manual':
            self.run_manual_calibration()
        elif self.mode == 'uncertainty':
            self.run_uncertainty()
            self.find_uncertainty_centiles()

    def run_manual_calibration(self):
        """
        Runs each of the scenarios a single time, starting from baseline scenario with fixed parameter values.
        """

        for scenario in self.scenarios_to_run:

            # name and initialise model
            self.model_dict[scenario] = tb_model.SimpleTbModel(self.fixed_parameters, self.inputs, scenario)

            # describe model and integrate
            print('Running scenario ' + str(scenario) + ' conditions for ' + self.country +
                  ' using single parameter set')
            self.model_dict[scenario].make_times(self.integration_times['start'],
                                                 self.integration_times['finish'],
                                                 self.integration_times['step'])
            self.model_dict[scenario].integrate(method='explicit')

            # find epidemiological model outputs for each scenario
            self.epi_outputs[scenario] = self.find_epi_outputs(scenario)

    ####################################
    ### Model interpretation methods ###
    ####################################

    def find_epi_outputs(self, scenario):
        """
        Method to extract requested epidemiological outputs from the models. Returns a data object (rather than adding
        data directly to self) in order that it can be used for both scenario running and uncertainty.
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
        total_denominator = prepare_denominator(epi_outputs['population'])

        # incidence
        if 'incidence' in self.epi_outputs_to_analyse:

            # fixed flows
            for from_label, to_label, rate in self.model_dict[scenario].fixed_transfer_rate_flows:
                if 'latent' in from_label and 'active' in to_label:
                    incidence_increment \
                        = self.model_dict[scenario].get_compartment_soln(from_label) * rate / total_denominator * 1e5
                    epi_outputs['incidence'] = elementwise_list_addition(incidence_increment, epi_outputs['incidence'])

            # variable flows (note that there are currently none to which this is applicable, but could be)
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

    def find_uncertainty_centiles(self):
        """
        Find percentiles from uncertainty dictionaries.
        """

        accepted_no_burn_in_indices = [i for i in self.accepted_indices if i > self.burn_in]
        for output in self.epi_outputs_to_analyse:
            self.epi_outputs_uncertainty_centiles[output] \
                = numpy.percentile(self.epi_outputs_uncertainty[output][accepted_no_burn_in_indices, :],
                                   [2.5, 50., 97.5], axis=0)

    ###########################
    ### Uncertainty methods ###
    ###########################

    def run_uncertainty(self):
        """
        Main method to run all the uncertainty processes.
        """

        # prepare for uncertainty loop
        print('Uncertainty analysis commenced')
        n_accepted = 0
        prev_log_likelihood = -1e10
        run = 0
        self.model_dict['uncertainty'] = tb_model.SimpleTbModel(self.fixed_parameters, self.inputs, 0)

        # find initial set of parameters
        new_param_list = []
        for i in range(len(self.param_ranges_unc)):
            new_param_list.append(self.param_ranges_unc[i]['start'])
            params = new_param_list

        # until a sufficient number of parameters are accepted
        while n_accepted < self.uncertainty_accepted_runs:

            # run baseline integration
            self.run_with_params(new_param_list)

            # store regardless of acceptance (provided model ran successfully)
            if self.is_last_run_success:

                # calculate uncertainty outputs and store results
                self.store_uncertainty()

                # calculate prior
                prior_log_likelihood = 0.
                for i in range(len(self.param_ranges_unc)):
                    param_val = new_param_list[i]
                    bound_low, bound_high \
                        = self.param_ranges_unc[i]['lower_bound'], self.param_ranges_unc[i]['upper_bound']

                    # normalise and find log PDF from appropriate distribution
                    if self.param_ranges_unc[i]['distribution'] == 'beta':
                        prior_log_likelihood += beta.logpdf((param_val - bound_low) / (bound_high - bound_low), 2., 2.)
                    elif self.param_ranges_unc[i]['distribution'] == 'uniform':
                        prior_log_likelihood += numpy.log(1. / (bound_high - bound_low))

                # calculate posterior
                epi_result = self.epi_outputs['uncertainty'][self.target['indicator']][
                    self.find_time_index(self.target['year'], 'uncertainty')]
                posterior_log_likelihood = norm.logpdf(self.epi_outputs['uncertainty'][self.target['indicator']][
                    self.find_time_index(self.target['year'], 'uncertainty')], self.target['estimate'], self.target['sd'])

                # determine acceptance
                log_likelihood = prior_log_likelihood + posterior_log_likelihood
                accepted = numpy.random.binomial(n=1, p=min(1., numpy.exp(log_likelihood - prev_log_likelihood)))

                # update likelihood, parameters and record acceptance for next run
                if bool(accepted):
                    n_accepted += 1
                    prev_log_likelihood = log_likelihood
                    params = new_param_list
                    self.accepted_indices += [run]

                # update run number
                run += 1

                # report on progress
                print('run')
                print(run)
                print('accepted')
                print(accepted)
                print('incidence')
                print(epi_result)
                for i in range(len(self.param_ranges_unc)):
                    print(self.param_ranges_unc[i]['name'])
                    print(new_param_list[i])
                print('\n')

            # obtain a new parameter list for the next run
            new_param_list = self.update_params(params)

    def run_with_params(self, params):
        """
        Run the model with the proposed parameter set.

        Args:
            params: The parameters to be set in the model.
        """

        param_dict = self.convert_param_list_to_dict(params)

        # set parameters and run
        for key in param_dict: self.model_dict['uncertainty'].set_param(key, param_dict[key])
        self.is_last_run_success = True
        try:
            self.model_dict['uncertainty'].make_times(self.integration_times['start'],
                                                      self.integration_times['finish'],
                                                      self.integration_times['step'])
            self.model_dict['uncertainty'].integrate()
        except:
            print "Warning: parameters=%s failed with model" % params
            self.is_last_run_success = False

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
        for i in range(len(self.param_ranges_unc)): param_dict[self.param_ranges_unc[i]['name']] = params[i]
        return param_dict

    def update_params(self, old_params):
        """
        Update all the parameter values being used in the uncertainty analysis.

        Args:
            old_params:

        Returns:
            new_params: The new parameters to be used in the next model run.
        """

        new_params = []

        # iterate through each uncertainty parameter
        for i in range(len(self.param_ranges_unc)):
            search_width = self.param_ranges_unc[i]['search_width']
            random = -100.

            # search for new parameters and add to list
            while random < self.param_ranges_unc[i]['lower_bound'] or random > self.param_ranges_unc[i]['upper_bound']:
                random = norm.rvs(loc=old_params[i], scale=search_width, size=1)
            new_params.append(float(random[0]))

        return new_params

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

    def find_time_index(self, time, model):
        """
        General method to find first time point in times list at or after a certain specified time point.

        Args:
            time: Float for the time point of interest.
        """

        return [i for i, j in enumerate(self.model_dict[model].times) if j >= time][0] - 1
        raise ValueError('Time not found')

