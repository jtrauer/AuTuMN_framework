
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

    def __init__(self, country, fixed_parameters, mode='manual', scenarios_to_run=[0], param_ranges_unc=[],
                 epi_outputs_to_analyse=[], uncertainty_accepted_runs=100):
        """
        Instantiation method for model runner.

        Args:
            country: String for country being simulated
            fixed_parameters: Dictionary of parameter set used to run manual calibration
            mode: Whether scenario or uncertainty being run, set to either 'manual' or 'uncertainty'
        """

        self.country = country
        self.fixed_parameters = fixed_parameters
        self.mode = mode
        self.scenarios_to_run = scenarios_to_run
        self.param_ranges_unc = param_ranges_unc
        self.epi_outputs_to_analyse = epi_outputs_to_analyse
        self.uncertainty_accepted_runs = uncertainty_accepted_runs

        self.accepted_indices = []

        self.inputs = data_processing.Inputs(self.country, self.scenarios_to_run, self.fixed_parameters)
        self.inputs.read_and_load_data()

        # preparing for model runs
        self.model_dict = {}

        # uncertainty-related attributes
        self.is_last_run_success = True

        # output-related attributes
        self.epi_outputs = {}
        self.epi_outputs_uncertainty = []

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
                = tb_model.SimpleTbModel(self.fixed_parameters, self.inputs, scenario)

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

        target_incidence = [150., 30., 2016]

        # prepare for uncertainty loop
        n_accepted = 0
        prev_log_likelihood = -1e10
        run = 0

        # instantiate uncertainty model object for baseline scenario
        self.model_dict['uncertainty'] = tb_model.SimpleTbModel(self.fixed_parameters, self.inputs, 0)

        # get initial set of parameters
        new_param_list = []
        for i in range(len(self.param_ranges_unc)):
            new_param_list.append(self.param_ranges_unc[i]['start'])
            params = new_param_list

        # until a sufficient number of parameters are accepted
        while n_accepted < self.uncertainty_accepted_runs:

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
                    bound_low, bound_high \
                        = self.param_ranges_unc[i]['lower_bound'], self.param_ranges_unc[i]['upper_bound']

                    # normalise value and find log of PDF from appropriate distribution
                    if self.param_ranges_unc[i]['distribution'] == 'beta':
                        prior_log_likelihood += beta.logpdf((param_val - bound_low) / (bound_high - bound_low), 2., 2.)
                    elif self.param_ranges_unc[i]['distribution'] == 'uniform':
                        prior_log_likelihood += numpy.log(1. / (bound_high - bound_low))

                # calculate posterior
                posterior_log_likelihood = 0.

                incidence_result \
                    = self.epi_outputs['uncertainty']['incidence'][self.find_time_index(target_incidence[2],
                                                                                        'uncertainty')]

                posterior_log_likelihood += norm.logpdf(incidence_result, target_incidence[0], target_incidence[1])

                # determine acceptance
                log_likelihood = prior_log_likelihood + posterior_log_likelihood
                accepted = numpy.random.binomial(n=1, p=min(1., numpy.exp(log_likelihood - prev_log_likelihood)))

                if bool(accepted):
                    n_accepted += 1

                    # update likelihood and parameter set for next run
                    prev_log_likelihood = log_likelihood
                    params = new_param_list
                    self.accepted_indices += [run]

                run += 1

                print('run')
                print(run)
                print('accepted')
                print(accepted)
                print('incidence')
                print(incidence_result)
                for i in range(len(self.param_ranges_unc)):
                    print(self.param_ranges_unc[i]['name'])
                    print(new_param_list[i])
                print('\n')

            new_param_list = self.update_params(params)

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
            param_dict[self.param_ranges_unc[i]['name']] = params[i]
        return param_dict

    def get_fitting_data(self, ):
        """
        Define the characteristics (mean and standard deviation) of the normal distribution for model outputs
        (incidence, mortality).

        Returns:
            normal_char: Dictionary with keys outputs and values dictionaries. Sub-dictionaries have keys years
                and values lists, with first element of list means and second standard deviations.
        """

        # dictionary storing the characteristics of the normal distributions
        normal_char = {}
        for output_dict in self.inputs.outputs_unc:
            normal_char[output_dict['key']] = {}

            # incidence
            if output_dict['key'] == 'incidence':
                for year in self.inputs.data_to_fit[output_dict['key']].keys():
                    low = self.inputs.data_to_fit['incidence_low'][year]
                    high = self.inputs.data_to_fit['incidence_high'][year]
                    sd = output_dict['width_multiplier'] * (high - low) / (2. * 1.96)
                    mu = (high + low) / 2.
                    normal_char[output_dict['key']][year] = [mu, sd]

        return normal_char

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
            search_width = self.param_ranges_unc[i]['search_width']
            random = -100.

            # search for new parameters
            while random < self.param_ranges_unc[i]['lower_bound'] or random > self.param_ranges_unc[i]['upper_bound']:
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

    def find_time_index(self, time, scenario):
        """
        Method to find first time point in times list at or after a certain specified time point.

        Args:
            time: Float for the time point of interest.
        """

        return [i for i, j in enumerate(self.model_dict[scenario].times) if j >= time][0] - 1
        raise ValueError('Time not found')




