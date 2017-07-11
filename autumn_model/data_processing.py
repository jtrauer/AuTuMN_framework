
import spreadsheet
import tool_kit
import numpy
import matplotlib.pyplot as plt
import curve


class Inputs:

    def __init__(self, country, scenarios_to_run):

        # GUI inputs
        self.country = country
        self.scenarios = scenarios_to_run

        # parameter structures
        self.original_data = None
        self.derived_data = {}
        self.time_variants = {}
        self.model_constants = {}
        self.scaleup_data = {}
        self.scaleup_fns = {}
        self.param_ranges_unc = []
        self.data_to_fit = {}
        # for incidence for ex, width of normal posterior relative to CI width in data
        self.outputs_unc = [{'key': 'incidence', 'posterior_width': None, 'width_multiplier': 2.}]

        # model structure
        self.available_strains = ['_ds', '_mdr', '_xdr']
        self.available_organs = ['_smearpos', '_smearneg', '_extrapul']
        self.agegroups = None
        self.riskgroups = []
        self.mixing = {}
        self.compartment_types \
            = ['susceptible_fully', 'susceptible_vac', 'susceptible_treated', 'latent_early', 'latent_late', 'active',
               'detect', 'missed', 'treatment_infect', 'treatment_noninfect']

        # interventions
        self.irrelevant_time_variants = []
        self.relevant_interventions = {}
        self.interventions_to_cost = {}
        self.intervention_startdates = {}
        self.potential_interventions_to_cost \
            = ['vaccination', 'xpert', 'treatment_support', 'smearacf', 'xpertacf', 'ipt_age0to5', 'ipt_age5to15',
               'decentralisation', 'improve_dst', 'bulgaria_improve_dst', 'intensive_screening', 'ipt_age15up',
               'ngo_activities', 'opendoors_activities']
        self.freeze_times = {}

        # miscellaneous
        self.mode = 'uncertainty'
        self.plot_count = 0
        self.treatment_outcome_types = []
        self.include_relapse_in_ds_outcomes = True

    #####################
    ### Master method ###
    #####################

    def read_and_load_data(self):
        """
        Master method of this object, calling all sub-methods to read and process data and define model structure.
        """

        # read all required data
        print('Reading Excel sheets with input data.\n')
        self.original_data = spreadsheet.read_input_data_xls(self.find_keys_of_sheets_to_read(), self.country)

        self.process_case_detection()

        for scenario in self.scenarios:
            self.scaleup_data[scenario]['prop_vaccination'] = {1921: 0., 1980: .8, 2015: .85}
            if scenario == 1:
                self.scaleup_data[scenario]['prop_vaccination'][2017] = 1.
            elif scenario == 2:
                self.scaleup_data[scenario]['program_prop_detect'][2017] = .9
            self.scaleup_fns[scenario] = {}
            for time_variant_parameter in ['program_prop_detect', 'prop_vaccination']:
                self.scaleup_fns[scenario][time_variant_parameter] \
                    = curve.function_creator(self.scaleup_data[scenario][time_variant_parameter])

        # x_values = numpy.linspace(1900., 2050., 10001)
        # result = [self.scaleup_fns['program_rate_detect'](x) for x in x_values]
        # plt.plot(x_values, result)
        # plt.show()
        # print()

        # # process constant parameters
        # self.process_model_constants()
        #
        # # process time-variant parameters
        # self.process_time_variants()
        #
        # # define model structure
        # self.define_model_structure()
        #
        # # find parameters that require processing
        # self.find_additional_parameters()
        #
        # # classify interventions as to whether they apply and are to be costed
        # self.classify_interventions()
        #
        # # calculate time-variant functions
        # self.find_scaleup_functions()
        #
        # # define compartmental structure
        # self.define_compartment_structure()
        #
        # # uncertainty-related analysis
        # self.process_uncertainty_parameters()
        #
        # # optimisation-related methods
        # self.find_intervention_startdates()  # currently sitting with intervention classification methods, though
        #
        # # perform checks (undeveloped still)
        # self.checks()

    def process_case_detection(self):

        for scenario in self.scenarios:
            self.scaleup_data[scenario] = {}
            self.scaleup_data[scenario]['program_prop_detect'] \
                = {i: j / 1e2 for i, j in self.original_data['tb']['c_cdr'].items()}
            self.scaleup_data[scenario]['program_prop_detect'][1950] = 0.

    #############################################
    ### Constant parameter processing methods ###
    #############################################

    # populate with first round of unprocessed parameters
    def process_model_constants(self):
        """
        Master method to call methods for processing constant model parameters.
        """

        # note ordering to list of sheets to be worked through is important for hierarchical loading of constants
        sheets_with_constants = ['country_constants', 'default_constants']
        if self.gui_inputs['riskgroup_diabetes']: sheets_with_constants += ['diabetes']
        self.add_model_constant_defaults(sheets_with_constants)

        # add "by definition" hard-coded parameters
        self.add_universal_parameters()

    # derive further parameters
    def find_additional_parameters(self):
        """
        Find additional parameters.
        Includes methods that require the model structure to be defined,
        so that this can't be run with process_model_constants.
        """

        # find risk group-specific parameters
        if len(self.riskgroups) > 1: self.find_riskgroup_progressions()

        # calculate rates of progression to active disease or late latency
        self.find_latency_progression_rates()

        # find the time non-infectious on treatment from the total time on treatment and the time infectious
        self.find_noninfectious_period()

        # derive some basic parameters for IPT
        self.find_ipt_params()

    def find_latency_progression_rates(self):
        """
        Find early progression rates by age group and by risk group status - i.e. early progression to active TB and
        stabilisation into late latency.
        """

        for agegroup in self.agegroups:
            for riskgroup in self.riskgroups:

                prop_early = self.model_constants['tb_prop_early_progression' + riskgroup + agegroup]
                time_early = self.model_constants['tb_timeperiod_early_latent']

                # early progression rate is early progression proportion divided by early time period
                self.model_constants['tb_rate_early_progression' + riskgroup + agegroup] = prop_early / time_early

                # stabilisation rate is one minus early progression proportion divided by early time period
                self.model_constants['tb_rate_stabilise' + riskgroup + agegroup] = (1. - prop_early) / time_early

    def find_riskgroup_progressions(self):
        """
        Code to adjust the progression rates to active disease for various risk groups - so far diabetes and HIV.
        """

        # initialise dictionary of additional adjusted parameters to avoid dictionary changing size during iterations
        risk_adjusted_parameters = {}
        for riskgroup in self.riskgroups:
            for param in self.model_constants:

                # start from the assumption that parameter is not being adjusted
                whether_to_adjust = False

                # for age-stratified parameters
                if '_age' in param:

                    # find the age string, the lower and upper age limits and the parameter name without the age string
                    age_string, _ = tool_kit.find_string_from_starting_letters(param, '_age')
                    age_limits, _ = tool_kit.interrogate_age_string(age_string)
                    param_without_age = param[:-len(age_string)]

                    # diabetes progression rates only start from age groups with lower limit above the start age
                    # and apply to both early and late progression.
                    if riskgroup == '_diabetes' and '_progression' in param \
                            and age_limits[0] >= self.model_constants['riskgroup_startage' + riskgroup]:
                        whether_to_adjust = True

                    # HIV applies to all age groups, but only late progression
                    elif riskgroup == '_hiv' and '_late_progression' in param:
                        whether_to_adjust = True

                    # shouldn't apply this to the multiplier parameters or non-TB-specific parameters
                    if '_multiplier' in param or 'tb_' not in param:
                        whether_to_adjust = False

                    # now adjust the age-stratified parameter values
                    if whether_to_adjust:
                        risk_adjusted_parameters[param_without_age + riskgroup + age_string] \
                            = self.model_constants[param] \
                              * self.model_constants['riskgroup_multiplier' + riskgroup + '_progression']
                    elif '_progression' in param:
                        risk_adjusted_parameters[param_without_age + riskgroup + age_string] \
                            = self.model_constants[param]

                # parameters not stratified by age
                else:

                    # explanation as above
                    if riskgroup == '_diabetes' and '_progression' in param:
                        whether_to_adjust = True
                    elif riskgroup == '_hiv' and '_late_progression' in param:
                        whether_to_adjust = True
                    if '_multiplier' in param or 'tb_' not in param:
                        whether_to_adjust = False

                    # adjustment as above, except age string not included
                    if whether_to_adjust:
                        risk_adjusted_parameters[param + riskgroup] \
                            = self.model_constants[param] \
                              * self.model_constants['riskgroup_multiplier' + riskgroup + '_progression']
                    elif '_progression' in param:
                        risk_adjusted_parameters[param + riskgroup] \
                            = self.model_constants[param]

        self.model_constants.update(risk_adjusted_parameters)

    #########################################
    ### Methods to define model structure ###
    #########################################

    def define_model_structure(self):
        """
        Master method to define all aspects of model structure.
        """

        self.define_riskgroup_structure()

    def define_riskgroup_structure(self):
        """
        Work out the risk group stratification.
        """

        # create list of risk group names
        for item in self.gui_inputs:
            if item[:10] == 'riskgroup_':
                if self.gui_inputs[item] and 'riskgroup_prop' + item[9:] in self.time_variants:
                    self.riskgroups += [item[9:]]
                elif self.gui_inputs[item]:
                    self.add_comment_to_gui_window(
                        'Stratification requested for %s risk group, but proportions not specified'
                        % tool_kit.find_title_from_dictionary(item[9:]))

        # add the null group
        if len(self.riskgroups) == 0:
            self.riskgroups += ['']
            self.vary_force_infection_by_riskgroup = False
            self.add_comment_to_gui_window(
                'Heterogeneous mixing requested, but not implemented as no risk groups are present')
        else:
            self.riskgroups += ['_norisk']

        # ensure some starting proportion of births go to the risk group stratum if value not loaded earlier
        for riskgroup in self.riskgroups:
            if 'riskgroup_prop' + riskgroup not in self.model_constants:
                self.model_constants['riskgroup_prop' + riskgroup] = 0.

    #################################################
    ### Time variant parameter processing methods ###
    #################################################

    def process_time_variants(self):
        """
        Master method to perform all preparation and processing tasks for time-variant parameters.
        Does not perform the fitting of functions to the data, which is done later in find_scaleup_functions.
        Note that the order of call is important and can lead to errors if changed.
        """

        if 'country_programs' in self.original_data: self.time_variants.update(self.original_data['country_programs'])
        self.add_time_variant_defaults()  # add any necessary time-variants from defaults if not in country programs
        self.load_vacc_detect_time_variants()
        self.convert_percentages_to_proportions()
        self.tidy_time_variants()

    def add_time_variant_defaults(self):
        """
        Populates time-variant parameters with defaults if those values aren't found in the manually entered
        country-specific data.
        """

        for program_var in self.original_data['default_programs']:

            # if the key isn't in available for the country
            if program_var not in self.time_variants:
                self.time_variants[program_var] = self.original_data['default_programs'][program_var]

            # otherwise if it's there and load_data is requested in the country sheet, populate for the missing years
            else:
                for year in self.original_data['default_programs'][program_var]:
                    if year not in self.time_variants[program_var] \
                            and 'load_data' in self.original_data['country_programs'][program_var] \
                            and self.original_data['country_programs'][program_var]['load_data'] == u'yes':
                        self.time_variants[program_var][year] \
                            = self.original_data['default_programs'][program_var][year]

    def convert_percentages_to_proportions(self):
        """
        Converts time-variant dictionaries to proportions if they are loaded as percentages in their raw form.
        """

        for time_variant in self.time_variants.keys():
            if 'perc_' in time_variant:  # if a percentage
                perc_name = time_variant.replace('perc', 'prop')
                self.time_variants[perc_name] = {}
                for year in self.time_variants[time_variant]:
                    if type(year) == int or 'scenario' in year:  # to exclude load_data, smoothness, etc.
                        self.time_variants[perc_name][year] = self.time_variants[time_variant][year] / 1e2
                    else:
                        self.time_variants[perc_name][year] = self.time_variants[time_variant][year]

    def tidy_time_variants(self):
        """
        Final tidying of time-variants, as described in comments to each line of code below.
        """

        for time_variant in self.time_variants:

            # add zero at starting time for model run to all program proportions
            if 'program_prop' in time_variant or 'int_prop' in time_variant:
                self.time_variants[time_variant][int(self.model_constants['start_time'])] = 0.

            # remove the load_data keys, as they have been used and are now redundant
            self.time_variants[time_variant] \
                = tool_kit.remove_specific_key(self.time_variants[time_variant], 'load_data')

            # remove keys for which values are nan
            self.time_variants[time_variant] = tool_kit.remove_nans(self.time_variants[time_variant])

    def find_data_for_functions_or_params(self):

        """
        Method to load all the dictionaries to be used in generating scale-up functions to a single attribute of the
        class instance (to avoid creating heaps of functions for irrelevant programs).

        Creates: self.scaleup_data, a dictionary of the relevant scale-up data for creating scale-up functions in
            set_scaleup_functions within the model object. First tier of keys is the scenario to be run, next is the
            time variant parameter to be calculated.
        """

        for scenario in self.scenarios:
            self.scaleup_data[scenario] = {}

            # Find the programs that are relevant and load them to the scaleup_data attribute
            for time_variant in self.time_variants:
                if time_variant not in self.irrelevant_time_variants:
                    self.scaleup_data[scenario][str(time_variant)] = {}
                    for i in self.time_variants[time_variant]:
                        if i == 'scenario_' + str(scenario):
                            self.scaleup_data[scenario][str(time_variant)]['scenario'] = \
                                self.time_variants[time_variant][i]
                        elif type(i) == str:
                            if 'scenario_' not in i:
                                self.scaleup_data[scenario][str(time_variant)][i] = self.time_variants[time_variant][i]
                        else:
                            self.scaleup_data[scenario][str(time_variant)][i] = self.time_variants[time_variant][i]

    ###################################
    ### Uncertainty-related methods ###
    ###################################

    def process_uncertainty_parameters(self):
        """
        Master method to uncertainty processing, calling other relevant methods.
        """

        # specify the parameters to be used for uncertainty
        if self.gui_inputs['output_uncertainty']:
            self.find_uncertainty_distributions()
            self.get_data_to_fit()

    def find_uncertainty_distributions(self):
        """
        Populate a dictionary of uncertainty parameters from the inputs dictionary in a format that matches code for
        uncertainty.
        """

        for param in self.model_constants:
            if '_uncertainty' in param and type(self.model_constants[param]) == dict:
                self.param_ranges_unc += [{'key': param[:-12],
                                           'bounds': [self.model_constants[param]['lower'],
                                                      self.model_constants[param]['upper']],
                                           'distribution': 'uniform'}]

    def get_data_to_fit(self):
        """
        Extract data for model fitting. Choices currently hard-coded above.
        """

        # decide whether calibration or uncertainty analysis is being run
        if self.mode == 'calibration':
            var_to_iterate = self.calib_outputs
        elif self.mode == 'uncertainty':
            var_to_iterate = self.outputs_unc

        inc_conversion_dict = {'incidence': 'e_inc_100k',
                               'incidence_low': 'e_inc_100k_lo',
                               'incidence_high': 'e_inc_100k_hi'}
        mort_conversion_dict = {'mortality': 'e_mort_exc_tbhiv_100k',
                                'mortality_low': 'e_mort_exc_tbhiv_100k_lo',
                                'mortality_high': 'e_mort_exc_tbhiv_100k_hi'}

        # work through vars to be used and populate into the data fitting dictionary
        for output in var_to_iterate:
            if output['key'] == 'incidence':
                for key in inc_conversion_dict:
                    self.data_to_fit[key] = self.original_data['tb'][inc_conversion_dict[key]]
            elif output['key'] == 'mortality':
                for key in mort_conversion_dict:
                    self.data_to_fit[key] = self.original_data['tb'][mort_conversion_dict[key]]
            else:
                self.add_comment_to_gui_window(
                    'Warning: Calibrated output %s is not directly available from the data' % output['key'])


    def find_keys_of_sheets_to_read(self):
        """
        Find keys of spreadsheets to read. Pretty simplistic at this stage, but expected to get more complicated as
        other sheets (like diabetes) are added as optional.
        """

        return ['tb']
