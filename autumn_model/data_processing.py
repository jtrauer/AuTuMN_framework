
import spreadsheet
import curve
import copy


class Inputs:
    """
    Object responsible for processing both user inputs and data read from the input spreadsheets through the
    spreadsheet module.
    """

    def __init__(self, country, scenarios_to_run, fixed_parameters, time_variant_parameters,
                 scenario_implementation=[]):

        self.country = country
        self.scenarios = scenarios_to_run
        self.fixed_parameters = fixed_parameters
        self.time_variant_parameters = time_variant_parameters
        self.scenario_implementation = scenario_implementation
        self.original_data = None
        self.derived_data = {}
        self.scaleup_data = {}
        self.scaleup_fns = {}

    #####################
    ### Master method ###
    #####################

    def read_and_load_data(self):
        """
        Master method of this object, calling all sub-methods to read and process data and define model structure.
        """

        # read all required external data
        print('Reading Excel sheets with input data.\n')
        self.original_data = spreadsheet.read_input_data_xls(self.find_keys_of_sheets_to_read(), self.country)

        # work through all remaining parameter processing methods
        self.process_parameters()
        self.process_case_detection()
        self.extract_time_variants()
        self.create_scenarios()
        self.find_scaleup_functions()

    ###########################
    ### Second-tier methods ###
    ###########################

    def find_keys_of_sheets_to_read(self):
        """
        Find keys of spreadsheets to read. Pretty simplistic at this stage, but can get more complicated as
        further sheets are added.
        """

        return ['tb']

    def process_parameters(self):
        """
        Convert raw parameters into interpreted parameters.
        """

        # find late (non-infectious) treatment period from total and early treatment periods
        self.fixed_parameters['time_late_treatment'] \
            = self.fixed_parameters['time_treatment'] - self.fixed_parameters['time_early_treatment']

        # convert treatment progress rates from proportions into rates
        derived_parameters = {}
        for param in self.fixed_parameters:
            if 'program_prop' in param and '_infect' in param:
                derived_parameters[param.replace('_prop', '_rate')] \
                    = self.fixed_parameters[param] / self.fixed_parameters['time_early_treatment']
            elif 'program_prop' in param and '_noninfect' in param:
                derived_parameters[param.replace('_prop', '_rate')] \
                    = self.fixed_parameters[param] / self.fixed_parameters['time_late_treatment']
        self.fixed_parameters.update(derived_parameters)

    def process_case_detection(self):
        """
        Extract case detection rates from data read in from Global TB Report and then add a zero at the start.
        """

        self.time_variant_parameters['program_prop_detect'] \
            = {i: j / 1e2 for i, j in self.original_data['tb']['c_cdr'].items()}
        self.time_variant_parameters['program_prop_detect'][1950] = 0.

    def extract_time_variants(self):
        """
        Extract the data for each individual time-variant parameter to a dictionary for each scenario and for each
        scaling parameter.
        """

        for scenario in self.scenarios:
            self.scaleup_data[scenario] = {}
            for parameter in self.time_variant_parameters:
                self.scaleup_data[scenario][parameter] = copy.copy(self.time_variant_parameters[parameter])

    def create_scenarios(self):
        """
        Use scenario dictionary to add additional values to scale-up data functions.
        """

        for scenario in self.scenarios[1:]:
            self.scaleup_data[scenario][self.scenario_implementation[scenario]['intervention']].update(
                {self.scenario_implementation[scenario]['year']: self.scenario_implementation[scenario]['coverage']})

    def find_scaleup_functions(self):
        """
        Create functions for scaling parameters, using the curve module.
        """

        for scenario in self.scenarios:
            self.scaleup_fns[scenario] = {}
            for time_variant_parameter in ['program_prop_detect', 'prop_vaccination']:
                self.scaleup_fns[scenario][time_variant_parameter] \
                    = curve.function_creator(self.scaleup_data[scenario][time_variant_parameter])

