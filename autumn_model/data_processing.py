
import spreadsheet
import curve


class Inputs:

    def __init__(self, country, scenarios_to_run, fixed_parameters):

        self.country = country
        self.scenarios = scenarios_to_run
        self.fixed_parameters = fixed_parameters

        # parameter structures
        self.original_data = None
        self.derived_data = {}
        self.time_variants = {}
        self.model_constants = {}
        self.scaleup_data = {}
        self.scaleup_fns = {}

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

        self.process_parameters()

        self.process_case_detection()

        for scenario in self.scenarios:
            self.scaleup_data[scenario]['prop_vaccination'] = {1921: 0., 1980: .8, 2015: .85}
            if scenario == 1:
                self.scaleup_data[1]['prop_vaccination'][2020] = .99
            if scenario == 2:
                self.scaleup_data[2]['program_prop_detect'][2017] = .9
            self.scaleup_fns[scenario] = {}
            for time_variant_parameter in ['program_prop_detect', 'prop_vaccination']:
                self.scaleup_fns[scenario][time_variant_parameter] \
                    = curve.function_creator(self.scaleup_data[scenario][time_variant_parameter])

    ###########################
    ### Second-tier methods ###
    ###########################

    def process_parameters(self):

        self.fixed_parameters['time_late_treatment'] \
            = self.fixed_parameters['time_treatment'] - self.fixed_parameters['time_early_treatment']

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

        for scenario in self.scenarios:
            self.scaleup_data[scenario] = {}
            self.scaleup_data[scenario]['program_prop_detect'] \
                = {i: j / 1e2 for i, j in self.original_data['tb']['c_cdr'].items()}
            self.scaleup_data[scenario]['program_prop_detect'][1950] = 0.

    def find_keys_of_sheets_to_read(self):
        """
        Find keys of spreadsheets to read. Pretty simplistic at this stage, but can get more complicated as
        further sheets are added.
        """

        return ['tb']
