
import numpy
import tool_kit

"""
Very simplistic curve fitting module. In reality, AuTuMN uses piecewise-composed polynomial spline functions to fit
time-variant parameters to data. However, this adds a number of levels of complexity, as the polynomial functions
frequently must have a restricted range. Therefore, a very simple sinusoidal curve-fitting module is presented here as
a basic example.
"""


def sinusoidal_scaleup(time, x_start, y_start, duration, magnitude):
    """
    Uses a cosine function to create a simple sinusional function scaling up over "magnitude" and period of time
    "duration" from starting value "y_start" at time "x_start".

    Args:
         time: Argument to be fed to the function
         x_start: Starting time
         y_start: Starting value
         duration: Period of time that scaling occurs over
         magnitude: Distance to increase
    """

    return (-numpy.cos(numpy.pi * ((time - x_start) / duration)) + 1) / 2 * magnitude + y_start


def function_creator(input_dict):
    """
    Creates a scaling function from the input data, returning the function as an object, on the assumption the parameter
    is scaling from an initial value of zero to a constant final value at the latest year of the input dictionary.
    Uses the previous function (sinusoidal_scaleup) in a piecewise manner to create the overall function.

    Args:
         input_dict: Dictionary with keys time in years (usually integers) and values the values of parameters to be fit
    """

    def scaleup_function(time):
        time_data, values = tool_kit.extract_ordered_list_from_dictionary(input_dict)
        if time < time_data[0]:
            return 0.
        elif time_data[-1] <= time:
            return values[-1]
        else:
            for k in range(len(time_data)):
                if time_data[k] <= time < time_data[k + 1]:
                    return sinusoidal_scaleup(time, time_data[k], values[k],
                                              float(time_data[k + 1] - time_data[k]), values[k + 1] - values[k])
    return scaleup_function
