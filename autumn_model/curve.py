
import numpy

def sinusiodal_scaleup(time, x_start, y_start, duration, magnitude):
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


def function_creator(data):
    def scaleup_function(time):
        (keys, values) = zip(*data.iteritems())
        if time < keys[0]:
            return 0.
        elif keys[-1] <= time:
            return values[-1]
        else:
            for k in range(len(keys)):
                if keys[k] <= time < keys[k + 1]:
                    return sinusiodal_scaleup(
                        time, keys[k], values[k], float(keys[k + 1] - keys[k]), values[k + 1] - values[k])

    return scaleup_function