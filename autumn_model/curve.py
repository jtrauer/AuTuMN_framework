
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