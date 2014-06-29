import math

#######################################
#
#       DEFAULTS AND CONSTANTS
#
#######################################

C = 340  # m/s
RHO_AIR = 1.204
K = 0.0020589
DEADNODE = -999

_SQRT2 = math.sqrt(2)
_EPSILON = 1e-13

#####################################
#
#  Lambda
#
#####################################


def nodesize(samplerate, c=C):
    """
    The size in meters of one node (one "pixel" of the simulation)
    
    PARAMETERS
    
    samplerate: rate of the simulation
    c         : speed of propagation
    
    Returns --> size (m)

    It is important to know that tube length, sound speed and 
    sampling frequency are closely linked together in lambda.
    """
    return c * _SQRT2 / samplerate
    

def timestep(samplerate, c=C):
    """
    The realtime elapsed time in one iteration
    
    PARAMETERS

    samplerate : the samplerate of the simulation
    c          : speed of propagation

    RETURNS -> the time in seconds of one iteration

    NB: the total time of the simulation is timestep * steps
    """
    return K * nodesize(samplerate, c)


def distancestep(samplerate, c=C):
    """
    the distance travelled by a wave with speed of
    propagation c in one time step
    
    PARAMETERS
    samplerate : samplerate of the simulation
    c          : speed of propagation
    
    RETURNS
    distance (m) : distance travelled by a wave in one 
                   time-step of the simulation
    """
    return timestep(samplerate, c) * c


def simulationduration(samplerate, steps, c=C):
    """
    Calculate the duration of a simulation in seconds
    """
    return timestep(samplerate, c) * steps


def steps2duation(steps, samplerate, c=C):
    return simulationduration((samplerate, steps, c))


def duration2steps(duration, samplerate, c=C):
    """
    convert a duration in seconds to the number of
    steps in the simulation
    """
    ts = timestep(samplerate, c)
    return int(duration / ts)


def canvassize(xpixels, ypixels, samplerate, c=C):
    """
    Return the real-world dimension of the canvas in meters
    
    xpixels    : width of the canvas in pixels
    ypixels    : height of the canvas in pixels
    samplerate : rate of the simulation
    c          : speed of propagation
    
    RETURNS
    (x, y) in meters
    """
    L = nodesize(samplerate, c)
    x = L * xpixels
    y = L * ypixels
    return x, y
    

def coord2pix(x, y, samplerate, c=C):
    """
    Convert a coordinate in meters within the canvas to pixels
    
    x, y: position in meters
    
    RETURNS
    (xpixel, ypixel): the pixel index corresponding to the position given
    """
    L = float(nodesize(samplerate, c))
    xpix = int(float(x) / L)
    ypix = int(float(y) / L)
    return xpix, ypix
    

def pix2coord(x, y, samplerate, c=C):
    """
    Convert (x, y) to realworlds coordinates

    x, y: pixel coordinates

    RETURNS
    (x, y) in meters
    """
    L = float(nodesize(samplerate, c))
    x *= L
    y *= L
    return x, y


def snaptogrid(x, y, samplerate, c=C):
    """
    Snap the point (x, y) (in meters) to a grid determined by samplerate

    x, y: real world coordinates (in meters)

    RETURNS 

    (x, y) in meters, snapped to the nearest representable point
    """
    x2, y2 = coord2pix(x, y, samplerate)
    x3, y3 = pix2coord(x2, y2, samplerate)
    return x3, y3


def nodesize2samplerate(ns, c=C):
    """
    convert nodesize to samplerate

    ns: size of each node in meters

    RETURNS
    samplerate (in Hz)
    """
    samplerate = _SQRT2 * c / ns
    return samplerate


def samplerate2nodesize(sr, c=C):
    """
    sr: samplerate, in Hz

    RETURNS --> nodesize: in meters
    """
    ns = _SQRT2 * c / sr
    return ns