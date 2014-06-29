from math import pi

############# helper functions ###############

def _isiterable(obj, exclude=basestring):
    return hasattr(obj, '__iter__') and not isinstance(obj, exclude)

def celcius2kelvin(temp):
    return temp + 273.15

############ definitions

_GASES = {
    'air':      {'k': 1.4,  'R': 286.9},
    'helium':   {'k': 1.66, 'R': 2077},
    'hydrogen': {'k': 1.41, 'R': 4124},
    'nitrogen': {'k': 1.4,  'R': 296.8, 'formula': 'N2'},
}

_LIQUIDS = {
}

_SOLIDS = {
}

############ acoustic conversions ###############

def freq2wavelen(freq, speed='air@20'):
    """
    calculate the wavelength of a given frequency based
    on the soundspeed given.

    freq (Hz)  : the frequency to be converted
    speed (m/s): the speed of sound for the given medium
                 You can pass the medium directly, the speed
                 will be calculated.


    Returns --> the wavelength (lambda) in m

    NB: use soundspeed(medium, temp) to calculate the speed of sound
        for a specific medium.
    """
    speed = soundspeed(speed)
    wavelength = speed / freq
    return wavelength

def wavelen2freq(wavelength, speed='air@20'):
    """
    calculate the frequency corresponding to a given wavelength

    wavelength (m) : the wavelength to be converted
    speed (m/s)    : the speed of sound for the given medium
                     You can pass the medium directly, the speed
                     will be calculated.

    Returns --> the frequency in Hz

    NB: use soundspeed(medium, temp) to calculate the speed of sound
        for a specific medium.
    """
    speed = soundspeed(speed)
    freq = speed / wavelength
    return freq

def distance2delay(distance, speed='air@20'):
    """
    calculate the delay in seconds for a sound to travel the given distance
    at the indicated speed of sound.

    distance (m): distance between source and receptor
    speed       : speed of sound (m/s) or medium

    Returns --> delay in seconds
    """
    speed = soundspeed(speed)
    time = distance / speed
    return time

def delay2distance(delay, speed='air@20'):
    """
    calculate the distance necessary for a sound to arrive with a
    given `delay` at the indicated `speed` of sound
    """
    speed = soundspeed(speed)
    distance = delay * speed
    return distance

def speed_of_sound_gas(k, R, T):
    """
    k : ratio of specific heats
    R : gas constant
    T : absolute temperature

    from http://www.engineeringtoolbox.com/speed-sound-d_519.html
    """
    return (k * R * T) ** 0.5

def _speed_of_sound_hooks_law(E, p):
    return (E / p) ** 0.5

def _medium_to_function(medium):
    medium = medium.lower()
    if medium in _GASES:
        medium_properties = _GASES.get(medium)
        return lambda t: speed_of_sound_gas(medium_properties['k'], medium_properties['R'], celcius2kelvin(t))
    return None

SOUNDSPEED_SUPPORTED_MEDIA = _GASES.keys() + _LIQUIDS.keys() + _SOLIDS.keys()

def _parse_medium(medium, temp=20):
    if isinstance(medium, basestring):
        if '@' in medium:
            medium, temp = medium.split('@')
            temp = float(temp)
    elif isinstance(medium, (tuple, list)):
        medium, temp = medium
    return medium, temp

def soundspeed(medium='air', temp=20):
    """
    return the speed of sound for the given temperature and medium

    SOUNDSPEED_SUPPORTED_MEDIA holds a list of valid media.

    Temperature only makes sense for gases.

    Formats:
        for usability, all these function calls mean the same:

        soundspeed('air', 20)
        soundspeed(('air', 20))
        soundspeed('air@20')

    """
    try:
        speed = float(medium)
        return speed
    except (ValueError, TypeError):
        pass
    medium, temp = _parse_medium(medium, temp)
    if _isiterable(medium):
        return soundspeed(*medium)
    func = _medium_to_function(medium)
    if func is None:
        raise ValueError("medium not supported, see speed_of_sound_supported_media for more info")
    return func(temp)

def phaseshift(frequency, distance, medium='air@20'):
    """
    calculate the phase shift in radians of a signal with a given
    frequency after distance

    NB: to calculate the phase-shift after a given time delay,
    convert it with

    >>> delay2distance(delay, medium)

    phase = lamba * dt = 2pi * freq * dt

    where lambda : wave_length
          dt     : time_difference

    via: http://www.sengpielaudio.com/calculator-timedelayphase.htm
    """
    timedelay = distance2delay(distance, medium)
    phase_shift = 2 * pi * frequency * timedelay
    return phase_shift