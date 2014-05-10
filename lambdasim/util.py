import itertools

def window(iterable, size=2, step=1):
    """
    iterate over subseqs of iterable
    
    Example
    =======
    
    >>> seq = range(6)
    
    >>> list(window(seq, 3, 1))
    [(0, 1, 2), (1, 2, 3), (2, 3, 4), (3, 4, 5)]
    
    >>> list(window(seq, 3, 2))
    [(0, 1, 2), (2, 3, 4)]
    
    # the same as pairwise
    >>> assert list(window(range(5), 2, 1)) == [(0, 1), (1, 2), (2, 3), (3, 4)]
    """
    iterators = itertools.tee(iterable, size)
    for skip_steps, itr in enumerate(iterators):
        for ignored in itertools.islice(itr, skip_steps):
            pass
    window_itr = itertools.izip(*iterators)
    if step != 1:
        window_itr = itertools.islice(window_itr, 0, 99999999, step)
    return window_itr

def soundpressure_to_soundlevel(soundpressure, p0=0.00002):
    """
    convert soundpressure in Pascal to sound level in dB (dBSPL)

    Lp(dBSPL) = 20 * log10(p/p0)

    p0: threshold of hearing, 0.00002 Pa (20uPa)
    """
    return 20 * _math.log10(soundpressure/p0)

def soundlevel_to_soundpressure(soundlevel, p0=0.00002):
    """
    convert sound-level in dB to sound-pressure in Pascal

    p = p0 * e^(1/20*Lp*log10(10))

    p0: threshold of hearing, 0.00002 Pa (20uPa) 
    """
    return p0 * _math.exp(1/20*soundlevel*_math.log(10))

L2p = soundlevel_to_soundpressure
p2L = soundpressure_to_soundlevel

def detect_lambda():
    pass
