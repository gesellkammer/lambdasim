import numpy
from numbers import Number

source_kinds = {
    'sin': 1,
    'square': 2,
    'deltapulse': 3,
    'expdecay': 4,
    'hannsin': 5,
    'vel-sin': 6,
    'vel-square': 7,
    'vel-deltapulse': 8,
    'vel-expdecay': 9,
    'vel-hannsin': 10,
    'whitenoise': 20,
    'pinknoise': 21,
    'sample': 30
}

source_kinds_inv = {v:k for k, v in source_kinds.items()}

def source_name_to_index(kind):
    return source_kinds.get(kind)

def source_index_to_name(n):
    return source_kinds_inv.get(n)


class Source(object):
    def __init__(self, xpix, ypix, kind='hannsin', amp=1, freq=440, phase=0, sampleidx=None):
        """
        A source is a point source, defined in pixel space

        amp      : amplitude (sound pressure in Pa)
        freq     : freq of source (unused in whitenoise, pinknoise)
        phase    : phase of source (unused for whitenoise, pinknoise, sample)
        sampleidx: idx of a Sample definition, used only for the 'sample' kind

        NB: convert dB (sound-pressure-level) to Pa with dB2Pa
        """
        self.xpix = xpix
        self.ypix = ypix
        if isinstance(kind, Number):
            assert kind in source_kinds_inv
            kind_index = kind
            kind_name = source_index_to_name(kind_index)
        else:
            kind_name = kind.lower()
            kind_index = source_name_to_index(kind_name)
            assert kind in source_kinds

        self.kind = kind_index
        if kind_name == 'sample':
            freq = sampleidx
        self.freq = freq
        self.amp = amp
        self.phase = phase

    @property
    def kind_name(self):
        return source_index_to_name(self.kind)

    def __repr__(self):
        return "source kind:%s pos:(x=%d y=%d) amp:%f freq:%f phase:%d" % (
            self.kind_name, self.xpix, self.ypix, self.amp, self.freq, int(self.phase)
        )

    def asarray(self):
        return numpy.array(
            [self.ypix, self.xpix, self.kind, self.amp, self.freq, self.phase],
            dtype=float
        )

    def asmatlab(self):
        mat = self.asarray()
        mat[0] += 1
        mat[1] += 1
        return mat

    @property
    def deadnodes(self):
        return self.envmat == _DEADNODE


class SourceList(list):
    def __str__(self):
        numpoints = len(self)
        freqs = set(source.freq for source in self)
        minx = miny = float('inf')
        maxx = maxy = float('-inf')
        miny = min(source.ypix for source in self)
        maxy = max(source.ypix for source in self)
        minx = min(source.xpix for source in self)
        maxx = max(source.xpix for source in self)
        coords = "({minx}, {miny}) - ({maxx}, {maxy})".format(**locals())

        s = ("SourceList\n"
             "----------\n"
             "numpoints: {numpoints}\n"
             "freqs    : {freqs}\n"
             "coords   : {coords}\n"
        ).format(**locals())

        return s

