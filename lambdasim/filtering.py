import numpy
from sndfileio import dsp


class Filter(object):
    def __init__(self, b, a, numid=None):
        """
        numid = an integer id
        b     = the sequence of b coefficients
        a     = the sequence of a coefficients
        """
        self.numid = numid
        self.bb = b
        self.aa = a

    @classmethod
    def butter(cls, filtertype, freq, samplerate, order=5, numid=None):
        """
        filtertype: one of 'low', 'high', 'band'
        freq      : cutoff freq -- (low, high) for bandpass filter
        samplerate: samplerate of the simulation
        order     : order of the filter

        Returns --> a Filter
        """
        assert filtertype in ('low', 'high', 'band')
        b, a = dsp.filter_butter_coeffs(filtertype, freq, samplerate=samplerate, order=order)
        return cls(b, a, numid=numid)

    def __repr__(self):
        return "filter\nbb=%s \naa=%s" % (str(self.bb), str(self.aa))

    def asarray(self):
        """
        Format
        =======

        ID:       double
        NUMCOEFS: double
        AA:       list of doubles
        BB:       list of doubles

        :rtype : numpy.array
        """
        out = [self.numid, len(self.aa), len(self.bb)]
        out.extend(self.aa)
        out.extend(self.bb)
        return numpy.array(out, dtype=float)