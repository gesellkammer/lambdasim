import sndfileio
import os


def split_channels(sndfile, suffixes=None):
    """
    Split the channels in `sndfile`, resulting in mono files
    with the same format as `sndfile`.

    sndfile  : the multichannel soundfile to be split
    suffixes : if given, it must match the number of channels
               in `sndfile`. They will be used as suffixes

    Returns --> the names of the resulting files

    Example
    =======

    Given a multichannel file 'snd.wav' with 3 channels,

    split_channels(sndfile, ['L', 'R', 'C'])

    Will produce three files
    'snd-L.wav', 'snd-R.wav' and 'snd-C.wav'
    """
    samples = sndfileio.sndread(sndfile)
    return sndwrite(samples, sndfile, suffixes)

def sndwrite(samples, samplerate, outfile, split_suffixes=None):
    """
    Write samples to one or more audiofiles files

    :param samples: a numpy array of samples, possibly multichannel
    :param outfile: the outfile to save the samples to
    :param split_suffixes: if given and samples is multichannel,
                           a list of suffixes to attach to outfile
                           or True to assign numeric suffixes
    :return: the list of soundfiles generates
    """
    base, ext = os.path.splitext(outfile)
    numchannels = sndfileio.numchannels(samples)
    outfiles = []
    if numchannels == 1:
        sndfileio.sndwrite(samples, samplerate, outfile, encoding='flt32')
        outfiles.append(outfile)
    else:
        def generate_labels(numchannels, suffixes):
            if suffixes is None or len(suffixes) < numchannels:
                suffixes = ["%2d" % i for i in range(1, numchannels+1)]
            assert len(suffixes) == numchannels
            return suffixes
        labels = generate_labels(numchannels, split_suffixes)
        for chan in range(numchannels):
            s = samples[:,chan]
            outfile = "{base}-{suffix}{ext}".format(base=base, suffix=labels[chan], ext=ext)
            print("saving channel {chan} to {outfile}".format(chan=chan, outfile=outfile))
            sndfileio.sndwrite(s, samplerate, outfile, encoding='flt32')
            outfiles.append(outfile)
    return outfiles