from __future__ import absolute_import, print_function, division

import os
import math
import numpy
import warnings
from collections import namedtuple

# 3rd party
import shapelib
import sndfileio

# local
from .conversions import *
from .io import *
from .filtering import *
from .source import *
from .config import config
from .util import *
from . import snd
from .errors import *
from .plotting import plot_array


Sample = namedtuple("Sample", "frames sr idx")


class Wall(object):
    def __init__(self, envmatrix, angmatrix=None, bounds=None):
        self.envmat = envmatrix
        self.angmat = angmatrix
        if bounds:
            self.x0, self.y0, self.x1, self.y1 = bounds


#################################################
#
#                  SIMULATION
#
#################################################



class Simulation(object):
    """
    Defines a lambda simulation
    """
    def __init__(self, samplerate=None, size=None, duration=1, resolution=None, c=343, rho=RHO_AIR):
        """
        samplerate [Hz] : the audio samplerate of the simulation. This sets the time/space
                          resolution of the grid. See `resolution`
        size [(m,m)]    : size of the simulation in meters (x, y), or None to define it later
        duration [s]    : duration of the simulation
        resolution [m]  : alternative to samplerate, set the resolution of the grid
                          by setting the distance of one step
        c [m/s]         : the sound propagation sample_define
        rho [kg/m3]     : the density of air

        Examples
        ========

        Simulation(samplerate=44100, size=(10, 6)) # 10x6 m

        s = Simulation(size=(10,6))
        s.nodesize = 0.001    # 1 mm nodesize

        s = Simulation(samplerate=44100).set_size_in_pixels(

        """
        if samplerate is not None and resolution is not None:
            raise ValueError("The resolution of the simulation should be set either through `samplerate` or `resolution`, but not both")
        if samplerate is None and resolution is None:
            samplerate = 44100
        elif resolution is not None:
            samplerate = nodesize2samplerate(resolution, c)
        self._samplerate = samplerate

        if size is not None:
            xpix, ypix = coord2pix(size[0], size[1], samplerate, c)
        else:
            xpix, ypix = 100, 100

        if xpix % 2 != 0:
            xpix += 1
        if ypix % 2 != 0:
            ypix += 1

        self._xsize = xpix
        self._ysize = ypix

        self.c = c
        self.set_duration(duration)
        self.rho = rho
        self.nfilters = 0
        self._lastwav = None
        self._lastsimfile = None
        self._rcefile = None
        self._autoadd = config['autoadd']
        self.sources = []
        self.samples = []
        self.filters = {}
        self.envmat = None
        self.angmat = None
        self.envgeom = None
        self.receivers = []
        self.reset()
        lambdabin = detect_lambda()
        if lambdabin is None:
            warnings.warn("Lambda has not been detected in your system")


    def reset(self):
        """
        Regenerates the internal state (environment matrix, angle matrix, etc.)
        The dimensions of the simulation (size, duration, etc.) remain the same
        All sources and walls are removed.
        """
        assert (self._ysize % 2) == (self._xsize % 2) == 0
        self.sources = []
        self.samples = []
        self.filters = {}
        self.envmat = numpy.zeros((self._ysize, self._xsize), dtype=float)
        self.angmat = numpy.ones((self._ysize, self._xsize), dtype=float) * 400.
        self.envgeom = shapelib.box(0, 0, 0, 0)
        self.receivers = []

    def __str__(self):
        if self.nodesize < 1:
            nodesize_str = "%d mm" % (self.nodesize * 1000)
        elif self.nodesize < 10:
            frac, i = math.modf(self.nodesize)
            nodesize_str = "%d m %d cm" % (int(i), int(frac*100))
        else:
            nodesize_str = "%.3f" % self.nodesize
        msg = (
            "Simulation",
            "----------",
            "nodesize    : %s" % nodesize_str,
            "samplerate  : %d Hz" % self.samplerate,
            "C           : %.1f m/s" % self.c,
            "rho         : %.3f kg/m3" % self.rho,
            "size        : %.2f m x %.2f m" % (self.size),
            "matrix      : %d pixels x %d pixels" % self.size_in_pixels,
            "num. sources: %d" % len(self.sources),
            "num. recs.  : %d" % len(self.receivers),
            "num. samples: %d" % len(self.samples),
        )
        return "\n".join(msg)

    def __repr__(self):
        x, y = self.size
        return "Simulation(samplerate={sr:.1f}, size=({x:.1f}, {y:.1f}), duration={dur:.3f})".format(
            sr=self.samplerate, x=x, y=y, dur=self.duration)

    @property
    def samplerate(self):
        return self._samplerate

    @samplerate.setter
    def samplerate(self, sr):
        self._samplerate = sr

    @property 
    def nodesize(self):
        """
        The size of each node, in m
        """
        return nodesize(self._samplerate, self.c)

    @nodesize.setter
    def nodesize(self, ns):
        self._samplerate = nodesize2samplerate(ns, self.c)

    @property
    def distancestep(self):
        """The distance between two nodes"""
        return distancestep(self.samplerate, self.c)

    @property
    def timestep(self):
        """The time resolution of the simulation"""
        return timestep(self.samplerate, self.c)

    @property
    def size_in_pixels(self):
        """ the size of the simulation in pixels """
        return (self._xsize, self._ysize)

    def set_size_in_pixels(self, x, y):
        """ set the size of the simulation in pixels """
        xold, yold = self.size_in_pixels
        if x % 2 != 0:
            x += 1
            warnings.warn("Setting x to an even number. x=%d" % x)
        if y % 2 != 0:
            y += 1
            warnings.warn("Setting y to an even number. y=%d" % y)
        self._xsize = x
        self._ysize = y
        env, src, ang = self.envmat, self.sourcesmat, self.angmat
        self.reset()
        y0, x0 = env.shape
        ymin, xmin = min(y, y0), min(x, x0)
        if xold != x or yold != y:
            if xold < x or yold < y:
                warnings.warn("The simulation space will be reduced. Anything outside the new boundaries will be cropped")
            self.envmat = numpy.zeros((y, x), dtype=float)
            self.envmat[:ymin,:xmin] = env[:ymin,:xmin]
        if ang is not None:
            self.angmat = numpy.ones((y, x), dtype=float) * 400.
            self.angmat[:ymin,:xmin] = ang[:ymin,:xmin]
        if self.sources:
            numoldsources = len(self.sources)
            self.sources = [source for source in self.sources if source.xpix < x and source.ypix < y]
            if numoldsources > len(self.sources):
                warnings.warn("The simulation space has been reduced."
                              "Sources were defined outside of the new boundaries and have been removed")
        return self

    @property
    def size(self):
        """ the real size of the simulation in meters """
        return pix2coord(self._xsize, self._ysize, self.samplerate, self.c)

    def set_size(self, x, y):
        """
        set the size of the simulation in meters.
        """
        X, Y = coord2pix(x, y, self.samplerate, self.c)
        return self.set_size_in_pixels(X, Y)

    def calculate_deadnode_geom(self, wallgeom, eps=1e-6):
        g00 = shapelib.tight_envelope(wallgeom)
        x, y = self.size
        dead = shapelib.box(0, 0, x + eps, y + eps).difference(g00)
        return dead

    def set_size_from_geom(self, geom, margin=0.1, deadnodes=False):
        _, _, maxx, maxy = geom.bounds
        maxx += margin
        maxy += margin
        self.set_size(maxx, maxy)
        if deadnodes:
            dead = self.calculate_deadnode_geom(geom, eps=0.5)
            self.dead_from_geometry(dead)
        return self

    @property
    def duration(self):
        """ 
        duration of this simulation in seconds. Always a multiple of timestep

        SEE ALSO: timestep, set_duration
         """
        return self.timestep * self.steps

    def set_duration(self, duration_in_secs):
        """ 
        set the duration of the simulation in seconds. 
        The resulting duration will be rounded down to the nearest 
        multiple of timestep

        SEE ALSO: timestep, duration
        """
        steps = int(duration_in_secs / self.timestep)
        self.steps = steps
        return self

    #def gridbetween(self, x0, x1):
    #    L = self.nodesize
    #    X0, _ = snaptogrid(x0, 0, self.samplerate)
    #    X1, _ = snaptogrid(x1, 0, self.samplerate)
    #    xs = numpy.arange(X0, X1, L)
    #    return xs

    def snaptogrid(self, n):
        n, _ = self.coord2pix(n, n)
        n1, _ = self.pix2coord(n, n)
        return n1

    def make_receiver(self, x, y, label=""):
        """
        Create a receiver to sense pressure at position (x, y)

        x, y : coord of the receiver
        label: a string identifying the receiver. This is useful
               (and necessary) when using multiple receivers
        """
        px, py = map(int, self.coord2pix(x, y))
        self.envmat[py, px] = -2
        self.angmat[py, px] = 0
        self.receivers.append((x, y, label))
        X, Y = self.size_in_pixels
        # sort them the way they will be interpreted later
        self.receivers.sort(key=lambda rec: rec[1] * Y + rec[0])


    #def geometry_to_pixels(self, geom):
    #    """
    #    convert a shapely geometry to the pixels that form it
    #    """
    #    x0, y0, x1, y1 = geom.bounds
    #    px0, py0 = self.coord2pix(x0, y0)
    #    px1, py1 = self.coord2pix(x1, y1)
    #    for py in range(max(0, py0 - 1), py1 + 1):
    #        for px in range(max(0, px0 - 1), px1 + 1):
    #            rx0, ry0 = self.pix2coord(px, py)
    #            rx1, ry1 = rx0 + self.nodesize, ry0 + self.nodesize
    #            pix = Polygon([(rx0, ry0), (rx1, ry0), (rx1, ry1), (rx0, ry1)])
    #            if geom.intersects(pix):
    #                yield (px, py)

    def source_point(self, x, y, kind, amp=1, freq=440, phase=0, idx=None):
        """
        Add a point source to the simulation

        x, y : coordinates of the point, in m
        kind : sin, square, deltapulse, expdecay, hannsin,
               vel-sin, vel-square, vel-deltapulse, vel-hannsin,
               whitenoise, pinknoise, sample

        amp  : the amplitude of the source.
        freq : the freq of the source, when applicable.
        idx  : when applicable ("sample" source), the index of the sample.
        """
        xpix, ypix = coord2pix(x, y, self.samplerate)
        source = Source(xpix, ypix, kind, amp, freq, phase, sampleidx=idx)
        if self._autoadd:
            self._add_source(source) 
        return source

    def source_from_geometry(self, geom, kind, amp, freq=None, phase=None, sampleidx=None):
        """
        geom: a shapely geometry
        kind: sin, square, deltapulse, expdecay, hannsin,
              vel-sin, vel-square, vel-deltapulse, vel-hannsin,
              whitenoise, pinknoise, sample
        amp: the amplitude, in Pa
        freq: the freq. of the source, if applicable
        phase: the phase of the source, if applicable
        sampleidx: the index of the sample, for a sample source (30)
        """

        mat = self._rasterize(geom)
        return self.source_from_array(mat, kind=kind, amp=amp, freq=freq, phase=phase, sampleidx=sampleidx)

    def source_from_array(self, mat, kind, amp, freq=None, phase=None, sampleidx=None, adapt='crop'):
        """
        Create a source where the array is non zero

        kind : sin, square, deltapulse, expdecay, hannsin,
               vel-sin, vel-square, vel-deltapulse, vel-hannsin,
               whitenoise, pinknoise, sample

        adapt: defines the behaviour when the matrix is bigger than
               the simulation space
               'crop' --> crops the array to the simulation space
               'grow' --> grows the simulation space to the dimensions of the array

        Example
        =======

        # Create a source from the red traces of an image
        mask = img_to_mask(imgpath, color=(255,0,0), distance=0.2)
        sim.source_from_array(mask, 'sin', 0.2, 220)
        """
        sources = SourceList()
        if phase is None:
            phase = 0
        Y, X = numpy.nonzero(mat)
        for x, y in zip(X, Y):
            sources.append(Source(x, y, kind, amp, freq, phase, sampleidx=sampleidx))
        if self._autoadd:
            for source in sources:
                self._add_source(source)
        return sources


    # def ang_from_geometry(self, geom, angle=None):
    #     """
    #     generate a matrix where for each pixel where the
    #     geom is defined the given (or, when implemented, the calculated)
    #     angle is given, or 400 (any value outside 360) to indicate an empty pixel
    #     """
    #
    #     mat = self._geomgrid.rasterize(geom)
    #     if angle is not None:
    #         mat_fore = mat * angle
    #         mat[mat == 0] *= 400
    #         mat += mat_fore
    #     else:
    #         raise ValueError("calculating angles out of the geometry is not implemented YET")
    #     if self._autoadd:
    #         self._add_ang(mat)
    #     return mat

    def _add_ang(self, mat):
        if self.angmat is None:
            self.angmat = numpy.ones_like(self.envmat, dtype=float) * 400.0
        i = mat != 400
        self.angmat[i] = mat[i]

    def filter_from_geometry(self, geom, numid):
        """
        geom  = a shapely geometry
        numid = the .numid of the filter
        """
        assert any(f.numid == numid for f_name, f in self.filters.iteritems())
        assert numid >= 2
        mat = self._rasterize(geom)
        I = mat > 0
        mat *= numid
        if self._autoadd:
            self.envmat[I] = mat[I]
        return mat

    def filter_define(self, name, filtertype, freq, param, dbgain=0):
        """
        To apply a filter you first define it here, then use it
        in the simulation with something like `filter_from_geometry`,
        where you pass the .numid of the Filter defined here.
        Once a filter is defined, you can retrieve it later by calling
        `simulation.filters.get('filtername')`

        name: a name or a number to identify this filter.
              If a number if given, it will be used a id
              and will be the same when the simulation is 
              loaded from a .sim file 
              If a name is given, a number will be assigned.
              The numerical-id can be found as .numid
        filtertype : one of 'lpf', 'hpf', 'bpf'
        freq       : cut-off freq
        param      : Depending on filter, acts as Q, bandwidth, or shelf slope.
        dbgain     : gain  in db for peaking/shelving filters (defaults to 0).
        """
        assert name not in self.filters
        self.nfilters += 1
        if isinstance(name, int):
            numid = name
        else:
            numid = self.nfilters + 1    
        assert numid >= 2
        f = Filter.butter(filtertype=filtertype, freq=freq, samplerate=samplerate, numid=numid)
        self.filters[name] = f
        return f

    def sample_define(self, source):
        """
        To use a sample as source, first you define it here, then use the index
        to create a 'sample' source

        source: either the path to a soundfile, or a tuple (frames, samplerate)
        """
        if isinstance(source, basestring):
            sample = sndfileio.sndread(source)
            frames, sr = sample.samples, sample.sr
            channels = sndfileio.numchannels(frames)
            if channels > 1:
                print("The sample has %d channels. Taking only the first one" % channels)
                frames = sndfileio.getchannel(frames, 0)
        else:
            frames, sr = source
        duration = len(frames) / sr
        if duration > self.duration:
            print("The sample is %0.2f seconds long. Shortening to match the duration of the simulation" % sample.duration)
            sample = sample[:self.duration]
        if sr != self.samplerate:
            print("The samplerate of the sample (%d) is different from the samplerate of the simulation (%d)" % (sample.samplerate, self.samplerate))
            print("--> Resampling")
            frames = sndfileio.resample(frames, sr, self.samplerate)
        self.samples.append(frames)
        idx = len(self.samples) - 1
        return Sample(frames, sr, idx)
 
    def wall_from_geometry(self, geom, param=1, calculate_angle=False):
        """
        Creates a Wall from a geometry defined with `shapelib`

        geom : a geometry created with shapelib (or with Shapely)
        param: # TODO
        calculate_angle: calculate the angles from the geometry

        Returns --> a Wall
        """
        # TODO: use self._rasterize, using shapelib.rasterize
        mat = geom_rasterize(geom, self.size, self.size_in_pixels)
        mat *= param
        if calculate_angle:
            ang = self.angles_from_geom(geom)
        else:
            ang = None
        wall = Wall(mat, ang)
        if self._autoadd: 
            self._add_wall(wall, geom=geom)
        return wall

    def wall_from_array(self, A, adapt='crop'):
        """
        If A is smaller than the simulation space, it will be padded
        with 0s.
        adapt -> 'crop':  if A is smaller than the simulation matrix, it will be padded
                          if bigger, it will be cropped.
                 'grow':  if A is smaller than the simulation matrix, it will be padded
                          if bigger, the simulation space will grow accordingly
        """
        adapt_options = ('crop', 'grow')
        if adapt not in adapt_options:
            raise ValueError("adapt should be one of %s" % str(adapt_options))
        nx, ny = self.size_in_pixels
        ay, ax = A.shape
        if adapt == 'crop':
            print("adapt: crop")
            zeros = numpy.zeros((ny, nx))
            zeros[:ay, :ax] += A[:ny,:nx]
            wall = Wall(zeros)
        elif adapt == 'grow':
            print("adapt: grow")
            if ax > nx or ay > ny:
                X = max(ax, nx)
                Y = max(ay, ny)
                self.set_size_in_pixels(X, Y)
            return self.wall_from_array(A, adapt='crop')
        if self._autoadd:
            self._add_wall(wall)
        return wall

    def wall_from_image(self, path, adapt='crop', color=(255, 255, 255), dist=0.98):
        """Create a wall from an image.

        :param path: the path of the image defining a new wall
        :param adapt: 'crop' or 'grow', determines what happends when the image
                      is bigger than the simulation space
        :param color: the color of the wall (r, g, b)
        :param dist: max distance to the color
        :return: adds a wall to the simulation, returns it
        """
        mask = img_to_mask(path, color, dist)
        return self.wall_from_array(mask, adapt=adapt)

    def dead_from_geometry(self, geom):
        return self.wall_from_geometry(geom, param=DEADNODE)

    def coord2pix(self, x, y):
        return coord2pix(x, y, self.samplerate)

    def pix2coord(self, x, y):
        return pix2coord(x, y, self.samplerate)

    def plot_walls(self):
        """
        plot the walls of the simulation
        """
        X, Y = self.size
        plot_array(self.envmat, X, Y)

    def plot(self):
        colors = {
            'env':0.3,
            'src':1.0
        }
        mat = self.envmat * colors['env']
        mat += self.sourcesmat  * colors['src']
        X, Y = self.size
        plot_array(mat, X, Y, cmap='spectral')

    def _add_wall(self, wall, geom=None):
        if wall.envmat.shape != self.envmat.shape:
            wall.envmat = wall.envmat[:self.envmat.shape[0], :self.envmat.shape[1]]
            if wall.angmat:
                wall.angmat = wall.angmat[:self.envmat.shape[0], :self.envmat.shape[1]]

        i = wall.envmat > self.envmat
        self.envmat[i] = wall.envmat[i]
        deadnodes = (wall.envmat == DEADNODE)
        if deadnodes.any():
            i = deadnodes * (self.envmat == 0)  # only empty spaces can be declared as dead
            self.envmat[i] = DEADNODE
        if wall.angmat is not None:
            self.angmat[i] = wall.angmat[i]
        if geom is not None:
            self.envgeom = self.envgeom.union(geom)

    def _add_source(self, source):
        """ source can be an individual Source or a list of Sources 
           (sources are always point sources) """
        assert isinstance(source, Source)
        self.sources.append(source)
        
    @property
    def sourcesmat(self):
        """
        generate a 2D array representing the position of sources
        """
        out = numpy.zeros((self._ysize, self._xsize), dtype=float)
        for source in self.sources:
            out[source.ypix, source.xpix] = 1
        return out

    @property
    def simfile(self):
        return self._lastsimfile

    def write(self, outfile=None):
        """
        write .sim file
        """
        if outfile is None:
            if self._lastsimfile is None:
                warnings.warn("this Simulation has not been saved before. saving to a temporary file")
                import tempfile
                outfile = tempfile.mktemp(suffix='.sim', dir=os.getcwd())
            else:
                outfile = self._lastsimfile
                print("overwriting sim file: %s" % outfile)
        if not self.validate():
            raise SimValidationError("Error when validating this simulation")

        outfile = os.path.splitext(outfile)[0] + ".sim"
        f = open(outfile, 'w')
        f.write('LAMBDASIM200')
        f.write('DEF')
        # fwrite(simFile,[YSIZE XSIZE STEPS C L RHO],'double');
        header = numpy.array([self._ysize, self._xsize, self.steps, self.c, self.nodesize, self.rho], dtype=float)
        header.tofile(f)
        if self.envmat is not None:
            f.write('ENV')
            self.envmat.tofile(f)
        if self.angmat is not None:
            f.write('ANG')
            self.angmat.tofile(f)
        if self.filters:
            f.write('FLT')
            write_double(f, len(self.filters))
            for name, filterdef in self.filters.iteritems():
                filterdef.asarray().tofile(f)
        if self.samples:
            f.write('SMP')
            write_double(f, len(self.samples))
            for i, sample in enumerate(self.samples):
                write_double(f, i)  # filter ID
                write_double(f, sample.samplerate)
                write_double(f, len(sample.samples))
                sample.samples.astype(float).tofile(f)
        if not self.sources:
            raise ValueError("no sources defined, can't write sim file!")
        f.write('SRC')
        write_double(f, len(self.sources))
        for source in self.sources:
            a = source.asmatlab()
            a.tofile(f)
        self._lastsimfile = os.path.abspath(outfile)
        f.close()
        return self

    def _get_receivers_labels(self):
        labels = []
        unnamed = 0
        for rec in self.receivers:
            x, y, label = rec
            if not label:
                xpix, ypix = self.coord2pix(x, y)
                label = "recv{num}-{xpix}x{ypix}".format(num=unnamed, xpix=xpix, ypix=ypix)
                unnamed += 1
            labels.append(label)
        return labels

    def rce2wav(self, rcefile=None, resample=None, split=None):
        """
        Convert the already rendered `rce` to a wav file
        """
        if rcefile is None and self._lastsimfile is None:
            raise ValueError("this simulation has not been run to produce an rce file")
        if rcefile is None:
            rcefile = os.path.splitext(self._lastsimfile)[0] + ".rce"
        if not os.path.exists(rcefile):
            raise IOError("rce file not found")
        labels = self._get_receivers_labels()
        wavfile = rce2wav(rcefile, self.samplerate, resample=resample, split=split, splitsuffixes=labels)
        self._lastwav = wavfile
        return wavfile

    @property
    def wavfile(self):
        return self._lastwav

    def split_rendered_audio(self, path=None):
        """
        Split a rendered file into its channels with the labels given to 
        the receivers.
        Leave path unspecified to use the last rendered .wav (.wavfile)
        """
        if path is None:
            path = self.wavfile
        if not path or not os.path.exists(path):
            raise IOError(".wav file not found! cannot split")
        labels = self._get_receivers_labels()
        return snd.split_channels(self.wavfile, labels)

    def render_receivers(self, duration=None, resample=None):
        """
        Render the receivers to a .rce file.

        :param duration: in secods, it overrides the duration of the simulation
        :param resample: in Hz, overrides the duration of the simulation
        :return: a Future holding the samples. If resample was given, they will be resampled
                 to the given samplerate. Otherwise, the samples will have the samplerate of
                 this simulation
        """
        if not self.simfile:
            self.write()
        rcefile = os.path.splitext(self.simfile)[0] + '.rce'
        self._rcefile = rcefile
        simfile = os.path.abspath(self.simfile)
        args = ['-rce', '-exit', '-file', simfile]
        if duration is not None:
            iterations = duration2steps(duration, self.samplerate, self.c)
            args.extend(['-iterations', iterations])

        def render_and_load():
            proc = call_lambda(args)
            proc.wait()
            if os.path.exists(rcefile):
                samples = rce2array(rcefile)
                if resample:
                    samples = sndfileio.resample(samples, self.samplerate, resample)
                return samples
            else:
                return None
        return scheduler.call_now(render_and_load)

    def render_wav(self, duration=None, resample=None, split=False):
        """
        Render the simulation as a .wav

        :param duration: if given, overrides the duration of the simulation (in seconds)
        :param resample: if goven, overrides the samplerate of the simulation
        :param split: if True, for each receiver a soundfile will be generated
                      if False, a soundfile with as many channels as receivers will be generated
        :return: a Future holding the samples as a multichannel numpy array
        """
        labels = self._get_receivers_labels() if split else None
        samplerate = resample if resample is not None else self.samplerate
        def func():
            samples = self.render_receivers(duration=duration, resample=resample).result()
            if samples is not None:
                outfile = os.path.splitext(self.simfile)[0] + '.wav'
                return snd.sndwrite(samples, samplerate, outfile, labels)
            else:
                raise RuntimeError("Receivers did not render correctly")
        return scheduler.call_now(func)

    def render_video(self, walls=True, duration=None, contrast=50, quality=100, fps=None, cmap=None):
        """
        Render the simulation as video (async)

        walls (bool): whether to render the walls or not
        duration (s): render a part of the simulation, or
                      None to render the whole duration
        contrast (0-100): sets the contrast
        quality (0-100) : the video quality
        fps : the fps of the video, None for default
        cmap (int): index of the colormap
                    0 - grayscale
                    1 - fire
                    2 - temperature

        Returns --> Future(subprocess)
        """
        if not self.simfile:
            self.write()
        simfile = os.path.abspath(self.simfile)
        if " " in simfile:
            simfile = '"%s"' % simfile
        args = ['-avi', '-exit', '-file', simfile, '-contrast', contrast, '-quality', quality]
        if fps is not None:
            args.extend(['-avifps', str(fps)])
        if walls:
            args.append("-walls")
        if duration is not None:
            iterations = duration2steps(duration, self.samplerate, self.c)
            args.extend(['-iterations', iterations])
        if cmap is not None:
            args.extend(['-colormap', cmap])
        subp = call_lambda(args)
        return scheduler.wrap_subproc(subp)

    def _rasterize(self, geom):
        return geom_rasterize(geom, self.size, self.size_in_pixels)

    def angles_from_geom(self, geom):
        from math import degrees
        sr = self.samplerate
        edge = shapelib.edge(geom)
        mat = self._rasterize(edge)
        xpixs, ypixs = numpy.nonzero(mat)
        angles = numpy.ones_like(mat) * 400.
        maxy, maxx = angles.shape
        for xpix, ypix in zip(xpixs, ypixs):
            if xpix < maxx and ypix < maxy:
                x, y = pix2coord(xpix, ypix, sr)
                angle = degrees(shapelib.angle_at(geom, (x, y)))
                angles[ypix, xpix] = angle
        return angles

    def opensim(self, vis=True, walls=True, contrast=50, cmap=None, pipe=None, fps=None, skip=None):
        """
        Open this Simulation in Lambda
        """
        if self.simfile is None:
            print("Simulation needs to be written")
            self.write()
        return open_sim_in_lambda(self.simfile, vis=vis, walls=walls, contrast=contrast, cmap=cmap, pipe=pipe, fps=fps, skip=skip)

    def validate(self):
        """
        Check that this Simulation makes sense and, if problems
        are found, tries to solve them.

        1) Check for collisions between sources and walls
        2) Check that the sample sources have corresponding
           samples defined
        """
        # check dimensions
        xpix, ypix = self.size_in_pixels
        if xpix % 2 or ypix % 2:
            print("Pixel size should be even!")
            return False
        # remove sources at walls
        envmat = self.envmat
        def wall_at_source(src):
            return envmat[src.ypix, src.xpix] > 0
        nsources = len(self.sources)
        self.sources = [source for source in self.sources if not wall_at_source(source)]
        if not self.sources:
            print("No sources defined!")
            return False
        nsources2 = len(self.sources)
        if nsources2 < nsources:
            print("%d sources found at walls were removed" % (nsources - nsources2))
        for source in self.sources:
            if source.kind == 'sample':
                if not (0 <= source.freq < len(self.samples)):
                    raise SimValidationError("A sampled source was defined but no corresponding samples were found")
        return True

    def export_walls(self, outfile, color=None, background=(0, 0, 0)):
        """Export the walls defined in this simulation.

        outfile: the path to export the walls to. The format will be
                 determined by the extension.
                 Formats allowed: PNG
        """
        ext = os.path.splitext(outfile)[1].lower()
        if ext == '.png':
            if color is None:
                color = config['pngcolors']['wall']
            mat = (abs(self.envmat) > 0).astype(int)
            if isinstance(color, tuple):
                if background == color:
                    background = (255-color[0], 255-color[1], 255-color[2])
                colormap = lambda x: interpolate_color(x, background, color)
            else:
                colormap = color
            png_save(mat, outfile, colormap)


def readsim(path):
    """
    read a .sim file

    Returns a Simulation

    #TODO: support SMP chunk
    """
    def readassert(fileobj, s):
        out = fileobj.read(len(s))
        return out == s
    def toint(n):
        intn = int(n)
        assert n == intn
        return intn

    with open(path, 'rb') as f:
        readassert(f, 'LAMBDASIM200')
        readassert(f, 'DEF')
        ysize, xsize, steps, c, l, rho = read_doubles(f, 6)
        ysize, xsize, steps = map(toint, (ysize, xsize, steps))
        samplerate = nodesize2samplerate(l, c)
        envmat, angmat, sources = None, None, None
        def read_SRC(f):
            numsources = toint(read_double(f))
            srcs = []
            for i in range(int(numsources)):
                y, x, sourcetype, amp, freq, phase = read_doubles(f, 6)
                srcs.append( Source(x - 1, y - 1, sourcetype, amp, freq, phase) )
            return srcs
        while True:
            chunkheader = f.read(3)
            if chunkheader == 'ENV':
                envmat = read_numpy_mat(f, ysize, xsize)
            elif chunkheader == 'ANG':
                angmat = read_numpy_mat(f, ysize, xsize)
            elif chunkheader == 'SRC':
                sources = read_SRC(f)
                # SRC is always the last chunk, so exit the loop
                break
            else:
                raise ValueError("chunk %s not supported!" % chunkheader)
    duration = simulationduration(samplerate, steps, c)
    sim = Simulation(samplerate=samplerate, duration=duration, c=c, rho=rho)
    sim.set_size_in_pixels(xsize, ysize)
    sim.envmat  = envmat
    sim.angmat  = angmat
    sim.sources = sources
    assert len(sources) > 0
    return sim


def simwrite(outfile, samplerate, dur_ms, c, rho, pixsize, sources,
             envmat=None, angmat=None, filters=None, samples=None):
    """
    write .sim file

    outfile: the path of the .sim to write
    samplerate: the samplerate of the simulation
    dur_ms : duration in ms
    c      : propagation speed
    rho    : density
    pixsize: the size of the canvas in pixels (x, y)
    sources: a list of Sources
    envmat : a matrix representing the env., or None
    angmat : a matrix representing the angles, or None
    filters: a list of Filters, or None
    samples: a list of Samples, or None
    """
    outfile = os.path.splitext(outfile)[0] + ".sim"
    f = open(outfile, 'w')
    f.write('LAMBDASIM200')
    f.write('DEF')
    xsize, ysize = pixsize
    steps = duration2steps(dur_ms/1000., samplerate, c)
    ns = nodesize(samplerate, c)
    # write header
    numpy.array([ysize, xsize, steps, c, ns, rho], dtype=float).tofile(f)
    if envmat is not None:
        assert envmat.shape == (ysize, xsize)
        assert envmat.shape[0]%2 == 0
        assert envmat.shape[1]%2 == 0
        f.write('ENV')
        envmat.tofile(f)
    if angmat is not None:
        assert angmat.shape == (ysize, xsize)
        if envmat:
            assert angmat.shape == envmat.shape
        assert angmat.shape[0] % 2 == 0
        assert angmat.shape[1] % 2 == 0
        f.write('ANG')
        angmat.tofile(f)
    if filters:
        f.write('FLT')
        write_double(f, len(filters))
        for filt in filters:
            filt.asarray().tofile(f)
    if samples:
        f.write('SMP')
        write_double(f, len(samples))
        for i, sample in enumerate(samples):
            write_double(f, i)  # filter ID
            write_double(f, sample.samplerate)
            write_double(f, len(sample.samples))
            sample.samples.astype(float).tofile(f)
    if not sources:
        raise ValueError("no sources defined, can't write sim file!")
    f.write('SRC')
    write_double(f, len(sources))
    for source in sources:
        a = source.asmatlab()
        a.tofile(f)
    f.close()


def deriv2angle(deriv):
    """
    convert a derivative to an angle in 360 form. 
    0 is N, 90 is E
    """
    tan_alpha = deriv
    alpha = math.degrees(math.atan(tan_alpha))
    alpha = (360 - alpha + 90) % 360
    return alpha


def rce2array(rcefile):
    """
    Read a .rce (raw pressure) file, returns a numpy array of shape (numframes, numchannels)
    for multichannel files, or shape = (numframes,) for mono files

    :param rcefile: the path to a .rce file
    """
    f = open(rcefile, 'rb')
    numch = int(read_double(f))
    raw = numpy.fromfile(f, dtype=float)
    assert len(raw) % numch == 0
    if numch > 0:
        raw.shape = (len(raw) / numch, numch)
    return raw


def samples2wav(samples, samplerate, outfile, resample=None, split=False, splitsuffixes=None):
    def name_with_suffix(origname, suffix, ext=None):
        base, origext = os.path.splitext(origname)
        if ext is None:
            ext = origext
        return "%s-%s.%s" % (base, suffix, ext)
    if resample is not None:
        samples = sndfileio.resample(samples, samplerate, resample)
        samplerate = resample
    numch = samples.shape[1] if len(samples.shape) > 1 else 1

    if not split or numch == 1:
        sndfileio.sndwrite(samples, samplerate, outfile)
        return outfile
    else:
        outfiles = []
        if splitsuffixes is None:
            splitsuffixes = [str(i) for i in range(numch)]
        else:
            splitsuffixes = splitsuffixes[:numch]
        for i, suffix in enumerate(splitsuffixes):
            channelfile = name_with_suffix(outfile, suffix=suffix)
            outfiles.append(channelfile)
            sndfileio.sndwrite(samples[:,i], samplerate, channelfile)
        return outfiles


def rce2wav(rcefile, samplerate, resample=None, split=False, splitsuffixes=None):
    """
    a rce file is a raw, float64 file containing pressure level at each frame.

    rcefile    : path to the .rce file
    samplerate : samplerate of the simulation (sim.samplerate)
    resample   : new samplerate or None to keep current samplerate
    split      : create individual soundfiles for each receiver or
                 a multichannel soundfile in the case of multiple
                 receivers
    splitsuffixes : if multiple receivers are present and `split`
                    is True, use these suffixes (the number must
                    match the number of receivers)
                    Otherwise, each soundfile will be suffixed
                    with the index of each receiver

    FORMAT:
    1 double: number of sources
    Each frame then contains the data for each source, interleaved

    NOTE:
    an .rce file is simply raw data. To load it: numpy.fromfile(path, dtype=float).
    The samplerate is not saved with the data, but it is the same used by the
    simulation which created it.

    Returns --> the names of the outfile(s) written
    """
    samples = rce2array(rcefile)
    outfile = os.path.splitext(rcefile)[0] + ".wav"
    return samples2wav(samples, samplerate=samplerate, outfile=outfile, resample=resample, split=split, splitsuffixes=splitsuffixes)






