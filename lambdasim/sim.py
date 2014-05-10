from __future__ import absolute_import, print_function, division

import os
import math
from collections import namedtuple
import numpy
import array
import struct
from numbers import Number
import warnings
from shapely.geometry import Polygon, LineString
from bpf4 import api as bpf, interp_linear, interp_expon
from .util import p2L, L2p
import shapelib

    

C = 340  # m/s
RHO_AIR = 1.204
SQRT2 = math.sqrt(2)
K = 0.0020589
EPSILON = 1e-13
AUTOADD = True
LAMBDAPATH = "/usr/local/bin/lambda"
DEADNODE = -999


def read_double(fileobj):
    return struct.unpack('d', fileobj.read(struct.calcsize('d')))[0]


def write_double(fileobj, value):
    fileobj.write(struct.pack('d', value))


def read_doubles(fileobj, n):
    out = array.array('d')
    out.read(fileobj, n)
    return out


def read_numpy_mat(fileobj, ysize, xsize):
    arr = read_doubles(fileobj, ysize * xsize)
    mat = numpy.frombuffer(arr, dtype=float)
    mat.shape = (ysize, xsize)
    return mat


def nodesize(samplerate, c=C):
    """
    the size in m(eters) of one node (pixel)
    
    PARAMETERS
    
    samplerate: rate of the simulation
    c: speed of propagation
    
    It is important to know that tube length, sound speed and 
    sampling frequency are closely linked together in lambda.
    """
    return c * SQRT2 / samplerate
    

def timestep(samplerate, c=C, k=K):
    """
    the realtime elapsed time in one iteration
    
    PARAMETERS
    samplerate : the samplerate of the simulation
    c          : speed of propagation
    k          : an internal constant

    RETURNS
    the time in seconds of one iteration

    the total time of the simulation is timestep * steps
    """
    return k * nodesize(samplerate, c)


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
    return timestep(samplerate, c) * steps
    

def canvassize(xpixels, ypixels, samplerate, c=C):
    """
    return the dimension of the canvas in meters
    
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
    return (x, y)
    

def coord2pix(x, y, samplerate, c=C):
    """
    convert a coordinate in meters within the canvas to pixels
    
    PARAMETERS
    x, y: position in meters
    
    RETURNS
    (xpixel, ypixel): the pixel index corresponding to the position given
    """
    L = float(nodesize(samplerate, c))
    xpix = int(float(x) / L)
    ypix = int(float(y) / L)
    return (xpix, ypix)
    

def pix2coord(x, y, samplerate, c=C):
    """
    convert (x, y) to realworlds coordinates

    x, y: pixel coordinates

    RETURNS
    (x, y) in meters
    """
    L = float(nodesize(samplerate, c))
    x = x * L
    y = y * L
    return (x, y)


def snaptogrid(x, y, samplerate, c=C):
    """
    snap the point (x, y) to a grid determined by samplerate

    PARAMETERS

    x, y: real world coordinates
    """
    x2, y2 = coord2pix(x, y, samplerate)
    x3, y3 = pix2coord(x2, y2, samplerate)
    return (x3, y3)


def nodesize2samplerate(nodesize, c=C):
    """
    convert nodesize to samplerate

    nodesize: size of each node in meters
    """
    samplerate = SQRT2 * c / nodesize
    return samplerate


def samplerate2nodesize(sr, c=C):
    """
    RETURNS
    nodesize: in meters
    """
    ns = SQRT2 * c / sr
    return ns


class Wall(object):
    def __init__(self, envmatrix, angmatrix=None, bounds=None):
        self.envmat = envmatrix
        self.angmat = angmatrix
        if bounds:
            self.x0, self.y0, self.x1, self.y1 = bounds


class Simulation(object):
    def __init__(self, samplerate=44100, steps=0, c=343, rho=RHO_AIR):
        self._samplerate = samplerate
        self._xsize = 400
        self._ysize = 400
        self.c = c
        self.steps = steps
        self.rho = rho
        self.nfilters = 0
        self._lastwav = None
        self._lastsimfile = None
        self.reset()

    def reset(self):
        self._cachedgrid = None
        self.sources = []
        self.samples = []
        self.filters = {}
        self.nfilters = 0
        self.envmat = numpy.zeros((self._ysize, self._xsize), dtype=float)
        self.angmat = numpy.ones((self._ysize, self._xsize), dtype=float) * 400.
        self.envgeom = shapelib.box(0, 0, 0, 0)
        self.receivers = []

    @staticmethod
    def fromfile(path):
        def readassert(fileobj, s):
            out = fileobj.read(len(s))
            return out == s
        def toint(n):
            intn = int(n)
            assert n == intn
            return intn
        f = open(path, 'rb')
        readassert(f, 'LAMBDASIM200')
        readassert(f, 'DEF')
        defchunk = array.array('d')
        defchunk.read(f, 6)
        ysize, xsize, steps, c, l, rho = defchunk
        ysize, xsize, steps = map(toint, [ysize, xsize, steps])
        samplerate = nodesize2samplerate(l, c)
        envmat = None
        angmat = None
        sources = []
        while True:
            chunkheader = f.read(3)
            if chunkheader == 'SRC':
                break
            elif chunkheader == 'ENV':
                envmat = read_numpy_mat(f, ysize, xsize)
            elif chunkheader == 'ANG':
                angmat = read_numpy_mat(f, ysize, xsize)
            else:
                f.close()
                raise ValueError("chunk %s not supported!" % chunkheader)
        # now parse sources
        numsources = read_double(f)
        # make sure it is an int
        assert int(numsources) == numsources 
        for i in range(int(numsources)):
            sourcearr = array.array('d')
            sourcearr.read(f, 6)
            source_y, source_x, source_type, source_amp, source_freq, source_phase = sourcearr
            source = Source(source_x - 1, source_y - 1, source_type, source_amp, source_freq, source_phase)
            sources.append(source)
        sim = Simulation(samplerate, steps, c, rho)
        sim.set_size_in_pixels(xsize, ysize)
        sim.envmat = envmat
        sim.angmat = angmat
        sim.sources = sources
        return sim

    @property
    def samplerate(self):
        return self._samplerate

    @samplerate.setter
    def samplerate(self, sr):
        self._samplerate = sr

    @property 
    def nodesize(self):
        return nodesize(self._samplerate, self.c)

    @nodesize.setter
    def nodesize(self, ns):
        self._samplerate = nodesize2samplerate(ns, self.c)

    @property
    def distancestep(self):
        return distancestep(self.samplerate, self.c)

    @property
    def timestep(self):
        return timestep(self.samplerate, self.c)

    def size_in_pixels(self):
        """ the size of the simulation in pixels """
        return (self._xsize, self._ysize)

    def set_size_in_pixels(self, x, y):
        """ set the size of the simulation in pixels """
        self._xsize = x
        self._ysize = y
        self.reset()
        return self

    def size(self):
        """ the real size of the simulation in meters """
        return pix2coord(self._xsize, self._ysize, self.samplerate, self.c)

    def set_size(self, x, y):
        """ set the size of the simulation in meters """
        self._xsize, self._ysize = coord2pix(x, y, self.samplerate, self.c)
        self.reset()
        return self

    def calculate_deadnode_geom(self, wallgeom, eps=1e-6):
        g00 = shapelib.tight_envelope(wallgeom)
        x, y = self.size()
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

    def gridbetween(self, x0, x1):
        L = self.nodesize
        X0, _ = snaptogrid(x0, 0, self.samplerate)
        X1, _ = snaptogrid(x1, 0, self.samplerate)
        xs = numpy.arange(X0, X1, L)
        return xs

    def snaptogrid(self, n):
        n, _ = self.coord2pix(n, n)
        n1, _ = self.pix2coord(n, n)
        return n1

    def set_receiver(self, x, y, label=""):
        px, py = map(int, self.coord2pix(x, y))
        self.envmat[py, px] = -2
        self.angmat[py, px] = 0
        self.receivers.append((x, y, label))
        X, Y = self.size_in_pixels()
        # sort them the way they will be interpreted later
        self.receivers.sort(key=lambda rec: rec[1] * Y + rec[0])  

    def geometry_to_pixels(self, geom):
        """
        convert a shapely geometry to the pixels that form it
        """
        x0, y0, x1, y1 = geom.bounds
        px0, py0 = self.coord2pix(x0, y0)
        px1, py1 = self.coord2pix(x1, y1)
        for py in range(max(0, py0 - 1), py1 + 1):
            for px in range(max(0, px0 - 1), px1 + 1):
                rx0, ry0 = self.pix2coord(px, py)
                rx1, ry1 = rx0 + self.nodesize, ry0 + self.nodesize
                pix = Polygon([(rx0, ry0), (rx1, ry0), (rx1, ry1), (rx0, ry1)])
                if geom.intersects(pix):
                    yield (px, py)

    def source_from_geometry(self, geom, kind, amp, freq, phase=0, sampleidx=None):
        sources = []
        mat = self._rasterize(geom)
        xpixs, ypixs = numpy.nomzero(mat)
        for xpix, ypix in zip(xpixs, ypixs):
            sources.append(Source(xpix, ypix, kind, amp, freq, phase, sampleidx=sampleidx))
        if AUTOADD:
            self._add_source(sources)
        return sources

    def ang_from_geometry(self, geom, angle=None):
        """
        generate a matrix where for each pixel where the 
        geom is defined the given (or, when implemented, the calculated) 
        angle is given, or 400 (any value outside 360) to indicate an empty pixel
        """
        mat = self._geomgrid.rasterize(geom)
        if angle is not None:
            mat_fore = mat * angle
            mat[mat == 0] *= 400
            mat += mat_fore
        else:
            raise ValueError("calculating angles out of the geometry is not implemented YET")
        if AUTOADD:
            self._add_ang(mat)
        return mat

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
        mat = self._geomgrid.rasterize(geom)
        I = mat > 0
        mat *= numid
        if AUTOADD:
            self.envmat[I] = mat[I]
        return mat

    def filter_define(self, name, filtertype, freq, param, dbgain=0):
        """
        name: a name or a number to identify this filter.
              If a number if given, it will be used a id
              and will be the same when the simulation is 
              loaded from a .sim file 
              If a name is given, a number will be assigned. T
              he numerical-id can be found as .numid
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
        f = Filter.biquad(numid, filtertype=filtertype, freq=freq, param=param, dbgain=dbgain, sr=self.samplerate, normalized=False)
        self.filters[name] = f
        return f

    def sample_define(self, samples_or_filename, samplerate=None):
        if isinstance(samples_or_filename, basestring):
            from e import audiosample
            sample = audiosample.Sample(samples_or_filename)
            if sample.channels > 1:
                print("The sample has %d channels. Taking only the first one" % sample.channels)
                sample = sample.get_channel(0)
        else:
            raise NotImplemented('only filenames accepted at this point')
        assert sample.channels == 1
        if sample.duration > self.duration():
            print("The sample is %0.2f seconds long. Shortening to match the duration of the simulation" % sample.duration)
            sample = sample[:self.duration()]
        if sample.samplerate != self.samplerate:
            print("The samplerate of the sample (%d) is different from the samplerate of the simulation (%d)" % (sample.samplerate, self.samplerate))
            print("resampling...")
            sample = sample.resample(self.samplerate)
        self.samples.append(sample)
        sample.idx = len(self.samples) - 1
        return sample
 
    @property
    def _geomgrid(self):
        if self._cachedgrid is not None:
            return self._cachedgrid
        x, y = self.size()
        self._cachedgrid = grid = shapelib.Grid(0, 0, x, y, self.nodesize, self.nodesize)
        return grid

    def wall_from_geometry(self, geom, param=1, calculate_angle=False):
        grid = self._geomgrid
        m = grid.rasterize(geom) * param
        if calculate_angle:
            ang = self.angles_from_geom(geom)
        else:
            ang = None
        wall = Wall(m, ang)
        if AUTOADD: 
            self._add_wall(wall, geom=geom)
        return wall

    def dead_from_geometry(self, geom):
        return self.wall_from_geometry(geom, param=DEADNODE)
    
    def wall_line(self, x0, y0, x1, y1, param=1):
        angle = deriv2angle((y1 - y0) / (x1 - x0))
        x0, y0 = self.coord2pix(x0, y0)
        x1, y1 = self.coord2pix(x1, y1)
        l = LineString([(x0, y0), (x1, y1)])
        envmat = numpy.zeros((self._ysize, self._xsize), dtype=float)
        angmat = numpy.zeros((self._ysize, self._xsize), dtype=float)
        xs = numpy.arange(x0, x1 + 1, dtype=int)
        ys = numpy.zeros_like(xs)
        for i, x in enumerate(xs):
            intersection = shapelib.line_at_x(l, x)
            if intersection:
                ys[i] = intersection.y
        envmat[ys, xs] = param
        angmat[ys, xs] = angle
        wall = Wall(envmat, angmat, bounds=(x0, y0, x1, y1))
        if AUTOADD:
            self._add_wall(wall)
        return wall

    def coord2pix(self, x, y):
        return coord2pix(x, y, self.samplerate)

    def pix2coord(self, x, y):
        return pix2coord(x, y, self.samplerate)

    def torealspace(self, pix):
        return pix2coord(pix, 0, self.samplerate)[0]

    def topix(self, m):
        return coord2pix(m, 0, self.samplerate)[0]

    def source_point(self, x, y, kind, amp=1, freq=440, phase=0, sampleidx=None):
        """
        kind: see Simulation.available_sources()
        amp:  the amplitude of the source.
        freq: the freq of the source, when applicable.
        sampleidx: when applicable ("sample" source), the index of the source.
        """
        xpix, ypix = coord2pix(x, y, self.samplerate)
        source = Source(xpix, ypix, kind, amp, freq, phase, sampleidx=sampleidx)
        if AUTOADD:
            self._add_source(source) 
        return source

    def source_line(self, x0, y0, x1, y1, kind, amp, freq, phase=0):
        x0, y0 = map(int, self.coord2pix(x0, y0))
        x1, y1 = map(int, self.coord2pix(x1, y1))
        l = LineString([[x0, y0], [x1, y1]])
        sources = []
        ys = range(max(0, min(y0, y1) - 1), max(y0, y1) + 2)
        xs = range(max(0, min(x0, x1) - 1), max(x0, x1) + 2)
        for y in ys:
            for x in xs:
                pix = Polygon([[x, y], [x + 1, y], [x + 1, y + 1], [x, y + 1]])
                if l.intersects(pix):
                    sources.append(Source(x, y, kind, amp, freq, phase))
        if AUTOADD:
            self._add_source(sources)
        return sources


    def _add_wall(self, wall, geom=None):
        if wall.envmat.shape != self.envmat.shape:
            wall.envmat = wall.envmat[:self.envmat.shape[0], :self.envmat.shape[1]]
            if wall.angmat:
                wall.angmat = wall.angmat[:self.envmat.shape[0], :self.envmat.shape[1]]

        i = wall.envmat > self.envmat
        self.envmat[i] = wall.envmat[i]
        deadnodes = wall.envmat == DEADNODE
        if deadnodes.any():
            i = deadnodes * (self.envmat == 0)  # only empty spaces can be declared as dead
            self.envmat[i] = DEADNODE
        if wall.angmat is not None:
            self.angmat[i] = wall.angmat[i]
        if geom is not None:
            self.envgeom = self.envgeom.union(geom)

    def add_wall(self, wall):
        if AUTOADD:
            warnings.warn("AUTOADD is on. wall was added when created, bypassing.")
            return
        self._add_wall(wall)

    def _add_source(self, source):
        """ source can be an individual Source or a list of Sources (sources are always point sources) """
        if isinstance(source, Source):
            self.sources.append(source)
        elif isinstance(source, (list, tuple)):
            sources = source
            for source in sources:
                self._add_source(source)

    def add_source(self, source):
        if AUTOADD:
            warnings.warn("AUTOADD is on. source was added when create, bypassing")
            return
        self._add_source(source)

    def sourcesmat(self):
        out = numpy.zeros((self._ysize, self._xsize), dtype=float)
        for source in self.sources:
            out[source.ypix, source.xpix] = 1
        return out

    @property
    def simfile(self):
        return self._lastsimfile

    def write(self, outfile=None):
        if outfile is None:
            if self._lastsimfile is None:
                warnings.warn("this Simulation has not been saved before. saving to a temporary file")
                import tempfile
                outfile = tempfile.mktemp(suffix='.sim', dir=os.getcwd())
            else:
                outfile = self._lastsimfile
                print("overwriting sim file: %s" % outfile)
        self.validate()
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
        self._lastsimfile = outfile
        f.close()
        return self

    def _get_receivers_labels(self):
        def getlabel(receiver):
            x, y, label = receiver
            return label
        labels = []
        unnamed = 0
        for rec in self.receivers:
            label = getlabel(rec)
            if not label:
                label = "noname%d" % unnamed
                unnamed += 1
            labels.append(label)
        return labels

    def rce2wav(self, rcefile=None, resample=None, split=None):
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

    def split_wav(self, path=None):
        """
        split a rendered file into its channels with the labels given to 
        the receivers. 
        leave path unspecified to use .wavfile
        """
        path = path if path is not None else self.wavfile
        if not path or not os.path.exists(path):
            raise IOError(".wav file not found! cannot split")
        labels = self._get_receivers_labels()
        audiosample.split_channels(self.wavfile, labels)

    def render_wav(self, resample=None, split=None, block=True):
        if not self.simfile:
            self.write()
        import subprocess
        args = "/Applications/lambda.app/Contents/MacOS/lambda -rce -exit -file".split()
        args.append('"%s"' % os.path.abspath(self.simfile))
        
        if block:
            cmd = " ".join(args)
            os.system(cmd)
            # subprocess.call(args)
            rcefile = os.path.splitext(self.simfile)[0] + '.rce'
            assert os.path.exists(rcefile)
            self._rcefile = rcefile
            return self.rce2wav(rcefile, resample=resample, split=split)
        else:
            import threading
            def func():
                return self.render_wav(resample=resample, split=split, block=True)
            t = threading.Thread(target=func)
            t.start()
            return t

    def _rasterize(self, geom):
        return self._geomgrid.rasterize(geom)

    def angles_from_geom(self, geom):
        from math import degrees
        samplerate = self.samplerate
        edge = shapelib.edge(geom)
        mat = self._rasterize(edge)
        xpixs, ypixs = numpy.nonzero(mat)
        angles = numpy.ones_like(mat) * 400.
        maxy, maxx = angles.shape
        for xpix, ypix in zip(xpixs, ypixs):
            if xpix < maxx and ypix < maxy:
                x, y = pix2coord(xpix, ypix, samplerate)
                angle = degrees(shapelib.angle_at(geom, (x, y)))
                angles[ypix, xpix] = angle
        return angles

    def opensim(self, new=False):
        if not self.simfile:
            self.write()
        options = []
        if new:
            options.apend('--new')
        optstr = " ".join(options)
        cmd = 'open %s -a lambda --args -file "%s"' % (
            optstr,
            os.path.abspath(self.simfile)
        )
        os.system(cmd)

    def validate(self):
        # remove sources at walls
        envmat = self.envmat
        def wall_at_source(source):
            return envmat[source.ypix, source.xpix] > 0
        nsources = len(self.sources)
        self.sources = [source for source in self.sources if not wall_at_source(source)]
        nsources2 = len(self.sources)
        if nsources2 < nsources:
            print("%d sources found at walls were removed" % (nsources - nsources2))
        for source in self.sources:
            if source.kind == 'sample':
                assert 0 <= source.freq < len(self.samples)
        return self

    def plot_geom(self, geom):
        x0 = y0 = 0
        x1, y1 = self.size()
        return plot_geom(geom, (x0, x1), (y0, y1))


def _aslist(seq):
    if isinstance(seq, list):
        return seq
    return list(seq)


class Filter(object):
    def __init__(self, numid, bb, aa):
        """
        numid = an integer id
        bb    = the sequence of b coefficients
        aa    = the sequence of a coefficients
        """
        self.numid = numid
        self.bb = _aslist(bb)
        self.aa = _aslist(aa)

    @classmethod
    def biquad(cls, numid, filtertype, freq, param, dbgain=0, sr=44100, normalized=False):
        from e import dsp
        coeffs = dsp.biquad_coeffs(filtertype, fc=freq, param=param, dbgain=dbgain, fs=sr, normalized=normalized)
        return cls(numid, (coeffs.b0, coeffs.b1, coeffs.b2), (coeffs.a0, coeffs.a1, coeffs.a2))

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
        """
        out = [self.numid, len(self.bb), len(self.aa)]
        out.extend(self.bb)
        out.extend(self.aa)
        return numpy.array(out, dtype=float)


class Source(object):
    def __init__(self, xpix, ypix, kind, amp, freq=440, phase=0, sampleidx=None):
        """
        amp: amplitude (sound pressure in Pa)
        """
        self.xpix = xpix
        self.ypix = ypix
        if isinstance(kind, Number):
            self.kind = kind
        else:
            self.kind = {
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
            }.get(kind, 5)
        if kind == 'sample':
            freq = sampleidx
        self.freq = freq
        self.amp = amp
        self.phase = phase

    def __repr__(self):
        return "source kind:%d pos:(y=%d x=%d) amp:%f freq:%f phase:%d" % (
            self.kind, self.ypix, self.xpix, self.amp, self.freq, int(self.phase)
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
    Read a .rce (raw pressure) file
    """
    f = open(rcefile, 'rb')
    numch = read_double(f)
    raw = numpy.fromfile(f, dtype=float)
    if numch > 0:
        raw.shape = (len(raw) / numch, numch)
    return raw


def rce2wav(rcefile, samplerate, resample=None, split=False, splitsuffixes=None):
    """
    a rce file is a raw, float64 file containing pressure level at each frame.

    FORMAT: 
    1 double: number of sources
    Each frame then contains the data for each source, interleaved

    rcefile    : path to the .rce file
    samplerate : samplerate of the simulation (sim.samplerate)
    resample   : INT (new samplerate). 
                 If given, file will be resampled to this sample-rate.

    NOTE:
    an .rce file is simply raw data. To load it: numpy.fromfile(path, dtype=float).
    The samplerate is not saved with the data, but it is the same used by the
    simulation.
    """
    def name_with_suffix(origname, ext='wav', suffix=""):
        if suffix:
            suffix = "-%s" % str(suffix)
        base = os.path.splitext(origname)[0]
        return "%s%s.%s" % (base, suffix, ext)
    raw = rce2array(rcefile)
    from e import audiosample
    if resample is not None:
        try:
            from scikits import samplerate
            raw = samplerate.resample(raw, resample / samplerate, 'sync_fast')
            samplerate = resample
        except ImportError:
            warnings.warn("resampling not available (install scikits.samplerate to enable). Keeping original samplerate")
    numch = raw.shape[1] if len(raw.shape) > 1 else 1
    if not split or numch == 1:
        out = name_with_suffix(rcefile)
        audiosample.Sample(raw, samplerate).write(out)
        return out
    else:
        outs = []
        if splitsuffixes is None:
            splitsuffixes = [str(i) for i in range(numch)]
        else:
            splitsuffixes = splitsuffixes[:numch]
        for i, suffix in enumerate(splitsuffixes):
            out = name_with_suffix(rcefile, suffix=suffix)
            outs.append(out)
            audiosample.Sample(raw[:, i], samplerate).write(out)
        return outs






