from __future__ import division, print_function
import math
import itertools
import os
import sys
import subprocess
import shapelib
import numpy


from .config import config
from . import scheduler

def window(iterable, size=2, step=1):
    """
    iterate over subseqs of iterable
    
    """
    iterators = itertools.tee(iterable, size)
    for skip_steps, itr in enumerate(iterators):
        for _ in itertools.islice(itr, skip_steps):
            pass
    window_itr = itertools.izip(*iterators)
    if step != 1:
        window_itr = itertools.islice(window_itr, 0, 99999999, step)
    return window_itr


def soundpressure_to_soundlevel(Pa, p0=0.00002):
    """
    convert soundpressure in Pascal to sound level in dB (dBSPL)

    Lp(dBSPL) = 20 * log10(p/p0)

    p0: threshold of hearing, 0.00002 Pa (20uPa)
    """
    return 20 * math.log10(Pa / p0)


def soundlevel_to_soundpressure(dB, p0=0.00002):
    """
    convert sound-level in dB to sound-pressure in Pascal

    p = p0 * e^(1/20*Lp*log10(10))

    p0: threshold of hearing, 0.00002 Pa (20uPa) 
    """
    return p0 * math.exp(1 / 20 * dB * _math.log(10))


dB2Pa = L2p = soundlevel_to_soundpressure
Pa2dB = p2L = soundpressure_to_soundlevel


def _findfile(f, possible_folders):
    """
    Returns one of the `possible_folders` where `f` is present,
    or None
    """
    for folder in possible_folders:
        if os.path.exists(os.path.join(folder, f)):
            return folder
    return None


def detect_lambda():
    """
    Returns the path of the Lambda binary or None if not found
    """
    lambdabin = None
    def check_which():
        """
        :rtype : str or None. The path of the file
        """
        try:
            path = subprocess.check_output(["which", "lambda"])
            if os.path.exists(path):
                return path
        except subprocess.CalledProcessError:
            pass

    if config['lambdabin'] is not None:
        return config['lambdabin']

    if sys.platform == 'darwin':
        lambdabin = check_which()
        if lambdabin is None:
            possible_folders = [
                "/Applications",
                os.path.expanduser("~/Applications")
            ]
            lambda_app = _findfile('Lambda.app', possible_folders)
            if lambda_app:
                lambdabin = os.path.join(lambda_app, 'Lambda.app', 'Contents', 'MacOS', 'Lambda')
                assert os.path.exists(lambdabin), (
                    "Found the lambda app (%s) but the lambda binary was not present" % lambda_app)

    elif sys.platform == 'linux2':
        lambdabin = check_which()
    elif sys.platform == 'window':
        lambdabin = None
    if lambdabin is not None:
        config['lambdabin'] = lambdabin
    return lambdabin


def geom_rasterize(geom, size_in_meters, size_in_pixels):
    """
    rasterize the geometry `geom`, defined in real world coords,
    to a matrix of pixels with size `size_in_pixels`

    If geom is a line (a linestring, a linering), you should consider
    applying a .buffer before rasterizing, to control the "linesize"

    Returns --> a numpy.dnarray of 2D, where shape=size_in_pixels
                The array is binary
    """
    x_meters, y_meters = size_in_meters
    x_pixels, y_pixels = size_in_pixels
    pixratio = x_pixels / x_meters
    a = shapelib.rasterize(geom, pixratio, xrange=(0, x_meters), yrange=(0, y_meters)).array
    ay, ax = a.shape
    assert (ax, ay) == size_in_pixels
    return a


def call_lambda(args, stdout=None, stderr=None, wrap=False):
    """
    Call Lambda with the given `args`, as a subprocess

    args: passed to the lambda binary
    stdout, stderr: None, 'pipe' or 'log'
    wrap: wrap the subprocess in a future, to be able to add
          done callbacks

    Returns --> a subprocess

    NB: to add a done_callback:

    def finished(future):
        print("lambda finished!")

    call_lambda([arg1, arg2, ...], wrap=True).add_done_callback(finished)
    """
    binary = detect_lambda()
    if binary is None:
        if sys.platform == 'darwin':
            msg = ("Could not find the 'lambda' binary."
                   "Make sure the Lambda.app was copied to your"
                   "/Applications folder. If you installed it in another"
                   "location, create a symbolic link to the binary inside:"
                   "$ ln -s /path/to/Lambda.app/Contexts/MacOS/Lambda /usr/local/bin"
                   "and make sure that the destination is in your PATH"
            )
        else:
            msg = ""
        print(msg)
        raise IOError("Could not find the 'lambda' binary")
    args = [str(arg) for arg in args]
    cmds = [binary] + list(args)
    print("Calling Lamda as: %s" % str(cmds))
    if stdout == 'pipe':
        stdout = subprocess.PIPE
    elif stdout == 'log':
        stdout = open("lambda.log", "w")
    if stderr == 'pipe':
        stdout = subprocess.PIPE
    elif stderr == 'log':
        stdout = open("lambda-error.log", "w")
    if wrap:
        return scheduler.subproc_call(cmds, stdout=stdout, stderr=stderr)
    else:
        return subprocess.Popen(cmds, stdout=stdout, stderr=stderr)


def open_sim_in_lambda(simfile, vis=False, walls=True, contrast=None, cmap=None, fps=None, skip=None, pipe=None):
    """
    open the simfile in the Lambda application

    pipe: 'pipe' --> will call the subproc. with stdout=subprocess.PIPE
          'log'  --> will log stdout to 'lambda.log'
          None   --> does not pipe stdout
    """
    simfile = os.path.abspath(simfile)
    if " " in simfile:
        simfile = '"%s"' % simfile
    args = ["-file", simfile]
    if vis:
        args.append("-vis")
    if walls:
        args.append("-walls")
    if contrast is not None:
        args.extend(['-contrast', contrast])
    if cmap is not None:
        args.extend(['-colormap', cmap])
    if fps is not None:
        args.extend(['-fps', fps])
    if skip is not None:
        args.extend(['-skip', skip])
    return call_lambda(args, stdout=pipe)


def color_distance_rgb_abs(color1, color2):
    """Calculate the euclidian distance between these two colors

    :param color1: (r, g, b) or an image array of shape (Y, X, 3)
    :param color2: (r, g, b)
    :return: the euclidian distance
    """
    if isinstance(color1, numpy.ndarray):
        return _img_color_distance(color1, color2)
    r1, g1, b1 = color1
    r2, g2, b2 = color2
    return numpy.sqrt((r1-r2)**2 + (g1-g2)**2 + (b1-b2)**2)


def color_distance_rgb(color1, color2, maxvalue=255):
    """Calculate the normalized distance between these two colors
    (the distance will be between 0-1)

    :param color1: (r, g, b) or an image of shape (Y, X, 3)
    :param color2: (r, g, b)
    :return: the ditance (0-1) between these two colors

    Example
    =======

    Calculate a mask indicating where the image is near to a certain
    color

    im = png_load(pngfile)
    distances = color_distance_rgb(im, (0, 0, 255))
    mask = (distances < 0.1).astype(int)
    selected_color = im*mask
    """
    distance = color_distance_rgb_abs(color1, color2)
    max_distance = color_distance_rgb_abs((maxvalue, maxvalue, maxvalue), (0, 0, 0))
    return distance/max_distance


def _img_color_distance(img, color):
    """calculate the distance of each pixel to a given color

    :param img: a numpy.array of shape (y, x, 3)
    :param color: a color tuple (r, g, b)
    :return: a numpy.array of shape (y, x) where each value is a float 0-1
             representing the distance to this color
    """
    return numpy.sqrt(numpy.sum((img - color)**2, 2))


def _pypng_load(pngfile):
    """
    read png, discard alpha channel
    """
    import png
    print("using backend: pyPNG")
    X, Y, rows_iter, info = png.Reader(filename=pngfile).asDirect()
    # check type
    if info['greyscale'] or info['bitdepth'] != 8:
        raise ValueError("only 24-bit color PNGs are supported")
    rows = [numpy.asarray(row, dtype=numpy.uint8) for row in rows_iter]
    alpha = info['alpha']
    rows2 = []
    for row in rows:
        if not alpha:
            row.shape = (X, 3)
        else:
            row.shape = (X, 4)
            row = row[:,0:3]
        rows2.append(row)
    mat = numpy.vstack(rows2)
    print(mat.shape)
    mat.shape = (Y, X, 3)
    return mat


def _pil_load(pngfile):
    from PIL import Image
    def fromimage(image, flatten=0):
        if not Image.isImageType(image):
            raise TypeError("NOT a PIL Image")
        if flatten:
            image = image.convert('F')
        elif image.mode == '1':
            image.convert('L')
        return numpy.array(image)
    im = Image.open(pngfile)
    Y, X, planes = im.shape
    assert planes == 3
    return fromimage(im)


def png_load(pngfile):
    """
    Load a PNG file. Returns a numpy matrix with shape (y, x, 3)

    png = png_load(pngfile)
    pix = png[4, 5]
    r, g, b = pix

    :param pngfile: the path of the png file
    :return: a numpy array of shape (y, x, 3)
    """
    #backends = [_pil_load, _pypng_load]
    backends = [_pypng_load, _pil_load]
    for backend in backends:
        try:
            img = backend(pngfile)
            return img
        except ImportError:
            pass
    raise ImportError("needs either scipy or pypng to load a png")


def png_create(x, y, path, color=(0, 0, 0)):
    """
    create a PNG of given size

    :param x: x size (pixels)
    :param y: y size (pixels)
    :param path: path of generated png
    :param color: color (r,g,b) to fill the generated png
    """
    if color is None:
        color = (0, 0, 0)
    mat = numpy.ones((y, x))
    r, g, b = color
    return png_save(mat, path, lambda x:(x*r, x*g, x*b))


def _get_colormap(colormap=None):
    wall_r, wall_g, wall_b = config['pngcolors']['wall']
    colormaps = {
        'greyscale': lambda v: (v*255, v*255, v*255),
        'wall' :     lambda v: (v*wall_r, v*wall_g, v*wall_b)
    }
    if colormap is None:
        out = _get_colormap('greyscale')
    elif callable(colormap):
        out = colormap
    elif isinstance(colormap, tuple):
        r, g, b = colormap
        out = lambda v: (v*r, v*g, v*b)
    elif colormap in colormaps:
        out = colormaps.get(colormap)
    else:
        raise ValueError("Colormap not understood"
                         "Expecting a color (r, g, b) tuple, a callable f(x) -> (r, g, b)"
                         "or a preset: 'greyscale', 'wall', etc.")
    return out


def png_save(mat, path, colormap=None):
    """Save a monochrome image matrix as a png

    :param mat: a numpy array with shape=(height, width), with float values from 0-1
    :param path: the path to save the matrix to
    :param colormap: a function(x) -> (r, g, b), or the name of one
                     a preset ('grey', 'wall')
                     a color (r, g, b)
    :return: None
    """
    if len(mat.shape) == 3:
        return _png_save_color(mat, path)
    colormap = _get_colormap(colormap)
    mat3 = apply_colormap(mat, colormap)
    return _png_save_color(mat3, path)


def _png_save_color(mat, path):
    import png
    Y, X, b = mat.shape
    bitdepth = 8 # TODO: detect bitdepth
    w = png.Writer(width=X, height=Y, greyscale=False, bitdepth=bitdepth)
    with open(path, "wb") as f:
        w.write_array(f, mat.flatten())


def img_to_mask(path, color, dist=0.1):
    """
    Convert a color image to a binary mask

    :param path: The path of the image
    :param color: The color to use as mask
    :param dist: The max. distance to the color, as a value between 0-1
    :return:
    """
    ext = os.path.splitext(path)[1].lower()
    if ext == '.png':
        img = png_load(path)
        mask = (color_distance_rgb(img, color) <= dist).astype(int)
    else:
        raise ValueError("Image format not supported")
    return mask


def mask_to_pixels(mask):
    """
    generate an iterator over the non-zero items of the mask

    Example: get the realworld coords. of blue pixels
             defined in an image

    mask = img_to_mask("myimage.png", (0,255,255))
    for xpix, ypix in mask_to_pixels(mask):
        x, y = sim.pix2coords(x, y)
        print(x, y)

    """
    Y, X = mask.nonzero()
    return zip(X, Y)


def png_info(path):
    """Returns a dict with info about the png"""
    import png
    r = png.Reader(filename=path)
    x, y, frames, info = r.read()
    return info

def apply_colormap(mat, colormap=None):
    """
    Apply a function returning (r, g, b) to a monochrome matrix defined between 0, 1

    :param mat: a matrix of shape (y, x) of values between 0-1
    :param colormap: a callable f(x) -> (r, g, b)
    :return: a matrix of shape (y, x, 3) of values between 0-255
    """
    colormap = _get_colormap(colormap)
    Y, X = mat.shape
    r, g, b = colormap(mat)
    mat3 = numpy.zeros((Y, X*3))
    mat3[:,0::3] = r
    mat3[:,1::3] = g
    mat3[:,2::3] = b
    mat3.shape = (Y, X, 3)
    return mat3


def interpolate_color(x, color1, color2):
    """

    :param x: (0-1)
    :param color1: (r, g, b)
    :param color2: (r, g, b
    :return:
    """
    r0, g0, b0 = color1
    r1, g1, b1 = color2
    r = x*(r1-r0) + r0
    g = x*(g1-g0) + g0
    b = x*(b1-b0) + b0
    return (r, g, b)

