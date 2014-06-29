import struct
import array
import numpy

###################################
#  Binary IO
###################################


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



