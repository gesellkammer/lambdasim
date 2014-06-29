import warnings

try:
    import matplotlib.pyplot as plt
    MATPLOTLIB = True
except ImportError:
    MATPLOTLIB = False

def matplotlib_available():
    if MATPLOTLIB:
        return True
    warnings.warn("Matplotlib not available. Plotting will not be available")
    return False

def plot_array(A, xsize=None, ysize=None, cmap=None):
    """
    A: a 2D numpy array
    xsize, ysize: the "realworld" size of the array
    """
    Y, X = A.shape
    if xsize is None:
        xsize = X 
    if ysize is None:
        ysize = Y
    if matplotlib_available():
        plt.imshow(A, origin='upper', extent=[0, xsize, 0, ysize], cmap=cmap)