.. image:: pics/lambdaicon128.png

=========
lambdasim
=========

Python library to interact with Lambda_ (FDTD acoustic simulations)

.. image:: pics/muffle.gif

Features
--------

* Read and write simulations
* Create complex simulation setups using vector shapes
* Embed sampled sources in a simulation
* Convert sensor data to audio
* Calculate real-world dimensions of simulation
* Preset shapes like tubes, bells, resonators
* Precalculation of dead nodes to speed-up calculations

Installation
------------

1. Install Lambda. Download binaries_ or install from source_
2. Clone this repo::

    $ git clone https://github.com/gesellkammer/lambdasim


3. Install it::


    $ pip install -r requirements.txt
    $ python setup.py install


Dependencies
------------

* shapely_
* shapelib_
* bpf4_
* sndfileio_
* rasterio_

Examples
--------

The simplest simulation:

.. code:: python

    import lambdasim
    # a simulation space of 4mx3m with resolution of 1 cm
    sim = lambdasim.Simulation((4, 3), resolution=0.01)
    # attach a source
    sim.source_point(2, 1.5, 'sin', freq=440, amp=0.3)
    sim.write("simplesine")
    sim.opensim()

See xxx for more examples

.. _Lambda: https://github.com/gesellkammer/lambda
.. _binaries: https://github.com/gesellkammer/lambda/tree/master/dist
.. _source: https://github.com/gesellkammer/lambda
.. _shapely: https://github.com/sgillies/shapely
.. _shapelib: https://github.com/gesellkammer/shapelib
.. _bpf4: https://github.com/gesellkammer/bpf4
.. _sndfileio: https://github.com/gesellkammer/sndfileio
.. _rasterio: https://github.com/mapbox/rasterio
