lcdblib
=======

.. image:: https://travis-ci.org/lcdb/lcdblib.svg?branch=master
    :target: https://travis-ci.org/lcdb/lcdblib

`lcdblib` is a toolbox of functions and classes for working with genomic data.
While developing code for other projects, if we decide something could be more
generally useful then we add it to this package.

Installation
------------

Installation is intended to be from `conda`, with the `bioconda
<https://bioconda.github.io/>`_ channel and dependencies set up::

    conda install -c lcdb lcdblib


You may want to permanently add the `lcdb` channel to your ``.condarc`` to
avoid specifiying it each time::

    conda config --add channels lcdb

Then::

    conda install lcdblib

or if it's already installed::

    conda update lcdblib

Automated tests are run on `travis-ci
<https://travis-ci.org/lcdb/lcdblib/builds>`_, and successful builds on the
`master` branch are automatically uploaded to the `lcdb conda channel
<https://anaconda.org/lcdb/lcdblib/files>`_.
