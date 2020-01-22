============
Installation
============

Iscard is being developed and tested under Linux with python3.7. 

Quick installation
------------------

The easiest way to install iscard is to use ``pip`` on the command line::

    pip install iscard

This will download the software from `PyPI (the Python packaging
index) <https://pypi.python.org/pypi/iscard/>.
If an old version of iscard exists on your system, the ``--upgrade`` parameter is required in order
to install a newer version. 

.. _dependencies:

Dependencies
------------

Iscard depends of the following on pysam and pysamstats to compute depth coverage.
Those libraries should be installed automatically with iscard. 
A C compiler is possible required. 


Uninstalling
------------

Type  ::

    pip uninstall iscard

and confirm with ``y`` to remove the package. Under some circumstances, multiple
versions may be installed at the same time. Repeat the above command until you
get an error message in order to make sure that all versions are removed.

