.. iscard documentation master file, created by
   sphinx-quickstart on Wed Jan 22 16:12:21 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to iscard's documentation!
==================================

Iscard is a tool to aid the detection of `Copy Number Variation <https://en.wikipedia.org/wiki/Copy-number_variation>`_ from Next Generation Sequencing data. It works by computing 2 models from a normal training data set:    

* **inter-model** computes the sample depth deviation from a normal distribution estimated from the training data set. 

* **intra-model** computes the correlation of depths at different position within a sample through linear regression. These correlations are also obtained from the training data set. 
  
After creating your model, you can test a new sample against it and get different metrics such as **z-scores** which describe the deviation from the model. CNVs can then be called using these metrics or displayed in a genome browser such as IGV.

.. figure:: pkd.png

   The figure above display 3 samples with CNVs of different sizes. Each sample comes with 2 tracks describing inter-model and intra-model z-scores.

  

  


.. toctree::
   :maxdepth: 2

   quickstart
   installation
   iscard

   

   ...





Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

