Quickstart
===========================


Iscard create a model by looking on differents normal sample bam file.
After this step, you can test a new sample against this model and get differents metrics like z-score deviation.

Installation 
------------

Install your package using pip or conda as it is explained on `installation`_ ::

  pip install iscard


Creating the model 
-----------------

Suppose you have all normal sample bam file under the *normal* folder.
To create the model, type the following command ::

  iscard learn -b normal/*.bam -r region.bed -o model.h5  

Test a sample
-------------

After the model creation, you can test a sample against the model ::


    iscard test -b sample.bam -m model.h5 -o test.h5

