
from iscard import core, Model
import pytest 
import common as cm 
import numpy as np 
import tempfile
import pandas as pd
#Create a common fake model 

model = cm.create_fake_model(noise_mean=0, noise_std=1)




def test_create_inter_samples_model():

    model.create_inter_samples_model()
    assert model.inter_model is not None    
    assert len(model.raw) == len(cm.get_bed_interval(cm.BEDFILE))
    assert len(model.inter_model) == len(cm.get_bed_interval(cm.BEDFILE))

    assert "chrom" in model.inter_model.columns
    assert "pos" in model.inter_model.columns
    assert "name" in model.inter_model.columns
    assert "mean" in model.inter_model.columns
    assert "std" in model.inter_model.columns
    assert "min" in model.inter_model.columns
    assert "max" in model.inter_model.columns


def test_create_intra_samples_model():

    model.create_intra_samples_model()
    assert model.intra_model is not None

    bed_interval = cm.get_bed_interval(cm.BEDFILE)
    bed_interval = bed_interval[bed_interval.index % model.step == 0]
    assert len(model.intra_model) == len(bed_interval)

    assert "chrom" in model.intra_model.columns
    assert "pos" in model.intra_model.columns
    assert "name" in model.intra_model.columns
    assert "idx" in model.intra_model.columns
    assert "corr" in model.intra_model.columns
    assert "coef" in model.intra_model.columns
    assert "intercept" in model.intra_model.columns
    assert "std" in model.intra_model.columns


def test_to_hdf5():
    
    model_file = tempfile.gettempdir() + "/test_iscard_model.h5"
    
    a_model = Model()
    bamlist = cm.create_fake_bamlist(2)
    a_model.learn(bamlist, cm.BEDFILE, show_progress = False)

    a_model.to_hdf5(model_file)

    hdf = pd.HDFStore(model_file)

    hdf_keys = [key[1:] for key in hdf.keys()]

    assert "bamlist" in hdf_keys
    assert "raw" in hdf_keys
    assert "metadata" in hdf_keys
    assert "inter_model" in hdf_keys
    assert "intra_model" in hdf_keys

    hdf.close()

    # Todo test read ... 

def test_from_hdf5():
    model_file = tempfile.gettempdir() + "/test_iscard_model.h5"
    a_model = Model()
    bamlist = cm.create_fake_bamlist(2)
    a_model.learn(bamlist, cm.BEDFILE, show_progress = False)
    a_model.to_hdf5(model_file)

    second_model = Model(model_file)

    assert second_model.bamlist == a_model.bamlist




def test_learn():
    a_model = Model()
    bamlist = cm.create_fake_bamlist(2)
    a_model.learn(bamlist, cm.BEDFILE, show_progress = False)

    assert a_model.intra_model is not None
    assert a_model.inter_model is not None


