
from iscard import core, Model
from iscard.model import call_test, plot_test
import pytest 
import common as cm 
import numpy as np 
import tempfile
import pandas as pd
import tempfile
#Create a common fake model 

model = cm.create_fake_model(size=5,noise_mean=0, noise_std=1)


def  test_get_group_names():
    assert set(["GENEA","GENEB"])  == set(model.get_group_names())


def test_create_inter_samples_model():

    model.create_inter_samples_model()
    assert model.inter_model is not None    
    assert len(model.raw) == len(cm.get_bed_interval(cm.BEDFILE))
    assert len(model.inter_model) == len(cm.get_bed_interval(cm.BEDFILE))

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

    assert "idx" in model.intra_model.columns
    assert "corr" in model.intra_model.columns
    assert "coef" in model.intra_model.columns
    assert "intercept" in model.intra_model.columns
    assert "std" in model.intra_model.columns


def test_to_hdf5():
    
    model_file = tempfile.gettempdir() + "/test_iscard_model.h5"
    
    bamlist = cm.create_fake_bamlist(2)
    model.learn(bamlist, cm.BEDFILE, show_progress = False)

    model.to_hdf5(model_file)

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
    model = Model()
    bamlist = cm.create_fake_bamlist(2)
    model.learn(bamlist, cm.BEDFILE, show_progress = False)
    model.to_hdf5(model_file)

    second_model = Model(model_file)

    assert second_model.bamlist == model.bamlist


def test_learn():
    bamlist = cm.create_fake_bamlist(2)
    model.learn(bamlist, cm.BEDFILE, show_progress = False)
    assert model.intra_model is not None
    assert model.inter_model is not None



def test_sample():
    output = model.test_sample(cm.BAMFILE)
    assert len(output) == len(model)

def test_call_test():
    test_data = model.test_sample(cm.BAMFILE)
    calling = list(call_test(test_data))
    assert calling == []

def test_plot_test():
    test_data = model.test_sample(cm.BAMFILE)
    #plot_test(tempfile.mkstemp("png"), test_data, model, "GENEA")
