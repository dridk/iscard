
from iscard import core
from iscard.model import Model, call_test, plot_test
import pytest 
import common as cm 
import numpy as np 
import tempfile
import pandas as pd
import tempfile
#Create a common fake model 



@pytest.fixture
def fake_model(scope='session'):
    model = cm.create_fake_model(size=5, sample_rate = 100, noise_mean=0, noise_std=1)
    return model



def test_fake_model(fake_model):
    assert fake_model is not None  

def  test_get_group_names(fake_model):
    assert set(["GENEA","GENEB"])  == set(fake_model.get_group_names())


def test_create_inter_samples_model(fake_model):

    fake_model.create_inter_samples_model()
    assert fake_model.inter_model is not None    

    assert len(fake_model.inter_model) == len(cm.get_bed_interval(cm.BEDFILE, sample_rate = fake_model.sample_rate ))

    assert "mean" in fake_model.inter_model.columns
    assert "std" in fake_model.inter_model.columns
    assert "min" in fake_model.inter_model.columns
    assert "max" in fake_model.inter_model.columns

def test_create_intra_samples_model(fake_model):

    fake_model.create_intra_samples_model()
    assert fake_model.intra_model is not None

    bed_interval = cm.get_bed_interval(cm.BEDFILE,  sample_rate = fake_model.sample_rate )
    #bed_interval = bed_interval[bed_interval.index % model.sampling == 0]
    assert len(fake_model.intra_model) == len(bed_interval)

    assert "idx" in fake_model.intra_model.columns
    assert "corr" in fake_model.intra_model.columns
    assert "coef" in fake_model.intra_model.columns
    assert "intercept" in fake_model.intra_model.columns
    assert "std2" in fake_model.intra_model.columns

def test_learn():
    bamlist = cm.create_fake_bamlist(2)
    model = Model()
    model.learn(bamlist, cm.BEDFILE, show_progress = False)
    assert model.intra_model is not None
    assert model.inter_model is not None


def test_to_hdf5(fake_model):
    
    model_file = tempfile.gettempdir() + "/test_iscard_model.h5"
    
    bamlist = cm.create_fake_bamlist(2)
    fake_model.learn(bamlist, cm.BEDFILE, show_progress = False)

    fake_model.to_hdf5(model_file)

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
    pass
    # model_file = tempfile.gettempdir() + "/test_iscard_model.h5"
    # model = Model()
    # bamlist = cm.create_fake_bamlist(2)
    # model.learn(bamlist, cm.BEDFILE, show_progress = False)
    # model.to_hdf5(model_file)

    # second_model = Model(model_file)

    # assert second_model.bamlist == model.bamlist


def test_sample(fake_model):

    output = fake_model.test_sample(cm.BAMFILE)
    assert len(output) == len(fake_model)

def test_call_test(fake_model):
    test_data = fake_model.test_sample(cm.BAMFILE, show_progress=False)
    calling = list(call_test(test_data))
    assert calling == []

def test_plot_test(fake_model):
    test_data = fake_model.test_sample(cm.BAMFILE, show_progress=False)
    #plot_test(tempfile.mkstemp("png"), test_data, model, "GENEA")
