from iscard import core

import pytest
import pandas as pd
from tempfile import mkdtemp
import shutil 
import os 
import glob 
import numpy as np
import common as cm
import contextlib

def create_fake_bamlist(count = 10):
	
	tmp = mkdtemp()

	bams = []
	for i in range(count):
		bamfile = tmp + os.path.sep + f"fake_{i}.bam" 

		shutil.copyfile(cm.BAMFILE, bamfile )
		shutil.copyfile(cm.BAMFILE +".bai", bamfile +".bai" )

		bams.append(bamfile)
	return bams


def test_bed_reader():
	
	reader = core.read_bed(cm.BEDFILE)
	assert set(reader.columns) == set(["chrom","start","end","name"])
	assert len(reader) == len(open(cm.BEDFILE).readlines())
	assert reader.isna().any().any() == False


def test_get_coverage():
	start = 1952
	end = 3230
	df = core.get_coverage(cm.BAMFILE, "chr1", 1952,3230)
	assert len(df) == end - start


@pytest.mark.parametrize("sample_rate", [1,10,100,200])
def test_get_coverages_from_bed(sample_rate):

	df = core.get_coverages_from_bed(cm.BAMFILE,cm.BEDFILE, sample_rate = sample_rate, show_progress= False)

	assert list(df.index.names) == ["name","chrom", "pos"]
	assert len(df.columns) == 1

	bed = core.read_bed(cm.BEDFILE)
	
	# Loop over each bed region and check if position fit with sample_rate	
	df = df.droplevel(0)
	for _, line in bed.iterrows():
		chrom, start, end = line["chrom"], line["start"], line["end"]
		pos_index = (df.loc[chrom, :].query(" `pos` >= @start & `pos` <= @end ")).index.to_series()
		assert ((pos_index.diff()[1:] == sample_rate).all()) == True

	
def test_compute_coverage():

	names = list(core.read_bed(cm.BEDFILE)["name"].unique())
	sample_count = 10
	bamlist = cm.create_fake_bamlist(sample_count)
	df = core.compute_coverage(bamlist,cm.BEDFILE, show_progress= False)

	assert list(df.columns) == [f"fake_{i}" for i in range(sample_count)]
	assert list(df.index.names) == ["name","chrom", "pos"]
	

def test_scale_dataframe():

	df = pd.DataFrame({
		"A":[1,2,3,4],
		"B":[10,20,30,40]
 		})

	#If scalled, A and B should be same 
	df_sum = list(core.scale_dataframe(df).sum())
	assert df_sum[0] == df_sum[1]


def test_call_region():

	serie = pd.Series([1,2,4,4,4,4,3,1])

	for region in core.call_region(serie, threshold=3, consecutive_count=3 ):
		begin, end = region 
		assert begin == 2 
		assert end == 5

def test_bedgraph():

	df = pd.DataFrame({"chrom":[1,2], "pos":[10,20], "depth":[33,44]})

	with contextlib.redirect_stdout(open("/dev/null","w")):
		core.print_bedgraph(df, "depth", "name")
