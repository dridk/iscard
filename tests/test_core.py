from iscard import core

import pandas as pd
from tempfile import mkdtemp
import shutil 
import os 
import glob 

import common as cm



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


def test_get_coverages_from_bed():
	
	names = list(core.read_bed(cm.BEDFILE)["name"].unique())
	sample_count = 10
	bamlist = cm.create_fake_bamlist(sample_count)
	df = core.get_coverages_from_bed(bamlist,cm.BEDFILE, show_progress= False)

	assert set(df.columns) == set([f"fake_{i}" for i in range(sample_count)])

def test_scale_dataframe():

	df = pd.DataFrame({
		"A":[1,2,3,4],
		"B":[10,20,30,40]
 		})

	#If scalled, A and B should be same 
	df_sum = list(core.scale_dataframe(df).sum())
	assert df_sum[0] == df_sum[1]
