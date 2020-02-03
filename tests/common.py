
import shutil 
import os 
from tempfile import mkdtemp
import csv 
import pandas as pd
import numpy as np
from iscard import core, Model  

BEDFILE = "data/region.bed"
BAMFILE = "data/test.bam"



def create_fake_bamlist(count = 10):
	
	tmp = mkdtemp()

	bams = []
	for i in range(count):
		bamfile = tmp + os.path.sep + f"fake_{i}.bam" 

		shutil.copyfile(BAMFILE, bamfile )
		shutil.copyfile(BAMFILE +".bai", bamfile +".bai" )

		bams.append(bamfile)
	return bams


def get_bed_interval(bedfile):

	with open(bedfile, "r") as file:
		reader = csv.reader(file, delimiter="\t")

		data = []
		for line in reader:
			chrom, start, end, name = line
			start = int(start)
			end = int(end)
			for i in range(start,end):
				data.append([chrom, start, end, name])

		df = pd.DataFrame(data, columns = ["chrom", "start","end", "name"])

		return df

def create_fake_model(size=10, noise_mean = 0, noise_std = 1):
	model = Model()
	bamlist = create_fake_bamlist(size)

	model.bamlist = bamlist
	model.bedfile = BEDFILE
	model.raw = core.get_coverages_from_bed(model.bamlist,model.bedfile, show_progress = True)
	# Add fake noise 
	noise = np.random.normal(noise_mean, noise_std, model.raw.shape)
	model.raw = model.raw.add(noise,axis=0).astype(int)

	model.create_inter_samples_model()
	model.create_intra_samples_model()

	return model
