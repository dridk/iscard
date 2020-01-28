from iscard import core

from tempfile import mkdtemp
import shutil 
import os 
import glob 

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


def test_bed_reader():
	
	reader = core.read_bed(BEDFILE)
	assert set(reader.columns) == set(["chrom","start","end","name"])
	assert len(reader) == len(open(BEDFILE).readlines())
	assert reader.isna().any().any() == False


def test_get_coverage():
	start = 1952
	end = 3230
	df = core.get_coverage(BAMFILE, "chr1", 1952,3230)
	assert len(df) == end - start


def test_get_coverages_from_bed():
	
	names = list(core.read_bed(BEDFILE)["name"].unique())


	sample_count = 10
	bamlist = create_fake_bamlist(sample_count)
	df = core.get_coverages_from_bed(bamlist,BEDFILE)


	assert set(df.columns) == set([f"fake_{i}" for i in range(sample_count)])
