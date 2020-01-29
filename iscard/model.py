import pandas as pd 
from sklearn.metrics import pairwise_distances_chunked, pairwise_distances
from itertools import product
import numpy as np
import os
from iscard import core


class Model(object):
    def __init__(self, modelfile = None):
        super().__init__()

        self.raw = None
        self.bamlist = []

        self.inter_model = None
        self.intra_model = None
        self.step = 100

        if modelfile:
            self.from_hdf5(modelfile)


    def learn(self, bamlist: list, bedfile: str, show_progress = True):
        """Create intrasample and intersample model
        
        Args:
            bamlist (list): List of bam files 
        """
        self.bamlist = bamlist
        self.bedfile = bedfile

        self.raw = core.get_coverages_from_bed(self.bamlist, self.bedfile, show_progress = show_progress)

        self.create_inter_samples_model()
        self.create_intra_samples_model()

    def create_inter_samples_model(self):

        self.norm_raw = core.scale_dataframe(self.raw)
        self.inter_model = pd.DataFrame(
            {
                "mean": self.norm_raw.mean(axis=1),
                "median": self.norm_raw.median(axis=1),
                "std": self.norm_raw.std(axis=1),
                "min": self.norm_raw.min(axis=1),
                "max": self.norm_raw.max(axis=1),
            }
        ).reset_index()

    def create_intra_samples_model(self):
        
        # Keep row every step line 
        # reset index because we are going to work on integer index
        sub_raw = self.raw.reset_index()
        sub_raw = sub_raw[sub_raw.index % self.step == 0]

        # Create Mask index 
        # This is used to avoid pairwise comparaison within same name   
        # For example, if name is = [A,A,A,B,B,C], it computes the following mask 

        #   A A A B B C
        # A 0 0 0 1 1 0
        # A 0 0 0 1 1 0 
        # A 0 0 0 1 1 0
        # B 1 1 1 0 0 1
        # B 1 1 1 0 0 1
        # C 1 1 1 1 1 0

        index = sub_raw.name
        mask = np.array([i[0] == i[1] for i in product(index,index)]).reshape(len(index),len(index))

        # return to multiindex 
        sub_raw = sub_raw.set_index(["name","chrom","pos"])
        
        def _reduce(chunk, start):
            """This function is called internally by pairwise_distances_chunked
            @see https://scikit-learn.org/stable/modules/generated/sklearn.metrics.pairwise_distances_chunked.html

            This function looks for the maximum correlation  value in the chunk matrix and return the id 
            Same name in pairwise are skiped by the mask 
            For example :  
                  A   B   C 
            A    NA  0.9  0.8
            B    0.5 NA  0.4
            C    0.3 0.7  NA
            
            Will return a dataframe :  
            
             id   idx  corr 
             A     B    0.9
             B     C    0.4
             C     B    0.9

            """
            # skip na value
            chunk[np.isnan(chunk)] = 1
            # correlation metrics from sklearn is 1 - corr 
            chunk = 1 - chunk 
            rows_size = chunk.shape[0]
            
            select_mask = mask[start:start+rows_size]
            # looks for id of maximum correlation value 
            idx  = np.argmax(np.ma.masked_array(chunk, select_mask),axis=1)    
            
            # We only get idx, let's get correlation value 
            corr = []
            for i, index in enumerate(idx):
                corr.append(chunk[i][index])
            
            # Create a dataframe 
            return pd.DataFrame({"idx": idx, "corr": corr}, index = range(start, start + rows_size))

        # Perform pairwise correlation by using pairwise_distances_chunked to avoid memory limit 
        
        all_reduce_chunk = []
        
        for chunk in pairwise_distances_chunked(sub_raw, metric="correlation",reduce_func = _reduce, n_jobs=10,working_memory=1):
            all_reduce_chunk.append(chunk)

        self.intra_model  = pd.concat(all_reduce_chunk)
        self.intra_model.index = sub_raw.index
        self.intra_model.reset_index(inplace=True)

        ## Let's now add linear regression parameter for each line in the model 
        for index, row in self.intra_model.iterrows():
            # create source multiindex and target multiindex  
            name,chrom,pos = tuple(row[["name","chrom","pos"]])
            target_name,target_chrom, target_pos = tuple(self.intra_model.iloc[row["idx"],:][["name","chrom","pos"]])
            
            if row["corr"] == 0:
                continue

            # Get values and compute linear regression 
            x = self.raw.loc[(name,chrom,pos)]
            y = self.raw.loc[(target_name,target_chrom,target_pos)]
            
            try:
                coef, intercept = tuple(np.polyfit(x,y,1))
            except:
                # Todo use debug
                print("error for ", name, chrom,pos)
                coef, intercept = 0, 0
            
            # Compute residual and standard deviation 
            yp = coef * x + intercept
            residual = yp - y
            std = residual.std()
            
            self.intra_model.loc[index,"coef"] = coef
            self.intra_model.loc[index,"intercept"] = intercept    
            self.intra_model.loc[index,"std"] = std    


    def to_hdf5(self, filename):
        self.raw.to_hdf(filename,"raw")

        self.inter_model.to_hdf(filename,"inter_model")
        self.intra_model.to_hdf(filename,"intra_model")
        
        pd.Series(self.bamlist).to_hdf(filename,"bamlist")
        pd.Series(
            {
            "step": self.step,
            "region":os.path.abspath(self.bedfile)

            }).to_hdf(filename, key="metadata")

        

    def from_hdf5(self, filename):
        self.raw = pd.read_hdf(filename,"raw")

        self.inter_model = pd.read_hdf(filename,"inter_model")
        self.intra_model = pd.read_hdf(filename,"intra_model")
        
        self.bamlist = list(pd.read_hdf(filename,"bamlist"))

        metadata = pd.read_hdf(filename, key="metadata")
        self.step = metadata["step"]

