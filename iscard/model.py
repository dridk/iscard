import pandas as pd 
from sklearn.metrics import pairwise_distances_chunked, pairwise_distances
from itertools import product
import numpy as np
import os
import matplotlib.pyplot as plt

from iscard import core, __version__


class Model(object):
    def __init__(self, modelfile = None):
        super().__init__()

        self.raw = None
        self.bamlist = []
        self.bedfile = None
        self.inter_model = None
        self.intra_model = None
        self.step = 100
        self.test_data = None

        if modelfile:
            self.from_hdf5(modelfile)

    def get_group_names(self):
        """Get group name from self.bedfile
        
        Returns:
            list: unique names list from self.bedfile
        """
        return list(core.read_bed(self.bedfile)["name"].unique())

    def learn(self, bamlist: list, bedfile: str, show_progress = True):
        """Create intrasample and intersample model
        
        Args:
            bamlist (list): List of bam files 
        """

        print("start learning")
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
        )

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

        self.intra_model = pd.concat(all_reduce_chunk)
        ss = sub_raw.reset_index(drop=True)


        for i, row in self.intra_model.iterrows():
            
            j = row["idx"]

            x  = ss.loc[i,:]
            y  = ss.loc[j,:]
            
            try:
                coef, intercept = tuple(np.polyfit(x,y,1))
                yp = x*coef + intercept
                error = yp - y 
                std = error.std()

            
            except:
                coef, intercept = 0,0 
                std = pd.np.NaN

         

         
            self.intra_model.loc[i,"coef"] = coef
            self.intra_model.loc[i,"intercept"] = intercept    
            self.intra_model.loc[i,"std"] = std    


        #self.intra_model.set_index(sub_raw.index, inplace=True)
        self.intra_model = self.intra_model.set_index(sub_raw.index)
        


    def to_hdf5(self, filename):
        self.raw.to_hdf(filename,"raw")

        self.inter_model.to_hdf(filename,"inter_model")
        self.intra_model.to_hdf(filename,"intra_model")
        
        pd.Series(self.bamlist).to_hdf(filename,"bamlist")
        pd.Series(
            {
            "step": self.step,
            "region":os.path.abspath(self.bedfile),
            "version": __version__

            }).to_hdf(filename, key="metadata")

        

    def from_hdf5(self, filename):
        self.raw = pd.read_hdf(filename,"raw")

        self.inter_model = pd.read_hdf(filename,"inter_model")
        self.intra_model = pd.read_hdf(filename,"intra_model")
        
        self.bamlist = list(pd.read_hdf(filename,"bamlist"))

        metadata = pd.read_hdf(filename, key="metadata")
        self.step = metadata["step"]
        self.bedfile =  metadata["region"]


    def test_sample(self,bamfile:str):
    
        del_coverage = core.get_coverages_from_bed([bamfile], self.bedfile)

        dd = del_coverage.copy()
        dd.columns = ["depth"]

        # Compute inter model 
        dd["depth_norm"] = core.scale_dataframe(dd)["depth"]
        dd["inter_z"] = (dd["depth_norm"] - self.inter_model["mean"])/ self.inter_model["std"]


        # Compute intra model 
        subset = dd.loc[self.intra_model.index]
        subset["depth_pair"] = subset.iloc[self.intra_model["idx"],:].iloc[:,0].values
        subset["depth_pair_predicted"] = (self.intra_model["coef"] * subset["depth"]) + self.intra_model["intercept"]
        subset["corr"] = self.intra_model["corr"]
        subset["error"] = subset["depth_pair_predicted"]  - subset["depth_pair"]
        subset["intra_z"] = subset["error"] / self.intra_model["std"]

        subset.drop(["depth","depth_norm","inter_z"], axis=1)
        test_data = dd.join(subset.drop(["depth","depth_norm", "inter_z"], axis=1))
        
        return test_data



    # def plot_test(self, filename:str, group_name:str, call = True,threshold = 2, consecutive_count = 1000 ):
        
    #     data = self.test_data.loc[group_name,:].reset_index()
    #     mm = self.inter_model[["min","max","mean"]].loc[group_name,:].reset_index()

    #     figure, ax = plt.subplots(3,1, figsize=(30,10))
        
    #     figure.suptitle("PKD1", fontsize=30)


    #     ax[0].grid(True)
    #     ax[0].set_xlabel('position')
    #     ax[0].set_ylabel('raw depth')
    #     ax[0].fill_between("pos","min", "max", color="lightgray", alpha=1, data = mm)
    #     ax[0].plot(data["pos"], data["depth_norm"], color="#32afa9")


    #     ax[1].grid(True)
    #     ax[1].set_xlabel('position')
    #     ax[1].set_ylabel('Inter z-score')
    #     ax[1].plot(data["pos"], data["inter_z"], color="lightgray")
    #     ax[1].plot(data["pos"], data["inter_z"].rolling(500).mean(), color="#32afa9")


    #     #ax[0].plot(data["pos"], data["inter_z"].rolling(500).mean())

    #     ax[2].grid(True)
    #     ax[2].set_ylabel('Intra z-score')
    #     ax[2].scatter(data["pos"],data["intra_z"], color="#32afa9")


    #     # Plot region 
    #     if call:
    #         for region in self.call_test(group_name,"inter_z" ,threshold, consecutive_count):
    #             x1, x2 = region 
    #             ax[0].axvspan(x1,x2, alpha=0.5, color='red')
    #             ax[1].axvspan(x1,x2, alpha=0.5, color='red')
    #             ax[2].axvspan(x1,x2, alpha=0.5, color='red')




    # def call_test(self, group_name, column = "inter_z", threshold = 1.96, consecutive_count = 1000):
    #     """Call region from self.test_data
        
    #     Args:
    #         group_name (str)
        
    #     Yields:
    #         TYPE: region 
    #     """
    #     data = self.test_data.loc[group_name,:].reset_index()
    #     for region in core.call_region(data[column], threshold, consecutive_count):

    #         first, last = region 
    #         x1 = data["pos"][first]
    #         x2 = data["pos"][last]
            
    #         yield (x1,x2)


    def __len__(self):
        return len(self.raw)

    def print_infos(self):

        print("Depth position counts: {}".format(len(self.inter_model)))
        print("bedfile: {}".format(self.bedfile))
        print("sample rate: {}".format(self.step))

        print("test data: {}".format(self.test_data is not  None))


        print("Bam(s) used: {}".format(len(self.bamlist)))
        for bam in self.bamlist:
            print("\t - "+bam)

        gps = self.get_group_names()

        print("Group names(s): {}".format(len(gps)))
        for g in gps:
            print("\t - "+g)

        print("Inter model shape: {}".format(self.inter_model.shape))
        print("Intra model shape: {}".format(self.intra_model.shape))





def call_test(test_data: pd.DataFrame, column = "inter_z", threshold = 1.96, consecutive_count = 1000):
    """Call region from self.test_data
    
    Args:
        group_name (str)
    
    Yields:
        TYPE: region 
    """

    for region in core.call_region(test_data[column], threshold, consecutive_count):
        first, last = region 
        chrom = test_data["chrom"][first]
        name = test_data["name"][first]

        x1 = test_data["pos"][first]
        x2 = test_data["pos"][last]
        
        yield (chrom,x1,x2, name)



def plot_test(outputfile:str, test_data:pd.DataFrame, model: Model, group_name:str, call = True,threshold = 2, consecutive_count = 1000 ):
    
    
    data = test_data.query("name == @group_name")
    mm = model.inter_model[["min","max","mean"]].loc[group_name,:].reset_index()

    figure, ax = plt.subplots(3,1, figsize=(30,10))
    
    figure.suptitle(group_name, fontsize=30)


    ax[0].grid(True)
    ax[0].set_xlabel('position')
    ax[0].set_ylabel('raw depth')
    ax[0].fill_between("pos","min", "max", color="lightgray", alpha=1, data = mm)
    ax[0].plot(data["pos"], data["depth_norm"], color="#32afa9")


    ax[1].grid(True)
    ax[1].set_xlabel('position')
    ax[1].set_ylabel('Inter z-score')
    ax[1].plot(data["pos"], data["inter_z"], color="lightgray")
    ax[1].plot(data["pos"], data["inter_z"].rolling(500).mean(), color="#32afa9")


    #ax[0].plot(data["pos"], data["inter_z"].rolling(500).mean())

    ax[2].grid(True)
    ax[2].set_ylabel('Intra z-score')
    ax[2].scatter(data["pos"],data["intra_z"], color="#32afa9")


    # Plot region 
    if call:
        for region in call_test(data):
            chrom, x1, x2, name = region 
            ax[0].axvspan(x1,x2, alpha=0.5, color='red')
            ax[1].axvspan(x1,x2, alpha=0.5, color='red')
            ax[2].axvspan(x1,x2, alpha=0.5, color='red')

    plt.savefig(outputfile)



