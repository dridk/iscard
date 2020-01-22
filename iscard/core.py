import numpy as np 
import pandas as pd 
import pysamstats
import pysam
from sklearn.preprocessing import MinMaxScaler, StandardScaler,MaxAbsScaler
import matplotlib.pyplot as plt

import os 
import re 
from tqdm import tqdm
import pprint 

def read_bed(filename:str) -> pd.DataFrame:
    """return a dataFrame from bed bedfile 

    Bedfile must have 4 colonnes : chr, start, end, name.
    Name can be the gene name or any other group name. 
    
    Args:
        filename (str): bed file
    
    Returns:
        pd.DataFrame: a Dataframe with 4 colonnes : chr, start,end, name
    """
    return pd.read_csv(filename, sep="\t", header=None, names=["chrom","start","end","name"])



def get_coverage(bamfile: list, chrom: str, start: int, end: int, window = 1, agg = np.mean) -> pd.DataFrame:
    '''Get read depth from from bam files according location  
    
    Args:
        bamfile (list): list of bam file with index
        chrom (str): chromosome name
        start (int): start position of interest
        end (int): end position of interest
        window (int, optional): window size for grouping. By default no group are performed
        agg (function, optional): Aggregate function for grouping. Can be a litteral like (mean,max,min,sum, mediam)
    
    Returns:
        pd.DataDrame: A dataframe with chromosom, position and depth for each bam file 
    
    '''
    mybam = pysam.AlignmentFile(bamfile)
    
    df = pd.DataFrame(pysamstats.load_coverage(mybam,chrom=chrom,start=start,end=end, truncate=True, pad=True, fields = ["chrom", "pos", "reads_all"]))
    
    df["chrom"] = df["chrom"].apply(lambda x : x.decode())
    df.rename({"reads_all": "depth"}, inplace=True, axis=1)
    
    if window and agg:
        if window > 1:
        # Group every window line  and aggregate depth 
            return df.groupby(df.index // window).agg({"chrom":"first", "pos":"first" , "depth":agg})

    return df 


def get_coverages_from_bed(bamlist, bedfile, window=1, agg="mean"):
    '''Get coverage dataframe from one of many bam file within a bedfile
    
    Args:
        bamlist (str): list of bam files
        bedfile (str): a bed file with 4 colonnes (chr,start,end,name)
        window (int, optional): window size 
        agg (str, optional): Aggregate function 
    
    Returns:
        pd.DataFrame: a Dataframe with depth for each bam file within region describe in bed 
    '''
    
    bed = read_bed(bedfile)
    
    all_df = []
    
    tqdm_range =  tqdm(bed.iterrows(), total = len(bed), desc="Read bam files", leave=True)

    for _, row in tqdm_range:

        chrom = row["chrom"]
        start = row["start"]
        end   = row["end"]
        name  = row["name"]

        tqdm_range.set_description(f"analysing region {chrom}:{start}-{end}")
        tqdm_range.refresh() 

        
        for index, bam in  enumerate(bamlist):

            sample_name = os.path.basename(bam).replace(".bam","")
      
            sample_df = get_coverage(bam, chrom, start, end, window, agg= lambda x: np.around(np.mean(x), 2))

            if index == 0:
                df = sample_df
                df.insert(2,"name",name)
                df = df.rename({"depth":sample_name}, axis=1)
            
            else:
                df[sample_name] =sample_df["depth"]
                
        all_df.append(df)
    
    final = pd.concat(all_df).reset_index(drop=True)   
    final = final.set_index(["name","chrom","pos"])
    return final
            

def scale_dataframe(df):
    new_df = df.copy()
    scaler = MaxAbsScaler()
    # scale data
    scale_df = pd.DataFrame(scaler.fit_transform(new_df), columns=new_df.columns)

    # Need to apply same index for scale_df before updating values
    new_df.update(scale_df.set_index(new_df.index))
    
    return new_df



def create_model(bamlist: list, bedfile:str, output:str, window = 1, agg="mean"):
    """ create hdf5 model """ 

    #write bam list 
    pd.Series([os.path.abspath(i) for i in bamlist]).to_hdf(output, key="bamlist")

    # write metadata
    pd.Series({"window": window, "agg":agg, "region":os.path.abspath(bedfile)}).to_hdf(output, key="metadata")

    print("compute model")
    # compute and write raw dataframe 
    raw = get_coverages_from_bed(bamlist,bedfile, window = 1, agg = agg )
    raw.to_hdf(output, "raw")

    # Scale 
    scale_raw = scale_dataframe(raw)
    model = pd.DataFrame(
    {
        "mean":scale_raw.mean(axis=1),
        "median":scale_raw.median(axis=1),
        "std": scale_raw.std(axis=1),
        "min": scale_raw.min(axis=1),
        "max": scale_raw.max(axis=1),
    })

    model.to_hdf(output,key="model")


def test_sample(bam: str, model: str, output: str):

    metadata = pd.read_hdf(model, "metadata")
    model = pd.read_hdf(model, "model")

    region = metadata["region"]
    window = metadata["window"]
    agg = metadata["agg"]

    sample_df = get_coverages_from_bed([bam],region, window = window, agg = agg )
    sample_df = scale_dataframe(sample_df)

    model = model.reset_index()

    model["sample"] = sample_df.reset_index(drop=True).iloc[:,0]
    model["sample_z"]  = ( model["sample"] - model["mean"]) / model["std"]
    model["sample_zsmooth"] = model["sample_z"].rolling(300).mean()

    model.to_hdf(output, key="sample")

def bedgraph_sample(testfile:str, column = "sample_z"):

    sample_df = pd.read_hdf(testfile, "sample")

    if column not in sample_df.columns:
        print("select another column")
        return 

    #sample_df[column] = sample_df[column].apply(round)

    with open("/dev/stdout","w") as file:
        header = f"track type=bedGraph name=\"Iscard sample\" description=\"Iscard sample\" visibility=full color=200,100,0 altColor=0,100,200 priority=20"
        file.write(header + "\n")
        sample_df[["chrom","pos","pos",column]].to_csv(file , sep="\t", header=False, index=False)


def detect_outside_region(sample_df, threshold = 2, times = 100):
    counter = 0
    valid = False
    for index, i in sample_df.items():
        if abs(i) > threshold:
            counter+= 1
            if counter == 1:
                begin = index[1]
        else:
            counter = 0 
            if valid == True:
                valid =  False
                yield (begin,end)
            
        if counter > times:
            end = index[1]
            valid = True


def plot_sample(testfile: str, name:str, coordinate:str, output:str):

    if not name:
        print("select a gene name ")
        return 

    match = re.search(r'(chr\w+)\:(\d+)-(\d+)', coordinate )

    if  match:
        (chrom,start, end) = match.groups()
        print(chrom,start, end)

    sample_df = pd.read_hdf(testfile, "sample")

    df = sample_df.query("name == @name")

    if match:
        df = query("pos >@start & pos < @end ")



    figure, ax = plt.subplots(2,1, figsize=(30,10))
    plt.subplots_adjust(hspace = 0.3)
    
    figure.suptitle("PKD1", fontsize=30)
    
    ax[0].grid(True)
    ax[0].set_xlabel('position')
    ax[0].set_ylabel('raw depth')

    ax[0].plot("pos", "sample", data = df, color="red")
    ax[0].plot("pos", "mean", data = df, color="#455c7c", ls='--', lw=1)
    ax[0].fill_between("pos","min", "max", color="blue", alpha=0.2, data = df)
    handles, labels = ax[0].get_legend_handles_labels()
    ax[0].legend(handles, labels)

    ax[1].grid(True)
    ax[1].set_xlabel('position')
    ax[1].set_ylabel('z-score')

    ax[1].plot("pos", "sample_z", data = df, color="#ff6961")
    ax[1].plot("pos", "sample_zsmooth", data = df, color="green", linewidth=2)

    ax[1].set_ylim(-10,10)
    handles, labels = ax[1].get_legend_handles_labels()
    ax[1].legend(handles, labels)

    d = sample_df.set_index(["chrom","pos"])


    # detect region 
    outside_regions = detect_outside_region(df.set_index(["chrom","pos"])["sample_z"], threshold=2, times=500)

    for region in outside_regions:
        ax[0].axvspan(*region, alpha=0.5, color='lightgray')
        ax[1].axvspan(*region, alpha=0.5, color='lightgray')



    figure.savefig(output)





def print_model_info(model_file:str):

    pp = pprint.PrettyPrinter(indent=4)

    bamlist = list(pd.read_hdf(model_file, key="bamlist"))

    model = pd.read_hdf(model_file,key="model")

    print("Regions scanned:\n")
    print("{} position(s)".format(model.shape[0]))

    print("\n")

    print("{} Bam(s) used for the model: \n".format(len(bamlist)))
    print("\n".join(list(pd.read_hdf(model_file, key="bamlist"))))
       
    print("\n")

    print("Model arguments: \n")
    for key, value in dict(pd.read_hdf(model_file, key="metadata")).items():
        print("{:<10}{:<10}".format(key+":",str(value)))



    # # write metadata
    # pd.Series({"window": window, "agg":agg, "region":os.path.abspath(bedfile)}).to_hdf(output, key="metadata")

    # print("compute model")
    # # compute and write raw dataframe 
    # raw = get_coverages_from_bed(bamlist,bedfile, window = 1, agg = agg )
    # raw.to_hdf(output, "raw")

    # # Scale 
    # scale_raw = scale_dataframe(raw)
    # model = pd.DataFrame(
    # {
    #     "mean":scale_raw.mean(axis=1),
    #     "median":scale_raw.median(axis=1),
    #     "std": scale_raw.std(axis=1),
    #     "min": scale_raw.min(axis=1),
    #     "max": scale_raw.max(axis=1),
    # })

    # model.to_hdf(output,key="model")