from iscard.core import * 
from pybedtools import BedTool
import pickle 
import numpy as np

class Model:
    def __init__(self):
        pass

    def fit(self, bedfile, bamfiles):
        self.bedfile   = bedfile
        self.bamfiles  = bamfiles
        self.matrix    = compute_coverage_matrix(self.bedfile, self.bamfiles)
        self.models    = compute_regressions_models(self.matrix, alpha = 0.1)
        self.amplicons = list(self.matrix.columns)

    def save(self, filename):
        with open(filename, "wb") as file:
            pickle.dump(self, file)

    def load(filename):
        with open(filename,"rb") as file:
            return pickle.load(file)

    def predict(self, amplicon, predictor):
        return self.models[amplicon].predict(predictor) 


    def observed_distribution(self, amplicon):
        return self.matrix[amplicon]

    def predicted_distribution(self, amplicon):
        other_amplicon = [i for i in self.amplicons if i != amplicon]
        predictor      = self.matrix[other_amplicon]
        return self.predict(amplicon,predictor)



    def diff_distribution(self, amplicon):
        other_amplicon = [i for i in self.amplicons if i != amplicon]
        target         = self.matrix[amplicon]
        predictor      = self.matrix[other_amplicon]
        y_predict      = self.predict(amplicon,predictor)
        y_observed     = target 
        return y_observed - y_predict


  