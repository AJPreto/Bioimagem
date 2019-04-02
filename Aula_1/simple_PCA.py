#!/usr/bin/env python

"""
This script performs Principal Component Analysis and retrieves information on the components
The dataset used was retrieved from "winequality_red.csv"
"""

__author__ = "A.J. Preto"
__email__ = "martinsgomes.jose@gmail.com"
__group__ = "Data-Driven Molecular Design"
__group_leader__ = "Irina S. Moreira"
__project__ = "Biomedical Imaging"

import os
import pandas as pd
import math
import sklearn

def PCA(input_data, header, retrieve_comp = True):
    
    """
    Use this function to perform PCA for dimensionality reduction.
    The output dataset will be reduced to the amount of dimensions you have defined on the end of the file
    """
    from sklearn.decomposition import PCA
    pca = PCA(n_components = PCA_components)
    if retrieve_comp == True:
        pca.fit(input_data)
        retrieve_PCA_components(input_data, header, pca)
    return pca.fit_transform(input_data)

def retrieve_PCA_components(input_array, header, current_pca):

    """
    This function can be called from the PCA function in order to write a new file
    in which the contributions of each feature to eah component are written.
    """
    PC_list = []
    for PC_lenght in range(PCA_components):
        current_PC = "PC-" + str(PC_lenght)
        PC_list.append(current_PC)
    components_data = pd.DataFrame(current_pca.components_, columns = pd.DataFrame(data = input_array,
                          columns = list(input_array)).columns,index = PC_list)
    components_data = pd.DataFrame.transpose(components_data)
    components_data.to_csv("PCA_results.csv", sep=',')

PCA_components = 11
opened_file = pd.read_csv("winequality_red.csv", header = 0, sep = ";")
no_class_table = opened_file.drop("quality", axis = 1)
PCA(no_class_table, list(opened_file))







