#!/usr/bin/env python

"""
This script performs simple statiscal analysis retrieving: mean, standard 
deviation, variance, maximum and minimum values on a dataset column wise. 
The dataset used was retrieved from "winequality_red.csv"
"""

__author__ = "A.J. Preto"
__email__ = "martinsgomes.jose@gmail.com"
__group__ = "Data-Driven Molecular Design"
__group_leader__ = "Irina S. Moreira"
__project__ = "Biomedical Imaging"

import os
import pandas as pd
import numpy as np
import math
import statistics

def process_csv(input_file, delimiter = ","):

	"""
	Simple generator to open an process the input csv file
	"""
	opened_file = open(input_file, "r").readlines()
	for row in opened_file:
		row = row.replace("\n","")
		row = row.replace('"','')
		new_row = row.split(delimiter)
		yield new_row

def write_csv_mean(input_data, header, output_name ,delimiter = ",", decimal_houses = 2):

	"""
	Open an output file, calculate the mean, meadian, standard deviation, variance, maximum
	and minimum values
	"""
	output_file = open(output_name, "w")
	output_file.write("attribute,mean,median,standard deviation, variance,maximum,minimum\n")
	for cols in range(len(header)):
		target_col = processed_file.iloc[:,cols]
		mean = round(statistics.mean(target_col), decimal_houses)
		median = round(statistics.median(target_col), decimal_houses)
		sd = round(statistics.stdev(target_col), decimal_houses)
		var = round(statistics.variance(target_col), decimal_houses)
		max_value = target_col.max()
		min_value = target_col.min()
		to_write_row = str(header[cols])+","+str(mean)+","+str(median)+","+str(sd)+","+str(var)+","+str(max_value)+","+str(min_value)+"\n"
		output_file.write(to_write_row)

target_file = "winequality_red.csv"
processed_file = pd.read_csv(target_file, sep = ";", header = 0)
header = list(processed_file)
write_csv_mean(processed_file, header, "wine_analysis.csv")


