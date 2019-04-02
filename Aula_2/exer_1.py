#!/usr/bin/env python

"""
This script performs ANOVA analysis
"""

import math
import os
import csv
import pandas as pd
import statistics
import scipy
from statsmodels.sandbox.stats.multicomp import multipletests

__author__ = "A.J. Preto"
__email__ = "martinsgomes.jose@gmail.com"
__group__ = "Data-Driven Molecular Design"
__group_leader__ = "Irina S. Moreira"
__project__ = "Biomedical Imaging"

def calculate_SSF(input_averages, global_average, input_groups):

	"""
	Calculate the squared sum factor
	"""

	output = 0
	for value, sample in zip(input_averages, input_groups):
		factor = (value - global_average)**2
		output += factor * len(sample)
	return output

def calculate_SSE(input_samples, input_averages):

	"""
	Calculate the squared sum error
	"""
	output = 0
	for group, group_average in zip(input_samples, input_averages):
		new_value = 0
		for sample in group:
			new_value += (sample - group_average) ** 2
		output += new_value
	return output

def calculate_factor(MSF, MSE):

	"""
	Calculate the factor
	"""
	return MSF/MSE

def calculate_MS(input_M, input_gl):

	"""
	Calculate the Mean Square
	"""
	return float(input_M)/float(input_gl)

def t_test(input_x, input_y):

	"""
	Perform t-test for method confirmation
	"""
	stat, p_value = scipy.stats.ttest_ind(input_x, input_y)
	return stat, p_value

def anova(*args):

	"""
	Perform anova and retrieve tau and p-values
	Note: this anova accepts an arbitrary number of groups
	"""
	tau, p_value = scipy.stats.f_oneway(*args)
	return tau, p_value

"""
Define input variable columns and groups
"""

opened_file = pd.read_csv("exer_1.csv", header = 0, sep = "\t")

parto = opened_file["Parto"]
peso = opened_file["Peso"]

grupo_1 = peso.loc[parto == 1]
grupo_2 = peso.loc[parto == 2]

grupos = [grupo_1.tolist(), grupo_2.tolist()]

"""
Calculate overal average and group averages
"""
media_global = statistics.mean(peso)
media_grupo_1 = statistics.mean(grupo_1)
media_grupo_2 = statistics.mean(grupo_2)
medias_grupos = [media_grupo_1, media_grupo_2]

"""
Calculate the squared sum factor and error in order to calculate
 the mean-squared factor and errors
"""
SSF = calculate_SSF(medias_grupos, media_global, grupos)
SSE = calculate_SSE(grupos, medias_grupos)
MSF = calculate_MS(SSF, 1)
MSE = calculate_MS(SSE, 28)

"""
Print desired outputs
"""
print(SSF)
print("mean-squared factors:",MSF)
print("mean-squared errors:",MSE)
print("factors",calculate_factor(MSF, MSE))
print("p-value:",anova(grupo_1, grupo_2)[1])
