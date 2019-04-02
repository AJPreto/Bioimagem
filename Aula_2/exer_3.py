#!/usr/bin/env python

"""
This script calculates two-way anova
"""

import math
import os
import csv
import pandas as pd
import statistics
import scipy
from statsmodels.sandbox.stats.multicomp import multipletests
from scipy import stats
from statsmodels.stats.anova import anova_lm

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

def calculate_SSI(input_samples, input_averages):

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
prov = opened_file["Proveniencia"]
peso = opened_file["Peso"]
parto = opened_file["Parto"]

"""
Define groups for "Proveniencia"
"""

grupo_1_prov = peso.loc[prov == 1]
grupo_2_prov = peso.loc[prov == 2]
grupo_3_prov = peso.loc[prov == 3]
grupos_prov = [grupo_1_prov.tolist(), grupo_2_prov.tolist(), grupo_3_prov.tolist()]

"""
Define groups for "Parto"
"""

grupo_1_parto = peso.loc[parto == 1]
grupo_2_parto = peso.loc[parto == 2]
grupos_parto = [grupo_1_parto.tolist(), grupo_2_parto.tolist()]

"""
Calcular graus de liberdade
"""
amostras = 30
gl_prov = len(prov.unique()) - 1
gl_parto = len(parto.unique()) - 1
gl_interaction = gl_prov*gl_parto
gl_interaction_E = amostras - (len(prov.unique())*len(parto.unique()))

"""
Calculate overal average and group averages
"""

media_global = statistics.mean(peso)
media_grupo_1_prov = statistics.mean(grupo_1_prov)
media_grupo_2_prov = statistics.mean(grupo_2_prov)
media_grupo_3_prov = statistics.mean(grupo_3_prov)
medias_grupos_prov = [media_grupo_1_prov, media_grupo_2_prov, media_grupo_3_prov]
media_grupo_1_parto = statistics.mean(grupo_1_parto)
media_grupo_2_parto = statistics.mean(grupo_2_parto)
medias_grupos_parto = [media_grupo_1_parto, media_grupo_2_parto]

"""
Calculate the squared sum factor and error in order to calculate
 the mean-squared factor and errors
"""
SSF_prov = calculate_SSF(medias_grupos_prov, media_global, grupos_prov)
SSE_prov = calculate_SSE(grupos_prov, medias_grupos_prov)
MSF_prov = calculate_MS(SSF_prov, gl_prov)

SSF_parto = calculate_SSF(medias_grupos_parto, media_global, grupos_parto)
SSE_parto = calculate_SSE(grupos_parto, medias_grupos_parto)
MSF_parto = calculate_MS(SSF_parto, gl_parto)

"""
Calculate interaction 
"""
SST = sum((peso - media_global)**2)
prov_1 = opened_file[opened_file["Proveniencia"] == 1]
prov_2 = opened_file[opened_file["Proveniencia"] == 2]
prov_3 = opened_file[opened_file["Proveniencia"] == 3]
prov_1_means = [prov_1[prov_1["Parto"] == entry].Peso.mean() for entry in prov_1["Parto"]]
prov_2_means = [prov_2[prov_2["Parto"] == entry].Peso.mean() for entry in prov_2["Parto"]]
prov_3_means = [prov_3[prov_3["Parto"] == entry].Peso.mean() for entry in prov_3["Parto"]]
SSI = sum((prov_1.Peso - prov_1_means)**2) +sum((prov_2.Peso - prov_2_means)**2) + sum((prov_3.Peso - prov_3_means)**2)
SSF_both = SST - SSF_prov - SSF_parto - SSI
MSF_interaction = calculate_MS(SSF_both, gl_interaction)
MSE_interaction = calculate_MS(SSI, gl_interaction_E)

"""
Calculate factors
"""
primeiro_fator = MSF_prov / MSE_interaction
segundo_fator = MSF_parto / MSE_interaction
inter_fator = MSF_interaction / MSE_interaction

"""
Print desired outputs
"""
p_fator_1 = stats.f.sf(primeiro_fator, gl_prov, gl_interaction_E)
p_fator_2 = stats.f.sf(segundo_fator, gl_parto, gl_interaction_E)
p_inter_fator = stats.f.sf(inter_fator, gl_interaction, gl_interaction_E)

opened_output_file = open("two_way_anova.csv","w")
header = ",Squared Sum, Mean-Squared Sum, Degrees of freedom, f-factor, p-value\n"
row_1 = "Proveniencia," + str(round(SSF_prov,2)) + "," + str(round(MSF_prov,2)) + "," + str(gl_prov) + "," + str(round(primeiro_fator,2)) + "," + str(p_fator_1) + "\n"
row_2 = "Proveniencia," + str(round(SSF_parto, 2)) + "," + str(round(MSF_parto,2)) + "," + str(gl_parto) + "," + str(round(segundo_fator,2)) + "," + str(p_fator_2) + "\n"
row_3 = "Interacao," + str(round(SSF_both, 2)) + "," + str(round(MSF_interaction,2)) + "," + str(gl_interaction) + "," + str(round(inter_fator,2)) + "," + str(p_inter_fator) + "\n"
row_4 = "Residual," + str(round(SSI,2)) + "," + str(round(MSE_interaction,2)) + "," + str(gl_interaction_E) + "\n"
opened_output_file.write(header)
opened_output_file.write(row_1)
opened_output_file.write(row_2)
opened_output_file.write(row_3)
opened_output_file.write(row_4)

anova_formula = 'Peso ~ C(Proveniencia) + C(Parto) + C(Proveniencia):C(Parto)'
model = ols(anova_formula, opened_file).fit()
anova_table = anova_lm(model, typ=2)
print(anova_table)

