#!/usr/bin/env python

"""
This script performs linear regression over a csv file
containing two columns (x and y)
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
from scipy.stats import f
from scipy.stats import t

__author__ = "A.J. Preto"
__email__ = "martinsgomes.jose@gmail.com"
__group__ = "Data-Driven Molecular Design"
__group_leader__ = "Irina S. Moreira"
__project__ = "Biomedical Imaging"

def read_csv(csv_name, delimiter = ","):

	"""
	Read the input csv file
	"""
	opened_csv = open(csv_name, "r").readlines()
	for row in opened_csv[1:]:
		row = row.replace("\n","")
		new_row = row.split(delimiter)
		yield new_row

def calculate_average(input_table, input_col):

	"""
	Calculate the average over one column.
	Define "input_col" to define the column.
	"""
	average = 0
	for field in input_table:
		average += float(field[input_col])
	return average / len(input_table)

def calculate_SS(input_table, input_col, input_average):

	"""
	Calculate the Sum square over a column of a table, with the average
	"""
	SS = 0
	for field in input_table:
		SS += (input_average - float(field[input_col])) ** 2
	return SS

def calculate_SS_single(input_row, input_average):

	"""
	Calculate the Sum Squared over a list, with the average
	"""
	SS = 0
	for field in input_row:
		SS += (input_average - float(field)) ** 2
	return SS

def calculate_SSxy(input_table, average_x, average_y):

	"""
	Calculate the Sum Squared over two simultaneous columns,
	considering their respective averages
	"""
	output_val = 0
	for field in input_table:
		output_val += (float(field[0]) - float(average_x)) * (float(field[1]) - float(average_y))
	return output_val

def estimate_y(input_table, input_col, slope, origin):

	"""
	Calculate an estimate of y on the basis of the equation:
	estimate_y = m*x + b
	"""
	estimated_y = []
	for field in input_table:
		y_val = slope * float(field[input_col]) + origin
		estimated_y.append(y_val)
	return estimated_y

def SS_res(estimated_y, input_table, input_col):

	"""
	Calculate the residual Sum Square
	"""
	output_val = 0
	for y_1, y in zip(estimated_y, input_table):
		output_val += (float(y_1) - float(y[input_col])) ** 2
	return output_val


def degrees_of_freedom(input_table):

	"""
	Yields degrees of freedom
	"""
	return len(input_table) - 2


def f_statistics(SS_reg, SS_res, df):

	"""
	Calculate the f statistics
	"""
	F = ((SS_reg) / 1.0) / (SS_res / (df))
	return F

def split_table(input_table):

	"""
	Split the table, useful for some further calculations
	"""
	x = []
	y = []
	for row in input_table:
		x.append(float(row[0]))
		y.append(float(row[1]))
	return x, y

def f_test(F, df1, df2):

	"""
	Calculate the p value over an f value
	"""
	p_value = f.cdf(F, df1, df2)
	return p_value

def t_test(F, df):

	"""
	Calculate the p value over a t value
	"""
	p_value = t.cdf(F, df)
	return p_value

def t_determination(slope, SS_res, df, SSx):

	"""
	Calculate the t statistic of a distribution
	"""
	SS_res_factor = (float(SS_res)/float(df)) ** (0.5)
	SSx_factor = float(SSx) ** (0.5)
	first_div = (float(slope)/float(SS_res_factor))*float(SSx_factor)
	return first_div

def R_squared(SS_res, SS_reg):

	"""
	Calculate correlation
	"""
	SS_total = SS_res + SS_reg
	return float(SS_reg) / float(SS_total)


"""
Input data, change the name of the table to apply over any
 two-column table with a header, might need to change delimiter.
The first column must be x and the second column must be the y
values
"""
start_table = list(read_csv("aula_4.csv", delimiter = ";"))
x, y = split_table(start_table)

"""
Calculate the averages, sum squares, slope and b value.
Estimate y
"""
average_x = calculate_average(start_table, 0)
average_y = calculate_average(start_table, 1)
SSx = calculate_SS(start_table, 0, average_x)
SSy = calculate_SS(start_table, 1, average_y)
SSxy = calculate_SSxy(start_table, average_x, average_y)
m = SSxy/SSx
b = average_y - m * average_x
estimated_y = estimate_y(start_table, 0, m, b)

"""
Calculate sum squared of regression and residue
and degrees of freedom
"""
SS_reg = calculate_SS_single(estimated_y, average_y)
df = degrees_of_freedom(start_table)
SS_res = SS_res(estimated_y, start_table, 1)

"""
Calculate the f statistics, t statistics,
the corresponding p-values and the correlation
"""
F = f_statistics(SS_reg, SS_res, df)
p_value_f = 1 - f_test(F, 1, df)
t_value = t_determination(m, SS_res, df, SSx)
p_value_t = t_test(t_value, df) * 2
coef_det = R_squared(SS_res, SS_reg)


