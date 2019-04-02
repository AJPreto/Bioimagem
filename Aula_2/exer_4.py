#!/usr/bin/env python

"""
This script performs ANCOVA analysis
"""

from statsmodels.formula.api import ols

__author__ = "A.J. Preto"
__email__ = "martinsgomes.jose@gmail.com"
__group__ = "Data-Driven Molecular Design"
__group_leader__ = "Irina S. Moreira"
__project__ = "Biomedical Imaging"

opened_file = pd.read_csv("exer_1.csv", header = 0, sep = "\t")
prov = opened_file["Proveniencia"]
peso = opened_file["Peso"]
parto = opened_file["Parto"]
ancova_formula = "Peso ~ Proveniencia * Parto"
lm = ols(ancova_formula, opened_file).fit()
print(lm.summary())
