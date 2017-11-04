#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 17:19:03 2017

@author: wandermol
"""

import sys
import os
import re
import math

# Invocation: python3 scriptname reference_molecule(s) database coefficient knn output
#
# reference_molecule(s) and database are expected to be in the format:
# mol1 25 32 33 42 51 52 53 55 58 59 60 61 64 65 67 69 73 75 [...]       
# mol2 25 52 56 62 63 65 70 71 77 78 80 83 94 96 98 99 102 105 [...]
#
# coefficient must one of the following similarity measures: 
# tanimoto, dice, cosine (Bajusz et al. 2015)
# or the distance measure:
# soergel (Bajusz et al. 2015)
#
# knn is a number k as in k-nn
# Open reference molecules and database. Sets coefficient, k-NN, and output file 

ref = open(sys.argv[1], "r")
db = open(sys.argv[2], "r")
coeff = sys.argv[3]
nn = int(sys.argv[4])
rank = open(sys.argv[5], "w")
ranktemp = open("ranktemp.tsv", "w+")

# Define Tanimoto coefficient

def tanimoto(mol1, mol2):

    fp1 = re.findall("\t(.+)", mol1)[0].split()
    a = len(fp1)

    fp2 = re.findall("\t(.+)", mol2)[0].split()
    b = len(fp2)

    c = 0
    for bit in fp1:
        if bit in fp2:
            c = c + 1

    tc = c / (a + b - c)

    return(tc)

# Define dice coefficient

def dice(mol1, mol2):

    fp1 = re.findall("\t(.+)", mol1)[0].split()
    a = len(fp1)

    fp2 = re.findall("\t(.+)", mol2)[0].split()
    b = len(fp2)

    c = 0
    for bit in fp1:
        if bit in fp2:
            c = c + 1

    dc = (2 * c) / (a + b)

    return(dc)

# Define cosine coeffient

def cosine(mol1, mol2):

    fp1 = re.findall("\t(.+)", mol1)[0].split()
    a = len(fp1)

    fp2 = re.findall("\t(.+)", mol2)[0].split()
    b = len(fp2)

    c = 0
    for bit in fp1:
        if bit in fp2:
            c = c + 1

    cc = c / math.sqrt(a + b)

    return(cc)

# Run the virtual screening

for dbmol in db:
    if len(dbmol.split()) < 2:
        continue
    name1 = re.findall("(.+)\t", dbmol)[0]
    simraw = [name1]
    for refmol in ref:
        if len(refmol.split()) < 2:
            continue
        if coeff == "tanimoto":
            simraw.append(tanimoto(dbmol, refmol))
        elif coeff == "dice":
            simraw.append(dice(dbmol, refmol))
        elif coeff == "cosine":
            simraw.append(cosine(dbmol, refmol))
        elif coeff == "soergel":
            simraw.append(1 - tanimoto(dbmol, refmol))
    ref.seek(0)
    if nn == 1:
        sim = max(simraw[1:])
        ranktemp.write(name1 + "\t" + str(sim) + "\n")
    else:
        sim = (sum(sorted(simraw[1:], reverse=True)[0:nn])/nn)
        ranktemp.write(name1 + "\t" + str(sim) + "\n")
ranktemp.seek(0)

# Write ranking

if coeff != "soergel":
    for line in sorted(ranktemp, key=lambda line: line.split()[1], reverse=True):
        rank.write(line)
else:
    for line in sorted(ranktemp, key=lambda line: line.split()[1]):
        rank.write(line)

# Remove temporary file

os.remove("ranktemp.tsv")                           
