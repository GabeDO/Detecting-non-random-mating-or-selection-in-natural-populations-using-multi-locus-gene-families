# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 15:19:33 2020

@author: z5035086
"""

import subprocess
import os, glob
import pandas as pd
import csv


path = os.getcwd()

with open("DATA.csv", "w", newline='') as fp:
                    wr = csv.writer(fp, dialect='excel')
                    #wr.writerow(FullList+SubList)
                    wr.writerow(['TreatmentGroup','Generations','AlleleDistrobution','PopulationSize','Number of alleles','RealNumberOfLoci','EstimatedNumberOfLoci','Fis','AverageShanPerIndi-1Hi','TotalShanOfPop-1Hs','Average1HiOver1Hs','AverageHighestAlleleCountPerIndividual','VarOfIndiSHAN', 'Number of Fixed loci','AverageRichPerIndi-0Hi','TotalRichOfPop-0Hs','VarOfIndiRICH','AverageHetPerIndi-2Hi','TotalHetOfPop-2Hs','VarOfIndiHET', 'EyePRime','Pee','Eye','Ell','Ohh'])



subprocess.call(["python", "MSHomo 3 Loci.py"])
print('69%, NICE')
subprocess.call(["python", "MSHomo 5 Loci.py"])
print('71.5%')
subprocess.call(["python", "MSHomo 10 Loci.py"])
print('77%')
subprocess.call(["python", "SSHomo 3 Loci.py"])
print('82.5%')
subprocess.call(["python", "SSHomo 5 Loci.py"])
print('88%')
subprocess.call(["python", "SSHomo 10 Loci.py"])
print('93.5%')
print('DONE!')
