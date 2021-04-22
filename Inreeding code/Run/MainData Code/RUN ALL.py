# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 15:19:33 2020

@author: z5035086
"""

import subprocess
import csv
import os
import pandas

path = os.getcwd()



with open("DATAsequence1.csv", "w", newline='') as fp:
                    wr = csv.writer(fp, dialect='excel')
                    #wr.writerow(FullList+SubList)
                    wr.writerow(['TreatmentGroup','Generations','AlleleDistrobution','PopulationSize','Number of alleles','RealNumberOfLoci','EstimatedNumberOfLociMEAN','EstimatedNumberOfLociMEDIAN','EstimatedNumberOfLociMAX','Fis','AverageShanPerIndi-1Hi','TotalShanOfPop-1Hs','Average1HiOver1Hs','AverageHighestAlleleCountPerIndividual','VarOfIndiSHAN', 'Number of Fixed loci','AverageRichPerIndi-0Hi','TotalRichOfPop-0Hs','VarOfIndiRICH','AverageHetPerIndi-2Hi','TotalHetOfPop-2Hs','VarOfIndiHET', 'EyePRime','Pee','Eye','Ell','Ohh'])

print('Starting data creation.')
print('0%')
subprocess.call(["python", "MSHomo 3 Loci.py"])
print('10.5%')
subprocess.call(["python", "MSHomo 5 Loci.py"])
print('16.5%')
subprocess.call(["python", "MSHomo 10 Loci.py"])
print('22%')
subprocess.call(["python", "RM 3 Loci.py"])
print('27.5%')
subprocess.call(["python", "RM 5 Loci.py"])
print('33%')
subprocess.call(["python", "RM 10 Loci.py"])
print('38.5%')
subprocess.call(["python", "SF 3 Loci.py"])
print('44%')
subprocess.call(["python", "SF 5 Loci.py"])
print('84ish I guess')
subprocess.call(["python", "SF 10 Loci.py"])
print('85%')
subprocess.call(["python", "SS 3 Loci.py"])
print('85.420%')
subprocess.call(["python", "SS 5 Loci.py"])
print('85.69%')
subprocess.call(["python", "SS 10 Loci.py"])
print('88%')
subprocess.call(["python", "MS 3 Loci.py"])
print('93.5%')
subprocess.call(["python", "MS 5 Loci.py"])
print('95%')
subprocess.call(["python", "MS 10 Loci.py"])
print('DONE!')
