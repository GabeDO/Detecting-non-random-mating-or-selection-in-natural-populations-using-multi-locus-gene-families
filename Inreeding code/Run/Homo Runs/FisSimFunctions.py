# -*- coding: utf-8 -*-
"""
Created on Wed Sep 12 16:36:43 2018

@author: Gabe O'Reilly
"""
def AbundData(pop, pop_size, Number_of_Alleles):
    count_geno_data = []
    ya_boy =[]
    geno_data= []
    
    for i in range(pop_size):    
        ind = pop.individual(i)
        geno_data.append([ind.genotype(0) ,ind.genotype(1)])

    for i in range(len(geno_data)):
        for l in range(Number_of_Alleles):
            C = l
            if l < Number_of_Alleles-1:
                ya_boy.append(geno_data[i][0].count(C) +  geno_data[i][1].count(C))
            else:
                ya_boy.append(geno_data[i][0].count(C) +  geno_data[i][1].count(C))
                count_geno_data.append(ya_boy)
                ya_boy=[]
    #this will convert the data to presence absence data
    return count_geno_data

def PresAbsData(count_geno_data):
    for i in range(len(count_geno_data)):
        for l in range(len(count_geno_data[i])):
           if count_geno_data[i][l] > 0:
                count_geno_data[i][l] = 1 

    return count_geno_data

def PresAbsDataOG(pop, pop_size, Number_of_Alleles):
    count_geno_data = []
    ya_boy =[]
    geno_data= []
    
    for i in range(pop_size):    
        ind = pop.individual(i)
        geno_data.append([ind.genotype(0) ,ind.genotype(1)])
    
    for i in range(len(geno_data)):
        for l in range(Number_of_Alleles):
            C = l
            if l < Number_of_Alleles-1:
                ya_boy.append(geno_data[i][0].count(C) +  geno_data[i][1].count(C))
            else:
                ya_boy.append(geno_data[i][0].count(C) +  geno_data[i][1].count(C))
                count_geno_data.append(ya_boy)
                ya_boy=[]
    #this will convert the data to presence absence data
    return count_geno_data
    for i in range(len(count_geno_data)):
        for l in range(len(count_geno_data[i])):
           if count_geno_data[i][l] > 0:
                count_geno_data[i][l] = 1 
   
    #return count_geno_data

def Selection4Simupop(loci_num, mating):


    if mating == 'inbred' or mating == 'i':
        X = 1
        Y = 0
        print("selection set for inbreeding")
    elif mating == 'outbred'or mating == 'o':
        X = 0
        Y = 1
        print("selection set for outbreeding")
    else:
        X = 1
        Y = 1
        print("selection set for random mating")
    
    selection_dic ={}
    selection_list = []    
    
    for x in range(loci_num):
        for y in range(loci_num):
            selection_list.append([x,y])
    
    for i in range(len(selection_list)):
        if int(selection_list[i][0]) >= int(selection_list[i][1]):
            if selection_list[i][0] == selection_list[i][1]:
                selection_dic[str(selection_list[i])] = X
            else:
                selection_dic[str(selection_list[i])] = Y
        else:
            yeet = 1
    dic = str(selection_dic)
    
    dic = dic.replace(' ', '')
    dic = dic.replace('[', ' (')        
    dic = dic.replace(']', ')')
    dic = dic.replace("'", '')
    return dic.strip('"\'')

def BillsMath(count_geno_data):
    
    import math
    import numpy as np
    #making it a matrix, then summing along each coloum. Removing all 0s, THEN looking at LEN. deals with lost alleles
    Matrix_count_geno_data = np.matrix(count_geno_data)
    SumCGD = Matrix_count_geno_data.sum(axis=0)
    ArrayCGD = np.squeeze(np.asarray(SumCGD))
    NoZeroMCGD = list(filter(lambda a: a != 0, ArrayCGD))
    P = len(NoZeroMCGD)
    
    data = count_geno_data
    #print(ArrayCGD)
    #for ever row in data (each individual)...
    for i in range(len(count_geno_data)):
        #sum up the total number of different varients for each individual  
        data[i] = sum(map(float,count_geno_data[i]))    

    #print(data)
    I = sum(data)/len(data)
    #if the max # of alleles is a decimal after being devided by 2, it's not even, so we add +1 to it to make it even. Otherwise, continue as normal
    if (max(count_geno_data)/2).is_integer() is True:
        L = max(count_geno_data) /2
    else:
        L = (max(count_geno_data)+1) /2
    
    O = I/P
    if int(L) == 0:
        Eye_Prime = 0
    else:
        Eye_Prime = (((3*L) - (2*I))/(L))
     
    
    #print("pee", P)
    #print("eye", I)
    #print("ell", L)
    #print("Oooh", O)
    #print("Eye_Prime",  Eye_Prime)
    
    return [Eye_Prime, str(P), str(I), L, str(O)]

#file = 'genepopfile.txt.DIV'
def readGENOPOPfis(file):
    f=open(file)
    lines=f.readlines()
    line_len=len(lines) - 12
    items = [lines[line_len]]
    Fis_List = []
    for items in items:
        Fis_List.extend(items.split())
    print("Fis =", Fis_List[3])
  
  
def measureFIS(file, NumLoci):
    import csv
    import numpy
    
    def func(a):
        if  a[0] == a[-1]:
            return 2
        else:
            return 1
    
     
    #EXPECTED HET PER LOCI
    a = numpy.genfromtxt(file, delimiter=',', skip_header=True)[:, 1:]
    HetE = []
    HetO = []
    Fis = []
    for i in range(NumLoci):
        res = numpy.unique(a[:, 1+(i*2):3+(i*2)],return_counts=True)
        y=0

        for i in range(len(res[1])):
            x = (res[1][i] / sum(res[1]))**2
            y = x + y
        HetE.append(1-y)
    

    #print("Expected Het per Loci", HetE)
    #OBSERVED HET PER LOCI
    for i in range(NumLoci):
        b = numpy.apply_along_axis( func , 1, a[:,1+(i*2):3+(i*2)])
        HetO.append(list(b).count(1)/len(list(b)))
    #print("Observed Het per Loci", HetO)
    
    #FIS per Loci
    for i in range(NumLoci):
        if HetE[i] == 0:
            Fis.append(0)
        else:
            Fis.append((HetE[i] - HetO[i])/HetE[i])
        
    #print("Fis per Loci", Fis)
    AvgFis = sum(Fis)/len(Fis)
    #print("Average Fis=", AvgFis)
    
    return AvgFis

def QProfile(file, NumLoci, total):
    import csv
    import numpy
    import math
    
    def func(a):
        if  a[0] == a[-1]:
            return 2
        else:
            return 1
    
     
    #EXPECTED HET PER LOCI
    a = numpy.genfromtxt(file, delimiter=',', skip_header=True)[:, 1:]
    DHetE = []
    DShan = []
    DAlRich = []
    for i in range(NumLoci):
        res = numpy.unique(a[:, 1+(i*2):3+(i*2)],return_counts=True)
        y=0
    
        for i in range(len(res[1])):
            x = (res[1][i] / sum(res[1]))**2
            y = x + y
        DHetE.append(1/(1-(1-y)))

    
    for i in range(NumLoci):
        res = numpy.unique(a[:, 1+(i*2):3+(i*2)],return_counts=True)
        y=0

        for i in range(len(res[1])):
            x = -(res[1][i] / sum(res[1]))*math.log((res[1][i] / sum(res[1])))
            y = x + y
        DShan.append(math.exp(y))
        
    for i in range(NumLoci):
        res = numpy.unique(a[:, 1+(i*2):3+(i*2)],return_counts=True)
        DAlRich.append(len(res[1]))
    
    
    avgDHetE = sum(DHetE) / len(DHetE)
    avgDShan = sum(DShan) / len(DShan)
    avgDAlRich = sum(DAlRich) / len(DAlRich)
    
    if total == 'NO':    
        return [avgDAlRich, avgDShan, avgDHetE]
    if total == 'YES':
        return [DAlRich,DShan,DHetE]

def AverageRich(resultt):
    import math
    import numpy as np
    ListOfRich = []
    ListOfMax = []
    #print(resultt)
    for r in resultt:      
        Prop = []
        total =  sum(r)
        for l in r:
            p = l/total
            if p == 0:
                yeet = 0
            else:
                Shp = 1                
                Prop.append(Shp)
        ListOfRich.append(1-sum(Prop))
        ListOfMax.append(max(r))
        
    AverageRich = sum(ListOfRich)/len(ListOfRich)
    RichVar = np.var(ListOfRich)
    #MaxEntropy = math.log(Number of alleles)
    c = resultt
    TotalAlleles = [sum(i) for i in zip(*c)]
    TotalRichProps = []
    for h in TotalAlleles:
        Prop = []
        Total = sum(TotalAlleles)
        p = h/Total
        if p == 0:
            yeet = 0
        else:
            TSP = 1
            TotalRichProps.append(TSP)

    return [AverageRich, sum(TotalRichProps), RichVar]
    
def AverageShan(resultt):
    import math
    import numpy as np
    ListOfShan = []
    ListOfMax = []
    #print(resultt)
    for r in resultt:      
        Prop = []
        total =  sum(r)
        for l in r:
            p = l/total
            if p == 0:
                yeet = 0
            else:
                Shp = -p*(math.log(p))                
                Prop.append(Shp)
        ListOfShan.append(sum(Prop))
        ListOfMax.append(max(r))
    MaxAlVal = sum(ListOfMax)/len(ListOfMax)
    AverageShan = sum(ListOfShan)/len(ListOfShan)
    ShanVar = np.var(ListOfShan)
    #MaxEntropy = math.log(Number of alleles)
    c = resultt
    TotalAlleles = [sum(i) for i in zip(*c)]
    TotalShanProps = []
    for h in TotalAlleles:
        Prop = []
        Total = sum(TotalAlleles)
        p = h/Total
        if p == 0:
            yeet = 0
        else:
            TSP = -p*(math.log(p))
            TotalShanProps.append(TSP)
    if sum(TotalShanProps) == 0:
        LstAvgIndiOvrTotShan=  [x*0 for x in ListOfShan]
    else:
        LstAvgIndiOvrTotShan=  [x/sum(TotalShanProps) for x in ListOfShan]
    AvgIndiOvrTotShan = sum(LstAvgIndiOvrTotShan)/len(LstAvgIndiOvrTotShan)
    return [AverageShan, sum(TotalShanProps), AvgIndiOvrTotShan, MaxAlVal, ShanVar]

def AverageHet(resultt):
    import math
    import numpy as np
    ListOfHET = []
    ListOfMax = []
    #print(resultt)
    for r in resultt:      
        Prop = []
        total =  sum(r)
        for l in r:
            p = l/total
            if p == 0:
                yeet = 0
            else:
                Shp = p**2                
                Prop.append(Shp)
        ListOfHET.append(1-sum(Prop))
        ListOfMax.append(max(r))
    AverageHET = sum(ListOfHET)/len(ListOfHET)
    HetVar = np.var(ListOfHET)
    #MaxEntropy = math.log(Number of alleles)
    c = resultt
    TotalAlleles = [sum(i) for i in zip(*c)]
    TotalHETProps = []
    for h in TotalAlleles:
        Prop = []
        Total = sum(TotalAlleles)
        p = h/Total
        if p == 0:
            yeet = 0
        else:
            TSP = p**2
            TotalHETProps.append(TSP)

    return [AverageHET, (1-sum(TotalHETProps)), HetVar]

def measureStarIS(file, NumLoci):
    import csv
    import numpy
    
    def func(a):
        if  a[0] == a[-1]:
            return 2
        else:
            return 1
    
     
    #EXPECTED HET PER LOCI
    a = numpy.genfromtxt(file, delimiter=',', skip_header=True)[:, 1:]
    HetE = []
    HetO = []
    Staris = []
    for i in range(NumLoci):
        res = numpy.unique(a[:, 1+(i*2):3+(i*2)],return_counts=True)
        y=0

        for i in range(len(res[1])):
            x = (res[1][i] / sum(res[1]))**2
            y = x + y
        HetE.append(1-y)
    

    print("Expected Het per Loci", HetE)
    #OBSERVED HET PER LOCI
    for i in range(NumLoci):
        b = numpy.apply_along_axis( func , 1, a[:,1+(i*2):3+(i*2)])
        HetO.append(list(b).count(1)/len(list(b)))
    print("Observed Het per Loci", HetO)
    
    #FIS per Loci
    for i in range(NumLoci):
        if HetE[i] == 0:
            Staris.append(0)
        else:
            Staris.append((HetE[i] - HetO[i])/(1-HetO[i]))
        
    print("Staris per Loci", Staris)
    AvgStaris = sum(Staris)/len(Staris)
    print("Average Staris=", AvgStaris)
    
    return Staris
    
    
def FisDistroMath(a):    
    AlleleList = []
    for i in range(len(a)):
        IndiAlleleList = []
        for l in range(len(a[i])):
            if a[i][l] == 1:
                IndiAlleleList.append(l+1)              
        AlleleList.append(IndiAlleleList)
    
    
    return AlleleList   
        
def SequSampling(data, Depth):
    #bills matlab code to stochastically sample data, as a sequencer would
    #it runs over the results from each individual and randomly takes their value
    #and if ti does, add it to the list. 
    import numpy as np
    plist = []
    output = []
    total = len(data)
    TDepth = total*Depth
    for i in data:
        p = i/sum(data)
        plist.append(p)
    for i in plist:
        output.append(np.random.binomial(TDepth, i))
    return output

def GabesMinEesultNo0(results):
        Minresultno0 = [] 
        for i in results:
            x = list(filter(lambda a: a != 0, i))
            Minresultno0.append(min(x))
        return Minresultno0

def GabesLowestComDom(results, Minresultno0):
    resultsLowComDom = []
    for i in range(len(results)):    
        resultsLowComDom.append([x / Minresultno0[i] for x in results[i]])
    return resultsLowComDom


def GabesLocEstimate(results):

    Minresultno0 = GabesMinEesultNo0(results)
    resultsLowComDom = GabesLowestComDom(results, Minresultno0)
    
    summedResLoCoDo = []
    for i in resultsLowComDom:
        summedResLoCoDo.append(sum(i)) 
    summedResLoCoDo.sort() 
    summedResLoCoDono1 = list(filter(lambda a: a != 1, summedResLoCoDo)) 
    #print(summedResLoCoDono1)
    
    if summedResLoCoDono1 == []:
        return 'Fixed' 
    else:    
        return (summedResLoCoDono1[-1]/2)

####https://stackoverflow.com/questions/1624883/alternative-way-to-split-a-list-into-groups-of-n
def grouper(n, iterable, fillvalue=None):
    import itertools
    "grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return itertools.zip_longest(*args, fillvalue=fillvalue)

def Make_LoComDom_Data_full_with_loc_estimate(LowestResults, NumberOfLoci):
    newResult = []
    for i in LowestResults:
        if (sum(i) == (NumberOfLoci)*2):
            newResult.append(i)
        else:
            y = [ round(elem, 0) for elem in [(NumberOfLoci*2)/sum(i)*x for x in i]]
            newResult.append(y)
    return newResult

def Expanding_Result_Data_nonstacked(StackedData):
    #stacked data either right from output, or from Make_LoComDom_Data_full_with_loc_estimate()
    FullNewResult = []
    for i in StackedData:
        A = 1   
        makingNewResult = []
        for x in i:
            for y in range(int(x)):
                makingNewResult.append(A)
            A = A + 1
        FullNewResult.append(makingNewResult)
    return FullNewResult         

def makingPermdata(inputResults):
    import itertools
    permutationdata = []
    for i in range(len(inputResults)):
        ####all possible permutations
        k = [list(p) for p in itertools.permutations(inputResults[i])]
        #### Removing duplicate permutaions (not dealing with a perloci basis yet)
        k.sort()
        individual_nonloci_permutaions = list(k for k,_ in itertools.groupby(k))
        grouped_loci_individuals_not_filtered = []
        for l in range(len(individual_nonloci_permutaions)):
            x = list(grouper(2, individual_nonloci_permutaions[l]))
            turn_tups_to_list = []
            ###Ordering the loci in accending order so I can just remove functional duplicates (like [2,1] and [1,2])
            for p in x:
                m = sorted(p)
                turn_tups_to_list.append(m)
            grouped_loci_individuals_not_filtered.append(turn_tups_to_list)
        #### Removing functional duplicates duplicate loci permutaions
        grouped_loci_individuals_not_filtered.sort()
        grouped_loci_individuals_filtered = list(grouped_loci_individuals_not_filtered for grouped_loci_individuals_not_filtered,_ in itertools.groupby(grouped_loci_individuals_not_filtered))
        #print("YEET", grouped_loci_individuals_filtered)
        permutationdata.append(grouped_loci_individuals_not_filtered)
        print("Progress on data generation:",((i/len(inputResults))*100))
    
    return permutationdata

def makingPermdata2(inputResults,HowManyReps):
    import itertools
    import csv
    import random
    permutationdata = []
    k = []
    progress = 0
    f = open("permutationdata.csv", "w+")
    f.close()
    for i in range(len(inputResults)):
        ####all possible permutations
        for l in range(HowManyReps): 
            random.shuffle(inputResults[i])
            k.append(inputResults[i])
        #### Removing duplicate permutaions (not dealing with a perloci basis yet)
        k.sort()
        individual_nonloci_permutaions = list(k for k,_ in itertools.groupby(k))
        grouped_loci_individuals_not_filtered = []
        for l in range(len(individual_nonloci_permutaions)):
            x = list(grouper(2, individual_nonloci_permutaions[l]))            
            grouped_loci_individuals_not_filtered.append(x)
        #######################################################################
        grouped_loci_individuals_not_filtered.sort()
        grouped_loci_individuals_filtered = list(grouped_loci_individuals_not_filtered for grouped_loci_individuals_not_filtered,_ in itertools.groupby(grouped_loci_individuals_not_filtered))
        permutationdata.append(grouped_loci_individuals_filtered)
        #######################################################################
        print("Progress on data generation:",((i/len(inputResults))*100))

    
    return permutationdata


#https://www.edureka.co/community/1245/splitting-a-list-into-chunks-in-python
def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i + n]  

#https://stackoverflow.com/questions/7207309/python-how-can-i-run-python-functions-in-parallel      
def runInParallel(*fns):
  from multiprocessing import Process
  proc = []
  for fn in fns:
    p = Process(target=fn)
    p.start()
    proc.append(p)
  for p in proc:
    p.join()
    
def CountFixedLoci(pop,NumerOfAlleles,LociNumber):
    
    FixedLoci = []
    print(LociNumber, NumerOfAlleles)
    for i in range(LociNumber):
        ListofalleleFreq = []
        for A in range(NumerOfAlleles):
                ListofalleleFreq.append(pop.dvars().alleleFreq[i][A])
        if max(ListofalleleFreq) == 1:
           FixedLoci.append(1)
        else:
            FixedLoci.append(0)    
    
    return sum(FixedLoci) 