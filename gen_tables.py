import pandas as pd
import re
import numpy as np
from collections import Counter


def mutDicts(muts,ref):
    mutStore = {}
    for m,i in zip(Counter(muts).keys(),Counter(muts).values()):
        if m != "-" and m != ref and m != "N":
            mutStore[m] = i
    return mutStore

inPile = pd.read_csv('./ncov_global.mpileup',sep='\t',header=None)
rawCalls = list(inPile.iloc[:,4].values)
rawSamples = list(inPile.iloc[:,6].values)
newCalls = []
newSamples = []

subs = [('[]^$]',''),("([.][-]\w*)",'.'),("([,][-]\w*)",',')]
for c in rawCalls:
    listOfCalls = []
    for o,r in subs:
        c = re.sub(o,r,c)
    for i in re.split("([.,][+]\w*)",c):
        if "+" in i:
            listOfCalls.append(i)
        else:
            for j in i:
                listOfCalls.append(j)
    if "" in listOfCalls:
        listOfCalls.remove("")

    newCalls.append(listOfCalls)
for s in rawSamples:
    n = s.split(',')
    newSamples.append(n)


cols = []
for l in newSamples:
    for col in l:
        if col not in cols:
            cols.append(str(col))

finTable = inPile.iloc[:,0:3]

callsList = []
indices = list(finTable.index)
refs = list(inPile.iloc[:,2].values)
for i,calls,samps,ref in zip(indices,newCalls,newSamples,refs):
    callsDict = {}
    for a,b in zip(calls,samps):
        allele = ""
        if a == '.' or a == ',':
            allele = ref
        else:
            allele = a
        callsDict[b] = allele
    callsList.append(callsDict)

testTable = pd.DataFrame.from_records(callsList,inPile.index.values)
final = finTable.join(testTable,how="outer")
final.fillna("-",inplace=True)
new_cols = final.columns.values
new_cols[0] = "Chromosome"
new_cols[1] = "Position"
new_cols[2] = "Reference Allele"
final.columns = new_cols
final.to_csv("./ncov_global_all.csv",index=False)
allelesOnly = final.iloc[:,3:]
allelesOnly.replace("N", "-", inplace=True)
allelesOnly["unique"] = allelesOnly.apply(lambda x: len(list(set([i for i in x.dropna().values if i != "-"]))),axis=1)
allelesOnly2 = allelesOnly[allelesOnly["unique"] > 1].drop("unique",axis=1)
outParsed = final.loc[allelesOnly2.index,:]
outParsed.to_csv("./ncov_global_parsed.csv",index=False)
finTable = outParsed.iloc[:,0:3]
refAlleles = outParsed["Reference Allele"]
sampAlleles = outParsed.iloc[:,2:]
finTable["Mut Call"] = sampAlleles.apply(lambda x: max(mutDicts(x.values,x["Reference Allele"]),key=mutDicts(x.values,x["Reference Allele"]).get),axis=1)
finTable.to_csv("./ncov_global_transformed.csv",index=False)
