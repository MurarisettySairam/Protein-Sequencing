"""
Protein Sequencing Project
Name:
Roll Number:
"""

from email.headerregistry import ContentDispositionHeader
from itertools import combinations
from posixpath import split
from typing import NewType
import hw6_protein_tests as test

project = "Protein" # don't edit this

### WEEK 1 ###

'''
readFile(filename)
#1 [Check6-1]
Parameters: str
Returns: str
'''
def readFile(filename):
    o=open(filename,"r").read()
    sl=''
    for i in o.splitlines():
        sl+=i
    return sl

    
'''
dnaToRna(dna, startIndex)
#2 [Check6-1]
Parameters: str ; int
Returns: list of strs
'''
def dnaToRna(dna, startIndex):
    l=[]
    s=dna.replace("T","U")
    for i in range(startIndex,len(s),3):
        l.append(s[i:i+3])
    for i in l:
        if i=="UAA" or i=="UGA" or i=="UAG":
            n=l.index(i)
            return l[:n+1]
    return l

'''
makeCodonDictionary(filename)
#3 [Check6-1]
Parameters: str
Returns: dict mapping strs to strs
'''
def makeCodonDictionary(filename):
    import json
    d1={}
    o=open(filename,"r")
    d=json.load(o)
    for i in d:
        for j in d[i]:
            j=j.replace('T','U')
            d1[j]=i
    return d1


'''
generateProtein(codons, codonD)
#4 [Check6-1]
Parameters: list of strs ; dict mapping strs to strs
Returns: list of strs
'''
def generateProtein(codons, codonD):
    L=[]
    for i in codons:
        for j in codonD:
            if i==j:
                L.append(codonD[j])
                L[0]="Start"
    return L


'''
synthesizeProteins(dnaFilename, codonFilename)
#5 [Check6-1]
Parameters: str ; str
Returns: 2D list of strs
'''
def synthesizeProteins(dnaFilename, codonFilename):
    list1=[]
    unused=0
    f= readFile(dnaFilename)
    s= makeCodonDictionary(codonFilename)
    i=0

    while i < len(f):
        n= f[i:i+3]
        if n == 'ATG':
            nca=dnaToRna(f,i)
            new_call= generateProtein(nca,s)
            list1.append(new_call)
            i+=3*len(nca)
        else:
            i=i+1
            unused=unused+1 
    return list1


def runWeek1():
    print("Human DNA")
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    print("Elephant DNA")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")


### WEEK 2 ###

'''
commonProteins(proteinList1, proteinList2)
#1 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs
Returns: 2D list of strs
'''
def commonProteins(proteinList1, proteinList2):
    l=[]
    for i in proteinList1:
        if i in proteinList2:
            l.append(i)
    return l
'''
combineProteins(proteinList)
#2 [Check6-2]
Parameters: 2D list of strs
Returns: list of strs
'''
def combineProteins(proteinList):
    l=[]
    for i in proteinList:
        for j in i:
            l.append(j)
    return l 


'''
aminoAcidDictionary(aaList)
#3 [Check6-2]
Parameters: list of strs
Returns: dict mapping strs to ints
'''
def aminoAcidDictionary(aaList):
    d={}
    for i in aaList:
        if i not in d:
            d[i]=aaList.count(i)   
    return d


'''
findAminoAcidDifferences(proteinList1, proteinList2, cutoff)
#4 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs ; float
Returns: 2D list of values
'''
def findAminoAcidDifferences(proteinList1, proteinList2, cutoff):
    n1=combineProteins(proteinList1)
    n2=combineProteins(proteinList2)
    n3=aminoAcidDictionary(n1)
    n4=aminoAcidDictionary(n2)
    l=[]
    for i,j in n3.items():
        if i!="Start" and i!="Stop" and i not in l:
            l.append(i)
    for i,j in n4.items():
        if i!="Start" and i!="Stop" and i not in l:
            l.append(i)
    f=[]
    for i in l:
        frequency_1=0
        frequency_2=0
        if i in n1:
            frequency_1=n3[i]/len(n1)
        if i in n2:
            frequency_2=n4[i]/len(n2)
        diff=frequency_1-frequency_2
        if diff>cutoff or diff<-cutoff:
            f.append([i,frequency_1,frequency_2])
    return f


'''
displayTextResults(commonalities, differences)
#5 [Check6-2]
Parameters: 2D list of strs ; 2D list of values
Returns: None
'''
def displayTextResults(commonalities, differences):
    l=[]
    for i in commonalities:
        i.remove('Start')
        i.remove('Stop')
        if len(i)>1:
            symbol='-'.join(i)
            l.append([symbol])
        else:
            if i not in l:
                l.append(i)
    s=sorted(l)
    for i in s:
        for j in i:
            print(j)
    for i in differences:
        print(i[0]+":"+str(round(i[1]*100,2))+"%"+ " in seq1 ,"+str(round(i[2]*100,2))+"%"+"in seq2")
    return 


def runWeek2():
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")

    commonalities = commonProteins(humanProteins, elephantProteins)
    differences = findAminoAcidDifferences(humanProteins, elephantProteins, 0.005)
    displayTextResults(commonalities, differences)


### WEEK 3 ###

'''
makeAminoAcidLabels(proteinList1, proteinList2)
#2 [Hw6]
Parameters: 2D list of strs ; 2D list of strs
Returns: list of strs
'''
def makeAminoAcidLabels(proteinList1, proteinList2):
    l=[]
    c1=combineProteins(proteinList1)
    c2=combineProteins(proteinList2)
    c=c1+c2
    for i in c:
        if i not in l:
            l.append(i)
    s=sorted(l)
    return s


'''
setupChartData(labels, proteinList)
#3 [Hw6]
Parameters: list of strs ; 2D list of strs
Returns: list of floats
'''
def setupChartData(labels, proteinList):
    return


'''
createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None)
#4 [Hw6] & #5 [Hw6]
Parameters: list of strs ; list of floats ; str ; list of floats ; str ; [optional] list of strs
Returns: None
'''
def createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None):
    import matplotlib.pyplot as plt
    return


'''
makeEdgeList(labels, biggestDiffs)
#5 [Hw6]
Parameters: list of strs ; 2D list of values
Returns: list of strs
'''
def makeEdgeList(labels, biggestDiffs):
    return


'''
runFullProgram()
#6 [Hw6]
Parameters: no parameters
Returns: None
'''
def runFullProgram():
    return


### RUN CODE ###

# This code runs the test cases to check your work
if __name__ == "__main__":
    # print("\n" + "#"*15 + " WEEK 1 TESTS " +  "#" * 16 + "\n")
    # test.week1Tests()
    # print("\n" + "#"*15 + " WEEK 1 OUTPUT " + "#" * 15 + "\n")
    # runWeek1()
    # test.testReadFile()
    # test.testDnaToRna()
    # test.testMakeCodonDictionary()
    # test.testGenerateProtein()
    # test.testSynthesizeProteins()
    # test.testCommonProteins()
    # test.testCombineProteins()
    # test.testAminoAcidDictionary()
    # test.testFindAminoAcidDifferences()
    # runWeek2()
    test.testMakeAminoAcidLabels()
    ## Uncomment these for Week 2 ##
    """
    print("\n" + "#"*15 + " WEEK 2 TESTS " +  "#" * 16 + "\n")
    test.week2Tests()
    print("\n" + "#"*15 + " WEEK 2 OUTPUT " + "#" * 15 + "\n")
    runWeek2()
    """

    ## Uncomment these for Week 3 ##
    """
    print("\n" + "#"*15 + " WEEK 3 TESTS " +  "#" * 16 + "\n")
    test.week3Tests()
    print("\n" + "#"*15 + " WEEK 3 OUTPUT " + "#" * 15 + "\n")
    runFullProgram()
    """
