'''
Created on Mar 4, 2014

@author: eotles
'''

from hw3_shared import *

#run tests                    
def runTests(testData, outFile, posModel, negModel, bases):
    for seq in testData:
        outFile.write(str(score(seq, posModel, negModel, bases))+"\n")

def main():     
    #bases
    bases = ['A','C','G','T']
    #get probability models
    #real
    realData = getData(open(sys.argv[1],"r"))   
    realMod = learnPWM(bases, realData)
    #false
    falseData = getData(open(sys.argv[2],"r"))
    falseMod = learnPWM(bases, falseData)
    #test out test sequences
    testData = getData(open(sys.argv[3],"r"))
    outFile = open(sys.argv[4], "w+")
    runTests(testData, outFile, realMod, falseMod, bases)
    outFile.close()

if __name__=='__main__':
    main()

