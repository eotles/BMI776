'''
Created on Mar 4, 2014

@author: eotles
'''

from MddNode import *

def runTests(testData, outFile, posMdd, negMdd, bases):
    for seq in testData:
        outFile.write(str(score(seq, posMdd.getClassPWM(seq), negMdd.getClassPWM(seq), bases))+"\n")

def main():     
    #bases
    bases = ['A','C','G','T']
    #get MDD models
    #real
    realData = getData(open(sys.argv[1],"r"))
    realMdd = MddNode(bases, range(len(realData[0])), realData)  
    #false
    falseData = getData(open(sys.argv[2],"r"))
    falseMdd = MddNode(bases, range(len(falseData[0])), falseData)
    #test out test sequences
    testData = getData(open(sys.argv[3],"r"))
    outFile = open(sys.argv[4], "w+")
    runTests(testData, outFile, realMdd, falseMdd, bases)

if __name__=='__main__':
    main()

