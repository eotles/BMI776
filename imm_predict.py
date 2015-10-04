'''
Created on Mar 6, 2014

@author: eotles
'''
from hw3_shared import *
from inhomMrkv import *


def runTests(testData, outFile, posIMM, negIMM, maxOrder, bases, chiVals):
    for seq in testData:
        posProb = 1
        negProb = 1
        for pos,_ in enumerate(seq):
                if(pos <= maxOrder):
                    start = 0
                else:
                    start = pos - maxOrder - 1

                subseq = ''.join(seq[start:pos+1])
                        
                posProb *= posIMM.immProb(pos,subseq,chiVals)
                negProb *= negIMM.immProb(pos,subseq,chiVals)
 
        score = math.log(posProb/negProb,2)
        #print(score)
        outFile.write(str(score)+"\n")
                
    


def main():
    #bases
    bases = ['A','C','G','T']
    maxOrder = 5
    #get MDD models
    #real
    realData = getData(open(sys.argv[1],"r"))
    realIMM = inhomMrkv(bases, maxOrder, realData)   
    #false
    falseData = getData(open(sys.argv[2],"r"))
    falseIMM = inhomMrkv(bases, maxOrder, falseData)
    
    chiValsFile = open("hw3_chisquare_df3_pvalues.txt")
    chiVals = [list(line.strip().split(" ")) for line in chiValsFile]
    chiValsFile.close()
    
    #test out test sequences
    testData = getData(open(sys.argv[3],"r"))
    outFile = open(sys.argv[4], "w+")
    runTests(testData, outFile, realIMM, falseIMM, maxOrder, bases, chiVals)    

if __name__=='__main__':
    main()