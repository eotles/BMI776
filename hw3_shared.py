'''
Created on Mar 4, 2014

@author: eotles
'''
import sys
import math

#open file, get data in list, close file
def getData(thisFile):
    return [list(line.strip()) for line in thisFile]
    thisFile.close()

#learn a model
def learnPWM(bases, data):
    #count with a 1 pseudo-count
    count = [[1.0 for _ in xrange(len(bases))] for _ in xrange(len(data[0]))]
    for seq in data:
        for i,char in enumerate(seq):
            count[i][bases.index(char)] += 1
    #normalize
    for i,posDist in enumerate(count):
        posSum = sum(posDist)
        for j,char in enumerate(posDist):
            count[i][j] = char/posSum
    return count

#calculate probability of a sequence based on a model
def calcProb(sequence, model, bases):
    prob = 1.0
    for i,char in enumerate(sequence):
        prob *= model[i][bases.index(char)]
    return prob

#score
def score(sequence, posModel, negModel, bases):
    return math.log(calcProb(sequence, posModel, bases)
                    /calcProb(sequence, negModel, bases),2);
                    
                    
#determine consensus sequence   
def getCons(pwm,bases):
    cons = [None for _ in xrange(len(pwm))]
    for i,charDist in enumerate(pwm):
        cons[i] = bases[charDist.index(max(charDist))]
    return cons

#calc pos i, j dependencies                  
def calcDependencies(pwm, cons, data, bases):
    dep = [[[[0 for _ in xrange(len(bases))] for _ in xrange(2)] for _ in xrange(len(data[0]))] for _ in xrange(len(data[0]))]
    for seq in data:
        for i,charI in enumerate(seq):
            for j,charJ in enumerate(seq):
                if(i!=j):
                    if(charI==cons[i]):
                        #consensus
                        dep[i][j][0][bases.index(charJ)] += 1
                    else:
                        #non-consensus
                        dep[i][j][1][bases.index(charJ)] += 1
    return(dep)

#calculate chiSq values for each position i - size n
#returns n vector of values
def chiSqN(pwm, cons, data, bases):
    chi = [0.0 for _ in xrange(len(data[0]))]
    nByN = chiSqNbyN(pwm, cons, data, bases)
    for n,nVals in enumerate(nByN):
        chi[n] = sum(nVals)
    return(chi)

#calculate chiSq values for each position pair (i,j) - size n^2
#returns n*n matrix of values
def chiSqNbyN(pwm, cons, data, bases):
    dep = calcDependencies(pwm, cons, data, bases)
    chi = [[0.0 for _ in xrange(len(data[0]))] for _ in xrange(len(data[0]))]
    for ni,chiI in enumerate(chi):
        for nj,chiJ in enumerate(chi):
            if(ni!=nj):
                chi[ni][nj] = chiSq(dep[ni][nj])
    return(chi)

# calculates the chiSq value for a consensus matrix
# return 1 value
def chiSq(consBaseMat):
    chi = 0.0    
    N = sum(consBaseMat[0]) + sum(consBaseMat[1])
    E_ij = 0.0
    for i, consState in enumerate(consBaseMat):
        R_i = sum(consState)
        for j,consStateBase in enumerate(consState):
            C_j = consBaseMat[0][j] + consBaseMat[1][j]
            E_ij = float(R_i*C_j)/N
            if(E_ij != 0):
                chi += ((consStateBase - E_ij)**2)/E_ij
    return(chi)



