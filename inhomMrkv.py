'''
Created on Mar 8, 2014

@author: eotles
'''

import copy
#import mrkvTrans
from mrkvTrans import mrkvTrans
from operator import concat
import hw3_shared
from hw3_shared import chiSq


class inhomMrkv(object):


    def __init__(self, bases, maxOrder, data):
        self.bases = bases
        self.maxOrder = maxOrder
        #template = mrkvTrans(bases, maxOrder)
        self.template = mrkvTrans(bases, maxOrder)
        self.allMrkvTrans = [copy.deepcopy(self.template) for _ in xrange(len(data[0]))]
        self.train(data)
        order = 0
        #for mrkv in self.allMrkvTrans:
            #print(order)
            #mrkv.toString()
            #order +=1
        
    def train(self, data):
        for seq in data:
            for pos,_ in enumerate(seq):
            #for pos in xrange(0,len(seq)):
                #print("POS: "+`pos`)
                for length in xrange(0,self.maxOrder+1):
                    #start = pos-len-1
                    start = pos-length
                    end = pos +1
                    if(start>=0):
                        subseq = ''.join(seq[start:end])
                        #print(subseq)
                        self.allMrkvTrans[pos].update(subseq)
        #for mrkvTrans in self.allMrkvTrans:
        #    mrkvTrans.normalize()
        
    def immProb(self,pos,subseq,chiVals):
        order = len(subseq)-1
        posMrkvTrans = self.allMrkvTrans[pos]
        #len 1 return straight probability
        #print("pos: "+`pos` +"|subseq: "+ `subseq`)
        lambdas = self.getLambdas(posMrkvTrans, subseq, 40, chiVals)
        probs = self.getProbs(posMrkvTrans, subseq)
        #print(lambdas)
        #print(probs)
        probability = probs[0]
        for i,l in enumerate(lambdas):
            if(i>0):
                probability = l*probs[i]+(1-l)*probability  
        return probability

    
    def getProbs(self,posMrkvTrans,subseq):
        probabilities = [0.0 for _ in xrange(self.maxOrder+1)]
        order = self.maxOrder+1
        for prob in reversed(probabilities):
            order-=1
            startPos = len(subseq)-order-1
            if(startPos>=0):
                lookupSeq = subseq[startPos:]
                probabilities[order] = self.calcProb(lookupSeq, posMrkvTrans)
        return probabilities
        
    
    def getLambdas(self,posMrkvTrans,subseq, thresh, chiVals):
        lambdas = [0.0 for _ in xrange(self.maxOrder+1)]
        i = self.maxOrder+1
        for lamb in reversed(lambdas):
            i-=1
            order = i
            orderDict = posMrkvTrans.getOrderDict(i)
            startPos = len(subseq)-order-1
            if(startPos>=0):
                lookupSeq = subseq[startPos:]
                count = orderDict.get(lookupSeq)
                if(count>=thresh):
                    lambdas[i] = 1
                else:
                    x = self.calcD(lookupSeq,posMrkvTrans)
                    xPos = int(round(10*x)+1)
                    d = 0
                    if(xPos<=502):
                        d = float(chiVals[xPos][1])
                    if(d>=0.5):
                        lambdas[i] = d*float(count)/thresh
        return lambdas
        
    def calcD(self,lookupSeq,posMrkvTrans):
        more = self.getPossibilities(lookupSeq[1:])
        less = self.getPossibilities(lookupSeq[2:])
        count = [[] for _ in xrange(2)]
        count[0] = self.getCounts(more, posMrkvTrans)
        count[1] = self.getCounts(less, posMrkvTrans)
        return chiSq(count)
        
    def calcProb(self, lookupSeq, posMrkvTrans):
        possibilities = self.getPossibilities(lookupSeq[1:])
        counts = self.getCounts(possibilities, posMrkvTrans)
        orderDict = posMrkvTrans.getOrderDict(len(lookupSeq)-1)
        count = orderDict.get(lookupSeq)
        return float(count)/sum(counts)
    
    
    def getCounts(self,seqList,posMrkvTrans):
        orderDict = posMrkvTrans.getOrderDict(len(seqList[0])-1)
        count = list()
        for seq in seqList:
            #print(seq)
            #print(orderDict.get(seq))
            count.append(orderDict.get(seq))
        return count
    
    
    def getZeroOrderProb(self,posMrkvTrans,char):
            posDict = posMrkvTrans.getOrderDict(0)
            #possibilities = self.getPossibilities(subseq[1:])
            totalCount = 0
            for base in self.bases:
                totalCount += posDict.get(base)
            return float(posDict.get(char))/totalCount
    
    def getPossibilities(self,subseq):
        poss = list()
        for base in self.bases:
            poss.append(base+subseq)
        #print(poss)
        return poss

                    

        