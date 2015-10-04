'''
Created on Mar 8, 2014

@author: eotles
'''

from operator import concat
import copy

class mrkvTrans(object):
    '''
    classdocs
    '''


    def __init__(self, bases, maxOrder):
        '''
        Constructor
        '''
        self.bases = bases;
        self.maxOrder = maxOrder;
        self.allOrders = [dict() for _ in xrange(maxOrder+1)]
        self.buildOrderDicts(bases, 0)
        
    def genNextOrderDict(self, seqDict,bases):
        newSeqDict = dict()
        for seq in seqDict:
            for base in bases:
                newSeq = concat(seq, base)
                #newSeq = copy.deepcopy(seq) + base
                #print(newSeq)
                newSeqDict.update({newSeq:1})
        return newSeqDict

    def buildOrderDicts(self, bases, order):
        if(order>self.maxOrder):
            return
        orderDict = dict()
        if(order == 0):
            for base in bases:
                #newSeq = list()
                #newSeq.append(base)
                orderDict.update({base: 1})
        else:
            orderDict = self.genNextOrderDict(self.allOrders[order-1], bases)
        self.allOrders[order] = orderDict
        return self.buildOrderDicts(bases, order+1)
    
    def getOrderDict(self,order):
        return(self.allOrders[order])
    
    def update(self, subseq):
        currMrkTrans = self.allOrders[len(subseq)-1]
        #print("     TRANS: " + `currMrkTrans`)
        currVal = currMrkTrans.get(subseq)
        #print("     Count: "+`currVal`)
        #print(self.allOrders[len(subseq)-1].get(subseq))
        currMrkTrans.update({subseq: currVal+1})
        
    def normalize(self):
        for orderDict in self.allOrders:
            count = 0.0
            for key,value in orderDict.iteritems():
                count += value;
            for key,value in orderDict.iteritems():
                orderDict.update({key: float(value)/count})
            print(orderDict)
            
    def toString(self):
        output = ""
        for orders in self.allOrders:
            print(orders)
        #return output
            