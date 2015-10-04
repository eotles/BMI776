'''
Created on Mar 5, 2014

@author: eotles
'''

from hw3_shared import *
import copy

class MddNode(object):


    #Constructor of MDD node
    #auto-magically creates children nodes if it has too many sequences
    def __init__(self, bases, availablePos, trainData):
        self.availablePos = availablePos
        self.trainData = trainData
        self.pwm = learnPWM(bases, trainData)
        self.cons  = getCons(self.pwm, bases)
        self.chiScores = chiSqN(self.pwm, self.cons, trainData, bases)
        self.availableChiScores = [self.chiScores[i] for i in self.availablePos]  
        self.maxChi = max(self.availableChiScores)
        self.splitPos = self.chiScores.index(self.maxChi)
        self.splitChar = self.cons[self.splitPos]
        self.hasChildren = False
        if(len(self.trainData)>400):
            self.hasChildren = True
            match = list()
            noMatch = list()
            #make a new list of available positions (need to remove the pos we just used)
            newAvailablePos = copy.deepcopy(availablePos)
            newAvailablePos.remove(self.splitPos)
            self.split(match,noMatch)
            self.childMatch = MddNode(bases,newAvailablePos,match)
            self.childNoMatch = MddNode(bases,newAvailablePos,noMatch)
    
    #creates sets of data from provided list of sequence
    #one set is of all of the sequences that meet the consensus
    #the other is all of the ones that do not
    def split(self, match, noMatch):
        for seq in self.trainData:
            if(seq[self.splitPos]==self.splitChar):
                match.append(seq)
            else:
                noMatch.append(seq)
    
    #this is a printing class - its dumb            
    def showStuff(self):
        print("Pos:" + `self.splitPos` + " | Char: "+ `self.splitChar` + "| TD: "+ `len(self.trainData)` +" | ChiScore: " + `self.chiScores`)
        if(len(self.trainData)>400):
            print("match:")
            #self.childMatch.showStuff()
            print("no match:")
            #self.childNoMatch.showStuff()
    
    #traverses the tree from a root node (self)
    #to find the correct PWM to use for a given sequence
    def getClassPWM(self, testSeq):
        #print("trying to find pwm...")
        if(not self.hasChildren):
            #print(self.pwm)
            return(self.pwm)
        else:
            if(testSeq[self.splitPos]==self.splitChar):
                return self.childMatch.getClassPWM(testSeq)
            else:
                return self.childNoMatch.getClassPWM(testSeq)
        