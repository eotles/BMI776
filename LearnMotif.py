import sys
import copy
import math

#Get the probability score for the window (inside motif)
def windowScore(pee, sequenceNumber, winRange):
    windowScore = 1
    for j in winRange:
        char = sequenceData[sequenceNumber][j]
        charScore = pee.get(char)[j-winRange[0]]
        windowScore = windowScore*charScore
    return windowScore

#Get the probability score for the background
def backgroundScore(back, sequenceNumber, backRange):
    backgroundScore = 1
    for j in backRange:
        char = sequenceData[sequenceNumber][j]
        charScore = back.get(char)
        backgroundScore = backgroundScore*charScore
    return backgroundScore

#E-Step: move "zee-window" (zwin) about
def eStep(back, pee, zee):
    probScore = 0
    for zwinI,zRow in enumerate(zee):
        zRowScore = 0
        winIJScore = 0 
        for zwinJ,zRowCol in enumerate(zRow):
            #inside window range & score
            zwinRange = range(zwinJ,zwinJ+motifWidth)
            winScore = windowScore(pee, zwinI, zwinRange)
            #background range &score
            zbackRange = range(0,zwinJ) + range(zwinJ+motifWidth,len(sequenceData[zwinI]))
            backScore=backgroundScore(back, zwinI, zbackRange)
            #total pos score
            winIJScore = winScore*backScore
            zee[zwinI][zwinJ] = winIJScore
            zRowScore = zRowScore+winIJScore
        for zwinJ, zRowCol in enumerate(zRow):
            zee[zwinI][zwinJ] = zee[zwinI][zwinJ]/zRowScore
        probScore = math.log(zRowScore) + probScore
    return probScore

#M-Step        
def mStep(back, pee, zee):   
    for zwinI,zRow in enumerate(zee):
        for zwinJ,zRowCol in enumerate(zRow):          
            for motPos in range(0,motifWidth):
                char = sequenceData[zwinI][zwinJ+motPos]
                count = pee.get(char)
                count[motPos] = count[motPos]+zee[zwinI][zwinJ]
                pee.update({char: count})
        
    peeTotal = [0]*motifWidth
    for motPos in range(0,motifWidth):
        for char in bases:
            peeTotal[motPos] = peeTotal[motPos] + pee.get(char)[motPos]
        for char in bases:
            peeNewVal = pee.get(char)
            peeNewVal[motPos] = peeNewVal[motPos]/peeTotal[motPos]
            pee.update({char: peeNewVal})
    
    back = dict()
    backTotal = 0
    for char in bases:
        backProb = backCounts.get(char)-sum(pee.get(char))
        back.update({char: backProb})
        backTotal = backTotal + backProb
    for char in bases:
        back.update({char: back.get(char)/backTotal})   

#Return max likelihood starting positions for motif
def getMotifStart(zee):
    motifStart = [0]*len(sequenceData)
    for zwinI, zRow in enumerate(zee):
        maxFound = 0
        for zwinJ, zRowCol in enumerate(zRow):
            if(zRowCol>maxFound):
                maxFound = zRowCol
                motifStart[zwinI] = zwinJ
    return(motifStart)

#Return a string of representation of the p-matrix
def toStringMotifMatrix(pee):            
    outString = ""
    for base in pee:
        outString +=(str(base)+": "+str(pee.get(base)) +"\n")
    return outString

#Return the max prob motif (string of bases)
def toStringMotif(bases, pee):
    outString = ""
    for motPos,throwAway in enumerate(pee.get(bases[0])):
        maxPee = 0
        maxBase = ""
        peeBasePos = 0
        for base in pee:
            peeBasePos = pee.get(base)[motPos]
            if(peeBasePos > maxPee):
                maxPee = peeBasePos
                maxBase = base
        outString += maxBase
    return(outString)

#generate a p/theta matrix from progenitor motif sequence
def wesP(probMax,propMotif):
    newPee = dict()
    for motPos,throwAway in enumerate(propMotif):
        probPos = probMax
        for base in bases:
            pArr = newPee.get(base)
            if(motPos == 0): pArr = [0]*len(propMotif)
            if(base == propMotif[motPos]):
                pArr[motPos] = (probPos)
            else:
                pArr[motPos] = ((1.0-probPos)/(len(bases)-1.0))
            newPee.update({base: pArr})
    return newPee

#OutputFile Functions
def createPositionsFile(motPos):
    positions_file = open('positions_file','w+')
    for seqPos in motPos:
        positions_file.write(str(seqPos) + "\n")

def createModelFile(back, pee):
    model_file = open('model_file','w+')
    line = list()
    for char in bases:
        line.append(char)
        line.append(back.get(char))
        line = line + pee.get(char)
        for col in line:
            model_file.write(str(col) + "\t")
        model_file.write("\n")
        line = list()
    
        
        

#info - load sequence data and motif width eventually have the bases be filled out by
#analyzing input file instead of assuming DNA nucleotides
sequences_file = open(sys.argv[1],"r")
motifWidth = int(sys.argv[2])
bases = ['A','C','G','T']
#load sequences
sequenceData = [list(line.strip()) for line in sequences_file]
# 
#initialize pheta (pee/theta) matrix
pheta = dict()
phetaZero = dict()
for base in bases: 
    pheta.update({base: [1.0/len(bases)]*motifWidth})
    phetaZero.update({base: [1.0/len(bases)]*motifWidth})

#initialize background
backCounts = dict()
backZero = dict()
totalChar = 0
for seq in sequenceData:
    for char in seq:
        totalChar = totalChar + 1.0
        if(backCounts.has_key(char)):
            backCounts.update({char: backCounts.get(char)+1.0})
        else:
            backCounts.update({char: 1.0})            
for base in bases:
    backZero.update({base: backCounts.get(base)/totalChar})
    

#initialize z
zZero = [[(1.0/(len(seq)-motifWidth+1))]*(len(seq)-motifWidth+1) for seq in sequenceData]
maxZpos = [1]*len(sequenceData)


#initialize maxima tacking stuff
maxima = -999999
maxPee = dict()
maxBack = dict()
maxMotPos = list()
epsilon = .1

#Since we are running OOPS we can operate under the assumption that each sequence
#will have 1 instance of the motif.  If we go through every set set of width n
#sequential motif sequences we are guaranteed to use a variation on the motif
#once in the enumeration 
for i,c in enumerate(zZero[0]):
    print("RR "+`i`+"------------------------------------------------'")
    z = copy.deepcopy(zZero)
    back = copy.deepcopy(backZero)
    pheta = wesP(.7,sequenceData[0][i:i+motifWidth])
    
    stepNo = 0
    tmpScore = oldScore = -999999
    dif = 9999999
    maximaDif = 0
    #loop motif test until inside convergence range (epsilon)
    while (dif>epsilon):
        stepNo = stepNo+1
        tmpScore = eStep(back, pheta, z)
        dif = math.fabs(tmpScore - oldScore)
        oldScore = tmpScore
        mStep(back, pheta, z)
        motifPos = getMotifStart(z)
        print("     |Step "+`'{0:02d}'.format(stepNo)`+" - Score: " +`'{0:.2f}'.format(tmpScore)` + " - " + `toStringMotif(bases, pheta)`)
    #store information if best found
    if(tmpScore > maxima):
        maxima = tmpScore
        maxPee = copy.deepcopy(pheta)
        maxBack = copy.deepcopy(back)
        maxMotPos = copy.deepcopy(motifPos)
    print("     +------------------------------------------------\n     |Maxima: " + `'{0:.02f}'.format(maxima)`+ " - " + `toStringMotif(bases,maxPee)`)

        
print("\nP/theta:\n"+toStringMotifMatrix(maxPee))
print("\nStart Positions:\n"+str(maxMotPos))
createPositionsFile(maxMotPos)
createModelFile(maxBack, maxPee)
