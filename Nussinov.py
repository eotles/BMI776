import sys
basePairs = {'A' : 'U', 'C' : 'G', 'G' : 'C', 'U' : 'A'}

#returns 1 if i j are base-pair-able
def delta(i,j):
    if(basePairs.get(i) == j):
        return 1
    else:
        return 0

#DP alg to build out gamma matrix    
def nussinov(seq):
    ELL = len(seq)
    gamma = [[0.0 for _ in xrange(ELL)] for _ in xrange(ELL)]
    for n in xrange(1,ELL):
        #get max # of paired bases
        for j in xrange(n,ELL):
            i = j-n
            down = gamma[i+1][j]
            left = gamma[i][j-1]
            diag = gamma[i+1][j-1] + delta(seq[i],seq[j])
            splitList = [0]
            #check splits
            for k in xrange(i+1,j):
                splitList.append(gamma[i][k]+gamma[k+1][j])
            split = max(splitList)
            gamma[i][j] = max(down, left, diag, split)
    return gamma

#traceback to find out where all the pairs are
def traceback(seq, i, j, gamma, pairs):
    if i<j:
        if gamma[i][j]==gamma[i+1][j]:
            traceback(seq, i+1, j, gamma, pairs)
        elif gamma[i][j]==gamma[i][j-1]:
            traceback(seq, i, j-1, gamma, pairs)
        #check if actual bair pairing if so add to pairs list
        elif gamma[i][j]==gamma[i+1][j-1]+delta(seq[i],seq[j]):
            pairs.append([i,j])
            traceback(seq, i+1, j-1, gamma, pairs)
    #check splits
    else:
        for k in xrange(i+1,j):
            if gamma[i][j]==gamma[i][k]+gamma[k+1][j]:
                traceback(seq, i, k, gamma, pairs)
                traceback(seq, k+1, j, gamma, pairs)
                break;
    return pairs

#display pair locations with format: i j
#indices i and j can be offset by a number
def disp(pairs,offSet):
    for _,pair in enumerate(pairs):
        print('%s %s' %(pair[0]+offSet,pair[1]+offSet))    

def main():
    seq = sys.argv[1]
    gamma = nussinov(seq)
    pairs = []
    traceback(seq, 0, len(seq)-1, gamma, pairs)
    disp(pairs,1)

if __name__=='__main__':
    main()