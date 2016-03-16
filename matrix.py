from optparse import OptionParser
from math import ceil, log

CROSS_OVER = 1

# standard matrix multiplication 
def matrixProduct(X, Y):
    n = len(X)
    C = [[0 for i in xrange(n)] for j in xrange(n)]
    for i in xrange(n):
        for k in xrange(n):
            for j in xrange(n):
                C[i][j] += X[i][k] * Y[k][j]
    return C

# Helper functions
def add(X, Y):
    C = [[X[i][j] + Y[i][j]  for j in range(len(X[0]))] for i in range(len(X))]
    return C

def subtract(X, Y):
    C = [[X[i][j] - Y[i][j]  for j in range(len(X[0]))] for i in range(len(X))]
    return C

# Strassen Algorithm 
def strassenAlg(X, Y):
    n = len(X)
    if n <= CROSS_OVER:
        return matrixProduct(X, Y)
    elif n % 2 == 0:
        # dimesnion of submatrices
        half_n = n/2

        # initialize sub-matrices of X
        A = [[0 for j in xrange(0, half_n)] for i in xrange(0, half_n)]
        B = [[0 for j in xrange(0, half_n)] for i in xrange(0, half_n)]
        C = [[0 for j in xrange(0, half_n)] for i in xrange(0, half_n)]
        D = [[0 for j in xrange(0, half_n)] for i in xrange(0, half_n)]
        
        # initialize sub-matrices of Y
        E = [[0 for j in xrange(0, half_n)] for i in xrange(0, half_n)]
        F = [[0 for j in xrange(0, half_n)] for i in xrange(0, half_n)]
        G = [[0 for j in xrange(0, half_n)] for i in xrange(0, half_n)]
        H = [[0 for j in xrange(0, half_n)] for i in xrange(0, half_n)]

        # dividing the matrices in 4 sub-matrices:
        for i in xrange(0, half_n):
            for j in xrange(0, half_n):
                A[i][j] = X[i][j]                      # top left
                B[i][j] = X[i][j + half_n]             # top right
                C[i][j] = X[i + half_n][j]             # bottom left
                D[i][j] = X[i + half_n][j + half_n]    # bottom right
 
                E[i][j] = Y[i][j]                      # top left
                F[i][j] = Y[i][j + half_n]             # top right
                G[i][j] = Y[i + half_n][j]             # bottom left
                H[i][j] = Y[i + half_n][j + half_n]    # bottom right

        # Calculating p1 to p7:
        p1 = strassenAlg(A, subtract(F, H))            # p1 = A(F - H)
        p2 = strassenAlg(add(A, B), H)                 # p2 = (A + B)H
        p3 = strassenAlg(add(C, D), E)                 # p3 = (C + D)E
        p4 = strassenAlg(D, subtract(G, E))            # p4 = D(G-E)
        p5 = strassenAlg(add(A, D), add(E, H))         # p5 = (A + D)(E + H)  
        p6 = strassenAlg(subtract(B, D),  add(G, H))   # p6 = (B - D)(G + H)
        p7 = strassenAlg(subtract(A, C), add(E, F))    # p7 = (A - C)(E + F)

        # calculating submatrices of C
        AE_plus_BG = subtract(add(add(p5, p4), p6), p2) 
        AF_plus_BH = add(p1, p2) 
        CE_plus_DG = add(p3, p4)  
        CF_plus_DH = subtract(add(add(p5, p1), p7), p3)
 
        # Grouping the results obtained in a single matrix: 
        C = [[0 for j in xrange(0, n)] for i in xrange(0, n)]
        for i in xrange(0, half_n):
            for j in xrange(0, half_n):
                C[i][j] = AE_plus_BG[i][j]
                C[i][j + half_n] = AF_plus_BH[i][j]
                C[i + half_n][j] = CE_plus_DG[i][j]
                C[i + half_n][j + half_n] = CF_plus_DH[i][j] 
        return C
    else:
        EvenX = [[0 for i in xrange(n+1)] for j in xrange(n+1)]
        EvenY = [[0 for i in xrange(n+1)] for j in xrange(n+1)]
        for i in xrange(n):
            for j in xrange(n):
                EvenX[i][j] = X[i][j]
                EvenY[i][j] = Y[i][j]
        EvenC = strassenAlg(EvenX, EvenY)
        C = [[0 for i in xrange(n)] for j in xrange(n)]
        for i in xrange(n):
            for j in xrange(n):
                C[i][j] = EvenC[i][j]
        return C


# Implementation that takes into account matrix dimension
# If we can get strassenAlg working above, we don't need this
def strassenImp(A, B):
    nextPowerOfTwo = lambda n: 2**int(ceil(log(n,2)))
    n = len(A)
    m = nextPowerOfTwo(n)
    APrep = [[0 for i in xrange(m)] for j in xrange(m)]
    BPrep = [[0 for i in xrange(m)] for j in xrange(m)]
    for i in xrange(n):
        for j in xrange(n):
            APrep[i][j] = A[i][j]
            BPrep[i][j] = B[i][j]
    CPrep = strassenAlg(APrep, BPrep)
    C = [[0 for i in xrange(n)] for j in xrange(n)]
    for i in xrange(n):
        for j in xrange(n):
            C[i][j] = CPrep[i][j]
    return C

# prints out a matrix
def printMatrix(matrix):
    for line in matrix:
        print "\t".join(map(str,line))

A = [[1,2,3,4],[1,2,5,6],[1,5,7,5],[1,5,7,5]]
B = [[1,2,3,4],[1,2,5,6],[1,5,7,5],[1,5,7,5]]
C = strassenAlg(A, B)
printMatrix(C)
C = strassenImp(A, B)
printMatrix(C)
C = matrixProduct(A, B)
printMatrix(C)
