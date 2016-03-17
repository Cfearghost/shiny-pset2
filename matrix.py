import time

# Could put all answers of normal in on array, then another
# run through whole array

def printMatrix(matrix):
    for line in matrix:
        print "\t".join(map(str,line))

def matrixproduct(A, B):
    n = len(A)
    C = [[0 for i in xrange(n)] for j in xrange(n)]
    for i in xrange(n):
        for k in xrange(n):
            for j in xrange(n):
                C[i][j] += A[i][k] * B[k][j]
    return C



start_time = time.time()
A = [[1,1,1],[1,1,1],[1,1,1]]
B = [[1,1,1],[1,1,1],[1,1,1]]
C = matrixproduct(A, B)
printMatrix(C)
end = time.time()
print("%s seconds" % (time.time() - start_time))
start_time = time.time()
C = matrixproduct(A, B)
C = matrixproduct(A, B)
C = matrixproduct(A, B)
C = matrixproduct(A, B)
C = matrixproduct(A, B)
end = time.time()
print("%s seconds" % (time.time() - start_time))



