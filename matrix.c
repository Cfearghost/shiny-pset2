#include <stdio.h>
#include <time.h>
#include <stdlib.h>

#define CROSS_OVER 1

// Makes matrix on the heap
int** makeMatrix(int n){
    int **C = (int **)malloc(n * sizeof(int *));
    for (int i=0; i<n; i++)
         C[i] = (int *)malloc(n * sizeof(int));
    return C;
}

// standard matrix multiplication 
void matrixProductHeap(int n, int** A, int** B) {
   int sums[n];
   for (int i = 0; i < n; i++){
       for (int j = 0; j < n; j++){
           int sum = 0;
           for (int k = 0; k < n; k++){
               sum = sum + A[i][k]*B[k][j];
           }
           sums[j] = sum;
       }
       for (int x = 0; x < n; x++){
           A[i][x] = sums[x];
       }
   }
}

void matrixProductStack(int n, int A[n][n], int B[n][n]) {
   int sums[n];
   for (int i = 0; i < n; i++){
       for (int j = 0; j < n; j++){
           int sum = 0;
           for (int k = 0; k < n; k++){
               sum = sum + A[i][k]*B[k][j];
           }
           sums[j] = sum;
       }
       for (int x = 0; x < n; x++){
           A[i][x] = sums[x];
       }
   }
}

void printMatrixHeap(int n, int** A){
    for(int i = 0; i < n; i++){
        printf("\n");
        for(int j = 0; j < n; j++){
            printf("%d\t",A[i][j]);
        }
    }
}

void printMatrixStack(int n, int A[n][n]){
    for(int i = 0; i < n; i++){
        printf("\n");
        for(int j = 0; j < n; j++){
            printf("%d\t",A[i][j]);
        }
    }
}

int** add(int n, int** X, int** Y){ 
    int** C = makeMatrix(n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
         C[i][j] = X[i][j] + Y[i][j];
    
    return C;
}

int** subtract(int n, int** X, int** Y){ 
    int** C = makeMatrix(n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
         C[i][j] = X[i][j] - Y[i][j];
    return C;
}


// Strassen Algorithm 
int** strassenAlg(int n, int** X, int** Y) {
    if (n <= CROSS_OVER){
        matrixProductHeap(n, X, Y);
        return X;
    }
    else if (n % 2 == 0){
        // dimesnion of submatrices
        int half_n = n/2;

        // initialize sub-matrices of X
        int** A = makeMatrix(half_n);
        int** B = makeMatrix(half_n);
        int** C = makeMatrix(half_n);
        int** D = makeMatrix(half_n);
              
        // initialize sub-matrices of Y
        int** E = makeMatrix(half_n);
        int** F = makeMatrix(half_n);
        int** G = makeMatrix(half_n);
        int** H = makeMatrix(half_n);
        
        for (int i = 0; i < half_n; i++){
          for (int j = 0; j < half_n; j++){
             A[i][j] = 0;
             B[i][j] = 0;
             C[i][j] = 0;
             D[i][j] = 0; 
             E[i][j] = 0;
             F[i][j] = 0;
             G[i][j] = 0;
             H[i][j] = 0;
          }
        }
        /*
        A = [[0 for j in xrange(0, half_n)] for i in xrange(0, half_n)]
        B = [[0 for j in xrange(0, half_n)] for i in xrange(0, half_n)]
        C = [[0 for j in xrange(0, half_n)] for i in xrange(0, half_n)]
        D = [[0 for j in xrange(0, half_n)] for i in xrange(0, half_n)]
        */

        
        /*
        E = [[0 for j in xrange(0, half_n)] for i in xrange(0, half_n)]
        F = [[0 for j in xrange(0, half_n)] for i in xrange(0, half_n)]
        G = [[0 for j in xrange(0, half_n)] for i in xrange(0, half_n)]
        H = [[0 for j in xrange(0, half_n)] for i in xrange(0, half_n)]
        */
        // dividing the matrices in 4 sub-matrices:
        
        for (int i = 0; i < half_n; i++){
            for (int j = 0; j < half_n; j++){
                A[i][j] = X[i][j];                     // top left
                B[i][j] = X[i][j + half_n];            // top right
                C[i][j] = X[i + half_n][j];             // bottom left
                D[i][j] = X[i + half_n][j + half_n];    // bottom right
 
                E[i][j] = Y[i][j];                      // top left
                F[i][j] = Y[i][j + half_n];            // top right
                G[i][j] = Y[i + half_n][j];             // bottom left
                H[i][j] = Y[i + half_n][j + half_n];    // bottom right
            }
        }
        /*
        for i in xrange(0, half_n):
            for j in xrange(0, half_n):
                A[i][j] = X[i][j]                      // top left
                B[i][j] = X[i][j + half_n]             // top right
                C[i][j] = X[i + half_n][j]             // bottom left
                D[i][j] = X[i + half_n][j + half_n]    // bottom right
 
                E[i][j] = Y[i][j]                      // top left
                F[i][j] = Y[i][j + half_n]             // top right
                G[i][j] = Y[i + half_n][j]             // bottom left
                H[i][j] = Y[i + half_n][j + half_n]    // bottom right
        */
        
        // Problem with half n???
        
        int** p1 = makeMatrix(half_n);
        int** p2 = makeMatrix(half_n);
        int** p3 = makeMatrix(half_n);
        int** p4 = makeMatrix(half_n);
        int** p5 = makeMatrix(half_n);
        int** p6 = makeMatrix(half_n);
        int** p7 = makeMatrix(half_n);
        
        p1 = strassenAlg(half_n, A, subtract(half_n, F, H));          //p1 = A(F - H)
        p2 = strassenAlg(half_n, add(half_n, A, B), H);                 //p2 = (A + B)H
        p3 = strassenAlg(half_n, add(half_n, C, D), E);                 //p3 = (C + D)E
        p4 = strassenAlg(half_n, D, subtract(half_n, G, E));         //p4 = D(G-E)
        p5 = strassenAlg(half_n, add(half_n, A, D), add(half_n, E, H));       //p5 = (A + D)(E + H)  
        p6 = strassenAlg(half_n, subtract(half_n, B, D),  add(half_n, G, H));   //p6 = (B - D)(G + H)
        p7 = strassenAlg(half_n, subtract(half_n, C, A), add(half_n, E, F));    //p7 = (A - C)(E + F)
        
        
        /*
        // Calculating p1 to p7:
        p1 = strassenAlg(A, subtract(F, H))            //p1 = A(F - H)
        p2 = strassenAlg(add(A, B), H)                 //p2 = (A + B)H
        p3 = strassenAlg(add(C, D), E)                 //p3 = (C + D)E
        p4 = strassenAlg(D, subtract(G, E))            //p4 = D(G-E)
        p5 = strassenAlg(add(A, D), add(E, H))         //p5 = (A + D)(E + H)  
        p6 = strassenAlg(subtract(B, D),  add(G, H))   //p6 = (B - D)(G + H)
        p7 = strassenAlg(subtract(C, A), add(E, F))    //p7 = (A - C)(E + F)
        */
        
        // calculating submatrices of C
        
        int** AE_plus_BG = makeMatrix(half_n);
        int** AF_plus_BH = makeMatrix(half_n);
        int** CE_plus_DG = makeMatrix(half_n);
        int** CF_plus_DH = makeMatrix(half_n);
        
        AE_plus_BG = subtract(half_n, add(half_n, add(half_n, p5, p4), p6), p2); 
        AF_plus_BH = add(half_n, p1, p2); 
        CE_plus_DG = add(half_n, p3, p4);  
        CF_plus_DH = subtract(half_n, add(half_n, add(half_n, p5, p1), p7), p3);
 
        // Grouping the results obtained in a single matrix: 
        int** C1 = makeMatrix(n);
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
               C1[i][j] = 0;
            }
        }
        for (int i = 0; i < half_n; i++){
            for (int j = 0; j < half_n; j++){
                C1[i][j] = AE_plus_BG[i][j];
                C1[i][j + half_n] = AF_plus_BH[i][j];
                C1[i + half_n][j] = CE_plus_DG[i][j];
                C1[i + half_n][j + half_n] = CF_plus_DH[i][j];
            }
         }
       
        /*for i in xrange(0, half_n):
            for j in xrange(0, half_n):
                C[i][j] = AE_plus_BG[i][j]
                C[i][j + half_n] = AF_plus_BH[i][j]
                C[i + half_n][j] = CE_plus_DG[i][j]
                C[i + half_n][j + half_n] = CF_plus_DH[i][j] */
        return C1;
    }
    
    else {
        int** EvenX = makeMatrix(n+1);
        int** EvenY = makeMatrix(n+1);
        
        for (int i = 0; i < n+1; i++)
          for (int j = 0; j < n+1; j++){
             EvenX[i][j] = 0;
             EvenY[i][j] = 0;
          }
    
        for (int i = 0; i < n; i++){
         for (int j = 0; j < n; j++){
             EvenX[i][j] = X[i][j];
             EvenY[i][j] = Y[i][j];
         }
        }    
        
        // Help Ben???????????????
        int** EvenC = strassenAlg(n+1,EvenX, EvenY);
        
        int** C = makeMatrix(n);
         for (int i = 0; i < n; i++)
          for (int j = 0; j < n; j++){
             X[i][j] = 0;
             Y[i][j] = 0;
          }
          
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                C[i][j] = EvenC[i][j];
            }
        }
        return C;
    /*
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
    */
    }
}

int main(){
    int n = 4;
    if (n == 4){
        int** A = makeMatrix(n);
        int** B = makeMatrix(n);
        int** C = makeMatrix(n);

        for (int i = 0; i < n; i++)
          for (int j = 0; j < n; j++){
             A[i][j] = 1;
             B[i][j] = 1;
          }
        clock_t t = clock();
        C = strassenAlg(n, A, B);
        t = clock() - t; 
        // Calculate the time 
        float time = ((float)t)/CLOCKS_PER_SEC;
        printf("%f seconds \n", time); 
        printMatrixHeap(n, C);
    }/*
    else if (n < 1024){
        int A[n][n];
        int B[n][n];

        for (int i = 0; i < n; i++)
          for (int j = 0; j < n; j++){
             A[i][j] = 1;
             B[i][j] = 1;
          }
        
        clock_t t = clock();
        matrixProductStack(n, A, B);
        t = clock() - t; 
        // Calculate the time 
        float time = ((float)t)/CLOCKS_PER_SEC;
        printf("%f seconds \n", time); 
       // printMatrixStack(n, A);
    } */
    else {
        int** A = makeMatrix(n);
        int** B = makeMatrix(n);

        for (int i = 0; i < n; i++)
          for (int j = 0; j < n; j++){
             A[i][j] = 1;
             B[i][j] = 1;
          }
        
        clock_t t = clock();
        matrixProductHeap(n, A, B);
        t = clock() - t; 
        // Calculate the time 
        float time = ((float)t)/CLOCKS_PER_SEC;
        printf("%f seconds \n", time); 
        printMatrixHeap(n, A);
    }
    return 0;
}

