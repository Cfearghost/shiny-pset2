#include <stdio.h>
#include <time.h>
#include <stdlib.h>

#define CROSS_OVER 1

// Makes matrix on the heap
int** makeMatrix(int n){
    int **C = (int **)malloc(n * sizeof(int *));
    for (int i=0; i<n; i++)
    // Don't we need to malloc every column/row??
         C[i] = (int *)malloc(n * sizeof(int));
    return C;
}

int** makeMatrix2(int** matrix, int** alloc, int n, int id, int matrix_size){
    matrix = (int**) (alloc + id*matrix_size);
    for (int i=0; i<n; i++)
        matrix[i] = (int *) (matrix + (i+1)*n);
    return matrix;
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

int** matrixProductHeap2(int n, int** A, int** B) {
   int** C = makeMatrix(n);
   int sum = 0; 
   for (int i = 0; i < n; i++){
       for (int j = 0; j < n; j++){
           sum = 0;
           for (int k = 0; k < n; k++){
               sum = sum + A[i][k]*B[k][j];
           }
           C[i][j] = sum;
       }
   }
   return C;
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
    int** C = NULL;
    int** allocated_memory = (int**) malloc(1*n*n * sizeof(int));
    C = makeMatrix2(C, allocated_memory, n, 0, n*n);;
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
         C[i][j] = X[i][j] + Y[i][j];
    return C;
}

int** subtract(int n, int** X, int** Y){ 
    int** C = NULL;
    int** allocated_memory = (int**) malloc(1*n*n * sizeof(int));
    C = makeMatrix2(C, allocated_memory, n, 0, n*n);;
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
         C[i][j] = X[i][j] - Y[i][j];
    return C;
}


// Strassen Algorithm 
int** strassenAlg(int n, int** X, int** Y) {
    if (n <= CROSS_OVER){
        int** C = makeMatrix(n);
        C = matrixProductHeap2(n, X, Y);
        return C;
    }
    else if (n % 2 == 0){
        // dimesnion of submatrices
        int half_n = n/2;

        int matrix_size = half_n + half_n * half_n;

        int** allocated_memory = (int**) malloc(20*matrix_size * sizeof(int));
        // initialize sub-matrices of X

        int** A = NULL;
        A = makeMatrix2(A, allocated_memory, half_n, 0, matrix_size);

        int** B = NULL;
        B = makeMatrix2(B, allocated_memory, half_n, 1, matrix_size);

        int** C = NULL;
        C = makeMatrix2(C, allocated_memory, half_n, 2, matrix_size);

        int** D = NULL;
        D = makeMatrix2(D, allocated_memory, half_n, 3, matrix_size);

        int** E = NULL;
        E = makeMatrix2(E, allocated_memory, half_n, 4, matrix_size);

        int** F = NULL;
        F = makeMatrix2(F, allocated_memory, half_n, 5, matrix_size);

        int** G = NULL;
        G = makeMatrix2(G, allocated_memory, half_n, 6, matrix_size);

        int** H = NULL;
        H = makeMatrix2(H, allocated_memory, half_n, 7, matrix_size);
      
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

        // Calculating p1 to p7:
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
        
        int** EvenC = strassenAlg(n+1,EvenX, EvenY);
        
        free(EvenX);
        free(EvenY);

        int** C = makeMatrix(n);
         for (int i = 0; i < n; i++)
          for (int j = 0; j < n; j++){
             C[i][j] = 0;
        }
          
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                C[i][j] = EvenC[i][j];
            }
        }
        free(EvenC);

        return C;
    }
}

int main(){
    int n = 5;
    if (n == 5){
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
        free(A);
        free(B);
        free(C);
        printMatrixHeap(n, C);
    }
    else if (n == 3){
        int** A = makeMatrix(n);
        int** B = makeMatrix(n);
        int** C = makeMatrix(n);
        for (int i = 0; i < n; i++)
          for (int j = 0; j < n; j++){
             A[i][j] = 1;
             B[i][j] = 1;
             C[i][j] = 0;
          }
        
        C = matrixProductHeap2(n, A, B);
        // Calculate the time 
        printMatrixHeap(n, C);
    }
    /*
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
        //printMatrixHeap(n, A);
    }
    return 0;
}

