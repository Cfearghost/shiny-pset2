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


int** add2(int n, int id, int** memory, int** X, int** Y){ 
    int** C = NULL;
    int sz = n*n +n;
    C = makeMatrix2(C, memory, n, id, sz);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
         C[i][j] = X[i][j] + Y[i][j];
    return C;
}

int** subtract2(int n, int id, int** memory, int** X, int** Y){ 
    int** C = NULL;
    int sz = n*n +n;
    C = makeMatrix2(C, memory, n, id, sz);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
         C[i][j] = X[i][j] - Y[i][j];
    return C;
}

// Strassen Algorithm 
int** strassenAlg(int n, int** X, int** Y) {
    if (n <= CROSS_OVER){
        int** C = matrixProductHeap2(n, X, Y);
        return C;
    }
    else if (n % 2 == 0){
        // dimesnion of submatrices
        int half_n = n/2;

        int matrix_size = half_n + half_n * half_n;

        int** allocated_memory = (int**) malloc(52*matrix_size * sizeof(int));

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
        int** p1 = strassenAlg(half_n, A, subtract2(half_n, 8, allocated_memory, F, H));          //p1 = A(F - H)
        int** p2 = strassenAlg(half_n, add2(half_n, 9, allocated_memory, A, B), H);                 //p2 = (A + B)H
        int** p3 = strassenAlg(half_n, add2(half_n, 10, allocated_memory, C, D), E);                 //p3 = (C + D)E
        int** p4 = strassenAlg(half_n, D, subtract2(half_n, 11, allocated_memory, G, E));         //p4 = D(G-E)
        int** p5 = strassenAlg(half_n, add2(half_n, 12, allocated_memory, A, D), add2(half_n, 13, allocated_memory, E, H));       //p5 = (A + D)(E + H)  
        int** p6 = strassenAlg(half_n, subtract2(half_n, 14, allocated_memory, B, D),  add2(half_n, 15, allocated_memory, G, H));   //p6 = (B - D)(G + H)
        int** p7 = strassenAlg(half_n, subtract2(half_n, 16, allocated_memory, C, A), add2(half_n, 17, allocated_memory, E, F));    //p7 = (A - C)(E + F)

        // calculating submatrices of C           
        int** AE_plus_BG = subtract2(half_n, 18, allocated_memory, add2(half_n, 19, allocated_memory, add2(half_n, 20, allocated_memory, p5, p4), p6), p2); 
        int** AF_plus_BH = add2(half_n, 21, allocated_memory, p1, p2); 
        int** CE_plus_DG = add2(half_n, 22, allocated_memory,p3, p4);  
        int** CF_plus_DH = subtract2(half_n, 23, allocated_memory, add2(half_n, 24, allocated_memory, add2(half_n, 25, allocated_memory, p5, p1), p7), p3);

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
        free(allocated_memory);
        //print("C -> %i", C1[0][0]);
        return C1;
    }
    
    else { 
        int matrix_size = (n+1) + (n+1) * (n+1);

        int** allocated_memory = (int**) malloc(4*matrix_size * sizeof(int));

        int** EvenX = makeMatrix2(EvenX, allocated_memory, (n+1), 0, matrix_size);
        int** EvenY = makeMatrix2(EvenX, allocated_memory, (n+1), 1, matrix_size);
        
        // DONT DELETE THIS, EVERYTHING FAILS WITHOUT IT
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

        int** C = makeMatrix(n);
          
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                C[i][j] = EvenC[i][j];
            }
        }

        free(allocated_memory);
        return C;
    }
}

int main(){
    int n = 2;
    int matrix_size = n + n*n;
    int** allocated_memory = (int**) malloc(4*matrix_size * sizeof(int));
    int** A = makeMatrix2(A, allocated_memory, n, 0, matrix_size);
    int** B = makeMatrix2(B, allocated_memory, n, 1, matrix_size);

    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++){
         A[i][j] = 1;
         B[i][j] = 1;
      }
    clock_t t = clock();
    int** C = strassenAlg(n, A, B);
    t = clock() - t; 
    // Calculate the time 
    float time = ((float)t)/CLOCKS_PER_SEC;
    printf("%f seconds \n", time); 

    free(allocated_memory);
    printMatrixHeap(n, C);
    return 0;
}

