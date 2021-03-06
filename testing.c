#include <stdio.h>
#include <time.h>
#include <stdlib.h>

#define CROSS_OVER 

int** makeMatrix2(int** matrix, int** alloc, int n, int id, int matrix_size){
    matrix = (int**) (alloc + id*matrix_size);
    for (int i=0; i<n; i++)
        matrix[i] = (int *) (matrix + (i+1)*n);
    return matrix;
}


void printMatrixHeap(int n, int** A){
    for(int i = 0; i < n; i++){
        printf("\n");
        for(int j = 0; j < n; j++){
            printf("%d\t",A[i][j]);
        }
    }
}

void sub(int n, int** a, int** b, int** c)
{
  int i, j;
  for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        c[i][j] = a[i][j] - b[i][j];
      }
  }
}

void add(int n, int** a, int** b, int** c)
{
  int i, j;

  for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        c[i][j] = a[i][j] + b[i][j];
      }
  }
}
void matrixProduct(int n, int** X, int** Y, int** C) { 
  int i, j, sum, k;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      for (sum = 0., k = 0; k < n; k++)
          sum += X[i][k] * Y[k][j];
      C[i][j] = sum;
    }
  }
}

// Strassen Algorithm 
void strassenAlg(int n, int** X, int** Y, int** C, int** D, int cs) {
    if (n <= cs){
        matrixProduct(n, X, Y, C);
    }
    else if (n % 2 == 0) {
        // dimesnion of submatrices
        int half_n = n/2;

        int matrix_size = half_n + half_n * half_n;

        int** allocated_memory = (int**) malloc(32*matrix_size * sizeof(int));
        // initialize sub-matrices of X

        int** a11 = NULL;
        a11 = makeMatrix2(a11, allocated_memory, half_n, 0, matrix_size);

        int** a12 = NULL;
        a12 = makeMatrix2(a12, allocated_memory, half_n, 1, matrix_size);

        int** a21 = NULL;
        a21 = makeMatrix2(a21, allocated_memory, half_n, 2, matrix_size);

        int** a22 = NULL;
        a22 = makeMatrix2(a22, allocated_memory, half_n, 3, matrix_size);

        int** b11 = NULL;
        b11 = makeMatrix2(b11, allocated_memory, half_n, 4, matrix_size);

        int** b12 = NULL;
        b12 = makeMatrix2(b12, allocated_memory, half_n, 5, matrix_size);

        int** b21 = NULL;
        b21 = makeMatrix2(b21, allocated_memory, half_n, 6, matrix_size);

        int** b22 = NULL;
        b22 = makeMatrix2(b22, allocated_memory, half_n, 7, matrix_size);

        int** c11 = NULL;
        c11 = makeMatrix2(c11, allocated_memory, half_n, 8, matrix_size);

        int** c12 = NULL;
        c12 = makeMatrix2(c12, allocated_memory, half_n, 9, matrix_size);

        int** c21 = NULL;
        c21 = makeMatrix2(c21, allocated_memory, half_n, 10, matrix_size);

        int** c22 = NULL;
        c22 = makeMatrix2(c22, allocated_memory, half_n, 11, matrix_size);

        int** d11 = NULL;
        d11 = makeMatrix2(d11, allocated_memory, half_n, 12, matrix_size);

        int** d12 = NULL;
        d12 = makeMatrix2(d12, allocated_memory, half_n, 13, matrix_size);

        int** d21 = NULL;
        d21 = makeMatrix2(d21, allocated_memory, half_n, 14, matrix_size);

        int** d22 = NULL;
        d22 = makeMatrix2(d22, allocated_memory, half_n, 15, matrix_size);
        
        // dividing the matrices in 4 sub-matrices:
        for (int i = 0; i < half_n; i++){
            for (int j = 0; j < half_n; j++){
                a11[i][j] = X[i][j];                     // top left
                a12[i][j] = X[i][j + half_n];            // top right
                a21[i][j] = X[i + half_n][j];             // bottom left
                a22[i][j] = X[i + half_n][j + half_n];    // bottom right
 
                b11[i][j] = Y[i][j];                      // top left
                b12[i][j] = Y[i][j + half_n];            // top right
                b21[i][j] = Y[i + half_n][j];             // bottom left
                b22[i][j] = Y[i + half_n][j + half_n];    // bottom right
            }
        }

        sub(half_n, a12, a22, d11);
        add(half_n, b21, b22, d12);
        strassenAlg(half_n, d11, d12, c11, d21, cs);
        sub(half_n, a21, a11, d11);
        add(half_n, b11, b12, d12);
        strassenAlg(half_n, d11, d12, c22, d21, cs);
        add(half_n, a11, a12, d11);
        strassenAlg(half_n, d11, b22, c12, d12, cs);
        sub(half_n, c11, c12, c11);
        sub(half_n, b21, b11, d11);
        strassenAlg(half_n, a22, d11, c21, d12, cs);
        add(half_n, c21, c11, c11);
        sub(half_n, b12, b22, d11);
        strassenAlg(half_n, a11, d11, d12, d21, cs);
        add(half_n, d12, c12, c12);
        add(half_n, d12, c22, c22);
        add(half_n, a21, a22, d11);
        strassenAlg(half_n, d11, b11, d12, d21, cs);
        add(half_n, d12, c21, c21);
        sub(half_n, c22, d12, c22);
        add(half_n, a11, a22, d11);
        add(half_n, b11, b22, d12);
        strassenAlg(half_n, d11, d12, d21, d22, cs);
        add(half_n, d21, c11, c11);
        add(half_n, d21, c22, c22);

        for (int i = 0; i < half_n; i++){
            for (int j = 0; j < half_n; j++){
                C[i][j] = c11[i][j];                     // top left
                C[i][j + half_n] = c12[i][j];            // top right
                C[i + half_n][j] = c21[i][j];             // bottom left
                C[i + half_n][j + half_n] = c22[i][j];    // bottom right
            }
        }

        free(allocated_memory);
    }
    else {
        int matrix_size = (n+1) + (n+1) * (n+1);

        int** allocated_memory = (int**) malloc(8*matrix_size * sizeof(int));

        int** EvenX = makeMatrix2(EvenX, allocated_memory, (n+1), 0, matrix_size);
        int** EvenY = makeMatrix2(EvenY, allocated_memory, (n+1), 1, matrix_size);
        int** EvenC = makeMatrix2(EvenC, allocated_memory, (n+1), 2, matrix_size);
        int** EvenD = makeMatrix2(EvenD, allocated_memory, (n+1), 3, matrix_size);
        
        // DONT DELETE THIS, EVERYTHING FAILS WITHOUT IT
        for (int i = 0; i < n+1; i++)
          for (int j = 0; j < n+1; j++){
             EvenX[i][j] = 0;
             EvenY[i][j] = 0;
             EvenC[i][j] = 0;
             EvenD[i][j] = 0;
        }
    
        for (int i = 0; i < n; i++){
         for (int j = 0; j < n; j++){
             EvenX[i][j] = X[i][j];
             EvenY[i][j] = Y[i][j];
             EvenC[i][j] = C[i][j];
             EvenD[i][j] = D[i][j];
         }
        }    
        
        strassenAlg(n+1, EvenX, EvenY, EvenC, EvenD, cs);
          
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                C[i][j] = EvenC[i][j];
            }
        }
        free(allocated_memory);
    }
}

int main(int argc, char *argv[]){
  int flag = atoi(argv[1]);
  int n = atoi(argv[2]);
  int cross = atoi(argv[3]);
  float array[cross];
  int count = 5;
  int count_duplicate = 5;
  int indexer = 0;
  while (cross <= 1024) {
    while(count > 0){
      if(flag == 0){
        int matrix_size = n + n*n;
        int** allocated_memory = (int**) malloc(8*matrix_size * sizeof(int));
        int** A = makeMatrix2(A, allocated_memory, n, 0, matrix_size);
        int** B = makeMatrix2(B, allocated_memory, n, 1, matrix_size);
        int** C = makeMatrix2(A, allocated_memory, n, 2, matrix_size);
        int** D = makeMatrix2(B, allocated_memory, n, 3, matrix_size);

        for (int i = 0; i < n; i++){
              for (int j = 0; j < n; j++){
                 A[i][j] = 1;
                 B[i][j] = 1;
                 C[i][j] = 0;
                 D[i][j] = 0;
              }
        }
        
        clock_t t = clock();
        // cs as cross-over
        strassenAlg(n, A, B, C, D, cross);
        t = clock() - t; 
        // Calculate the time 
        float time = ((float)t)/CLOCKS_PER_SEC;
        array[indexer] = time;
        //printf("%f \n", time); 
        free(allocated_memory);
        //printMatrixHeap(n, C);
      } 
      else {
        int matrix_size = n + n*n;
        int** allocated_memory = (int**) malloc(6*matrix_size * sizeof(int));
        int** A = makeMatrix2(A, allocated_memory, n, 0, matrix_size);
        int** B = makeMatrix2(B, allocated_memory, n, 1, matrix_size);
        int** C = makeMatrix2(A, allocated_memory, n, 2, matrix_size);

        for (int i = 0; i < n; i++)
          for (int j = 0; j < n; j++){
             A[i][j] = 1;
             B[i][j] = 1;
             C[i][j] = 0;
        }

            clock_t t = clock();
            matrixProduct(n, A, B, C);
            t = clock() - t; 
            // Calculate the time 
            float time = ((float)t)/CLOCKS_PER_SEC;
            array[count] = time;
            printf("%f \n", time); 
            free(allocated_memory);
            printMatrixHeap(n, C);

          }
          counter++
      }
      count--;
      indexer++; 
   }
   float sum = 0;
      for (int i = 0; i < count_duplicate; i++){
          sum = sum + array[i];
      }
    float average = sum/5.0;
    printf("%f\n", average);
    count = 5;
    indexer = 0;
    cross = cross * 2;
  }
  return 0;
}

