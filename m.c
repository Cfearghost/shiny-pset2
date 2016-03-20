#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include "m.h"

#define CROSS_OVER 1

matrix newmatrix(int n)
{
    matrix a;

    a = (matrix)malloc(sizeof(*a));
    if (n <= CROSS_OVER) {
	    int i;

	    a->d = (double **)calloc(n, sizeof(double *));
        for (i = 0; i < n; i++) {
            a->d[i] = (double *)calloc(n, sizeof(double));
        }
    } 
    else {
	    n /= 2;
	    a->p = (matrix *)calloc(4, sizeof(matrix));
	    a11 = newmatrix(n);
	    a12 = newmatrix(n);
	    a21 = newmatrix(n);
	    a22 = newmatrix(n);
    }
    return a;
}

void randomfill(int n, matrix a){
    if (n <= CROSS_OVER) {
	    int i, j;
	    double **p = a->d;

	    for (i = 0; i < n; i++)
	        for (j = 0; j < n; j++)
		    p[i][j] = 1;
    } 
    else {
	    n /= 2;
	    randomfill(n, a11);
	    randomfill(n, a12);
	    randomfill(n, a21);
	    randomfill(n, a22);
    }
}


// Strassen Algorithm 
void multiply(int n, matrix a, matrix b, matrix c, matrix d){
    if (n <= CROSS_OVER) {
	    double sum, **p = a->d, **q = b->d, **r = c->d;
	    int i, j, k;

	    for (i = 0; i < n; i++) {
	        for (j = 0; j < n; j++) {
		    for (sum = 0., k = 0; k < n; k++)
		        sum += p[i][k] * q[k][j];
		    r[i][j] = sum;
	        }
	    }
    } 
    else {
	    n /= 2;
	    sub(n, a12, a22, d11);
	    add(n, b21, b22, d12);
	    multiply(n, d11, d12, c11, d21);
	    sub(n, a21, a11, d11);
	    add(n, b11, b12, d12);
	    multiply(n, d11, d12, c22, d21);
	    add(n, a11, a12, d11);
	    multiply(n, d11, b22, c12, d12);
	    sub(n, c11, c12, c11);
	    sub(n, b21, b11, d11);
	    multiply(n, a22, d11, c21, d12);
	    add(n, c21, c11, c11);
	    sub(n, b12, b22, d11);
	    multiply(n, a11, d11, d12, d21);
	    add(n, d12, c12, c12);
	    add(n, d12, c22, c22);
	    add(n, a21, a22, d11);
	    multiply(n, d11, b11, d12, d21);
	    add(n, d12, c21, c21);
	    sub(n, c22, d12, c22);
	    add(n, a11, a22, d11);
	    add(n, b11, b22, d12);
	    multiply(n, d11, d12, d21, d22);
	    add(n, d21, c11, c11);
	    add(n, d21, c22, c22);
    }
}

/* c = a+b */
void add(int n, matrix a, matrix b, matrix c){
    if (n <= CROSS_OVER) {
	    double **p = a->d, **q = b->d, **r = c->d;
	    int i, j;

	    for (i = 0; i < n; i++) {
	        for (j = 0; j < n; j++) {
		    r[i][j] = p[i][j] + q[i][j];
	        }
	    }
    } 
    else {
	    n /= 2;
	    add(n, a11, b11, c11);
	    add(n, a12, b12, c12);
	    add(n, a21, b21, c21);
	    add(n, a22, b22, c22);
    }
}

/* c = a-b */
void sub(int n, matrix a, matrix b, matrix c){
    if (n <= CROSS_OVER) {
	    double **p = a->d, **q = b->d, **r = c->d;
	    int i, j;

	    for (i = 0; i < n; i++) {
	        for (j = 0; j < n; j++) {
		    r[i][j] = p[i][j] - q[i][j];
	        }
	    }
    } 
    else {
	    n /= 2;
	    sub(n, a11, b11, c11);
	    sub(n, a12, b12, c12);
	    sub(n, a21, b21, c21);
	    sub(n, a22, b22, c22);
    }
}

void printma(int n, matrix a){
    if (n <= CROSS_OVER){
	    int i, j;
	    double **p = a->d;

	    for (i = 0; i < n; i++){
	        printf("\n");
            for (j = 0; j < n; j++)
		        printf("%f", p[i][j]);
        }
	     
	     printf("\n");
	 }
	 else{
	    n /= 2;
	    printma(n, a11);
	    printma(n, a12);
	    printma(n, a21);
	    printma(n, a22);
	 }
}

int main()
{
    int n = 8;
    matrix a, b, c, d;
    
    a = newmatrix(n);
    b = newmatrix(n);
    c = newmatrix(n);
    d = newmatrix(n);
    randomfill(n, a);
    randomfill(n, b);
    clock_t t = clock();
    multiply(n, a, b, c, d);	/* strassen algorithm */
    printma(n, c);
    t = clock() - t; 
    // Calculate the time 
    float time = ((float)t)/CLOCKS_PER_SEC;
    printf("%f seconds \n", time); 
    return 0;
}

