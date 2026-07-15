#include "matrix.h"
#include <math.h>
#include <stdio.h>

// #define M_PI (double) 22/7

int main(void) {
    // Multiply A by 2 vectors
    Matrix *A = matrix_new(2, 2);

    MAT_AT(A,0,0)=1;
    MAT_AT(A,0,1)=2;
    MAT_AT(A,1,0)=3;
    MAT_AT(A,1,1)=4;
    
    Matrix *v1 = matrix_new(2,1);
    MAT_AT(v1,0,0)=1;
    MAT_AT(v1,1,0)=0;
    
    Matrix *v2 = matrix_new(2,1);
    MAT_AT(v2,0,0)=0;
    MAT_AT(v2,1,0)=1;
    
    Matrix *r1 = matrix_mul(A, v1);
    Matrix *r2 = matrix_mul(A, v2);
    
    printf("A * [1 0]^T\n");
    matrix_print(r1);
    
    printf("\nA * [0 1]^T\n");
    matrix_print(r2);
    
    
   // Multiply 2 rotation matrixes 
    double c30 = cos(M_PI / 6.0);
    double s30 = sin(M_PI / 6.0);
    
    double c60 = cos(M_PI / 3.0);
    double s60 = sin(M_PI / 3.0);
    
    Matrix *R30 = matrix_new(2,2);
    Matrix *R60 = matrix_new(2,2);
    

    // Use the formula 
    // [ cos x, -sin x]
    // [ sin x, cos x]
    MAT_AT(R30,0,0)=c30;
    MAT_AT(R30,0,1)=-s30;
    MAT_AT(R30,1,0)=s30;
    MAT_AT(R30,1,1)=c30;
    
    MAT_AT(R60,0,0)=c60;
    MAT_AT(R60,0,1)=-s60;
    MAT_AT(R60,1,0)=s60;
    MAT_AT(R60,1,1)=c60;
    
    Matrix *R90 = matrix_mul(R30,R60);
   
    printf("\n");
    matrix_print(R90);



    // Multiply a 2x2 matrix with identity matrix
    Matrix *I = matrix_new(2,2);
    
    MAT_AT(I,0,0)=1;
    MAT_AT(I,0,1)=0;
    MAT_AT(I,1,0)=0;
    MAT_AT(I,1,1)=1;
    
    
    Matrix *T = matrix_new(2, 2);

    MAT_AT(T,0,0)=5;
    MAT_AT(T,0,1)=7;
    MAT_AT(T,1,0)=1;
    MAT_AT(T,1,1)=3;

    Matrix *C = matrix_mul(T,I);
    
    printf("\n");
    matrix_print(C);
    printf("\n");

    //Determinants
    printf("det(A) = %.2f\n", matrix_det(A)); //Reusing the A matrix, As the values are same

    Matrix *B = matrix_new(2,2);
    
    MAT_AT(B,0,0)=2;
    MAT_AT(B,0,1)=4;
    MAT_AT(B,1,0)=1;
    MAT_AT(B,1,1)=2;
    
    printf("det(B) = %.2f\n", matrix_det(B));


// Inverse and multiply it with inverse
    Matrix *Inv = matrix_inverse(A);
    
    printf("Inverse\n");
    matrix_print(Inv);
    
    Matrix *Check = matrix_mul(A, Inv);
    
    printf("\nA * A^-1\n");
    matrix_print(Check);
    printf("\n");



    // printf("%s\n", matrix_is_invertible(A) ? "Yes" : "No");


    // Matrix with Rank 1
    Matrix *R = matrix_new(2,2);

    MAT_AT(R,0,0)=1;
    MAT_AT(R,0,1)=2;
    MAT_AT(R,1,0)=2;
    MAT_AT(R,1,1)=4;
    
    matrix_print(R);
    printf("\ndet = %.2f\n", matrix_det(R));
}
