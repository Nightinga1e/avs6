#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

clock_t start, stop;
long double *X;
long double *Y;
long double *R;

void dgemm_line(long double *X, long double *Y, long double *R, long int n){

	for( int i =0; i<n; i++)
		for(int j =0; j<n; j++){
			for(int k =0; k<n; k++)
				*(R+i * n+j) += *(X+i * n+k) * *(Y+k * n+j);
	}
}

void dgemm_block(long double *X, long double *Y, long double *R, long int n, long int BS){
    long int i, j, k, i0, j0, k0;
    long double *c0, *a0, *b0;

    for (i = 0; i < n; i += BS) {
        for (j = 0; j < n; j += BS) {
            for (k = 0; k < n; k += BS) {
                for (i0 = 0, c0 = (R + i * n + j), a0 = (X + i * n + k); i0 < BS; ++i0, c0 += n, a0 += n) {
                    for (k0 = 0, b0 = (Y + k * n + j); k0 < BS; ++k0, b0 += n) {
                        for (j0 = 0; j0 < BS; ++j0) {
                            c0[j0] += a0[k0] * b0[j0];
                        }
                    }
                }
            }
        }
    }
}

void init_matrix(long double *X, long double *Y, long double *R, long int n)
{
    long int i, j, k;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            *(X + i*n + j) = rand()%100 + 0.001 * (rand()%100);
            *(Y + i*n + j) = rand()%100 + 0.001 * (rand()%100);
            *(R + i*n + j) = 0.0;
        }
    }
}

void reset(long double *R, long int n){
	int i,j;
	for (i=0; i<n; i++){
		for(j=0; j<n; j++){
			*(R+i*n+j) = 0.0;
		}
	}
}

void print_matrix(long double *MAT, long int n)
{
    long int i, j;

    printf("Matrix:\n");
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("%12.2llf", *(MAT + i*n + j));
        }
        printf("\n");
    }
}


int main(){
	long int BS = 0;

	long int n;
	srand(time(0));
	printf("\n Enter block size: \n");
	scanf("%lld", &BS);
	printf("\n Enter matrix size: \n");
	scanf("%lld", &n);
	X = (long double*)malloc(n*n*sizeof(long double));
	Y = (long double*)malloc(n*n*sizeof(long double));
	R = (long double*)malloc(n*n*sizeof(long double));

	init_matrix(X,Y,R,n);
//	print_matrix(X,n);
//	print_matrix(Y,n);
//	print_matrix(R,n);
	start = clock();
	dgemm_line(X,Y,R,n);
	stop = clock();
	double timeresult_clock = (double)(stop-start)/CLOCKS_PER_SEC;
	printf("\n\n timeresult = %6.10lf \n", timeresult_clock);
//	print_matrix(R,n);
	reset(R,n);
	start = clock();
	dgemm_block(X,Y,R,n,BS);
	stop = clock();
	timeresult_clock = (double)(stop-start)/CLOCKS_PER_SEC;
	printf("\n\n timeresult = %6.10lf \n", timeresult_clock);
//	print_matrix(R,n);

	return 0;
}
