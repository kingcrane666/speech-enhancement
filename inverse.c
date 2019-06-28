#include "inverse.h"

void INV_3(double A[3][3], int n) // three-order matrix inversion  //blockinv1
{
	long i, j;
	double tp, tmp, temp[2];
	double a[2][2] ;

    temp[0]=A[0][2]/A[2][2];temp[1]=A[1][2]/A[2][2];
	for (i=0; i<2; i++)
	{
		for (j=0; j<2; j++)
		{
			a[i][j] = A[i][j]-temp[i]*A[2][j];
		}
	}
    //INV_2(a, 2);
	tp=a[0][1]/a[1][1];
	a[0][0]=1.0/(a[0][0]-tp*a[1][0]);
    a[0][1]=-a[0][0]*tp;
	a[1][1]=(1.0-a[1][0]*a[0][1])/a[1][1];
	a[1][0]=a[0][1];
	for (tmp=i=0; i<2; i++)
	{
		for (A[i][2]=j=0; j<2; j++)
		{
			A[i][2] -= a[i][j]*temp[j];//璐熷彿鍙?=涓?=
		}
        tmp += A[2][i]*A[i][2];	
		A[2][i] = A[i][2];
	}	
    A[2][2] = (1-tmp)/A[2][2];
    for (i=0; i<2; i++)
	{
		for (j=0; j<2; j++)
		{
			A[i][j] = a[i][j];
		}
	}
}
void Matrix_Mul_3(double a[3][3], double b[3][3], double c[3][3], int n) //three-order matrix multiplication,c=a*b
{
	int i, j, k;

	for (i=0; i<n; i++)
	{
		for (j=0; j<n; j++)
		{
			for (c[i][j]=k=0; k<n; k++) c[i][j] += a[i][k]*b[k][j];
		}
	}
}
void INV(double A[NUM1][NUM1], int n)
{
    int i, j;
	double a[3][3], b[3][3], c[3][3], d[3][3], temp1[3][3], temp2[3][3],a1[3][3], b1[3][3], d1[3][3];
	for (i=0; i<3; i++)
	{
		for (j=0; j<3; j++)
		{
			a[i][j] = A[i][j];
			b[i][j] = A[i][j+3];
			c[i][j] = A[i+3][j];
			d[i][j] = A[i+3][j+3];
		}
	}
	INV_3(d, 3);
	Matrix_Mul_3(b, d, temp1, 3);
    Matrix_Mul_3(temp1, c, a1, 3);
    for (i=0; i<3; i++)
	{
		for (j=0; j<3; j++)
		{
			a1[i][j] = a[i][j] - a1[i][j];
		}
	}
	INV_3(a1, 3);
	Matrix_Mul_3(a1, temp1, b1, 3);
	for (i=0; i<3; i++)
	{
		for (j=0; j<3; j++)
		{
			b1[i][j] = -b1[i][j];
		}
	}
	Matrix_Mul_3(c, b1, temp2,  3);		
	for (i=0; i<3; i++)
	{
		temp2[i][i] = 1 - temp2[i][i];
		for (j=0; j<3; j++)
		{
			if(j!=i) temp2[i][j] = -temp2[i][j];
		}
	}
	Matrix_Mul_3(d, temp2, d1, 3);
	for (i=0; i<3; i++)
	{
		for (j=0; j<3; j++)
		{
			A[i][j] = a1[i][j];
			A[i][j+3] = b1[i][j];
			A[i+3][j] = b1[j][i];
			A[i+3][j+3] = d1[i][j];
		}
	}
}

void Matrix_Mul(double a[NUM1][NUM1], double b[NUM1][NUM1], double c[NUM1][NUM1], int n) 
{
	int i, j, k;

	for (i=0; i<n; i++)
	{
		for (j=0; j<n; j++)
		{
			for (c[i][j]=k=0; k<n; k++) c[i][j] += a[i][k]*b[k][j];
		}
	}
}

void inverse(double A[NUM1][NUM1], double B[NUM1][NUM1], double C_re[NUM1][NUM1], double C_im[NUM1][NUM1], int n) 
{
	int i, j;
	double  temp_A[NUM1][NUM1], temp_AB[NUM1][NUM1], temp_BAB[NUM1][NUM1];
	
	for (i=0; i<n; i++)
	{
		for (j=0; j<n; j++) temp_A[i][j] = A[i][j];
	}
	INV(A, NUM1);
	Matrix_Mul(A, B, temp_AB, NUM1);
	Matrix_Mul(B, temp_AB, temp_BAB, NUM1);
	for (i = 0; i < NUM1; i++)
	{
		for (j = 0; j < NUM1; j++)	C_re[i][j] = temp_A[i][j] + temp_BAB[i][j];
	}
	INV(C_re, NUM1);
	Matrix_Mul(temp_AB, C_re, C_im, NUM1);
	for (i=0; i<n; i++)
	{
		for (j=0; j<n; j++) C_im[i][j] = -C_im[i][j];
	}
}

