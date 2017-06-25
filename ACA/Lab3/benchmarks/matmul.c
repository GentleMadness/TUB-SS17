#include <stdio.h>

#define SIZE 50


int dot(int a[SIZE][SIZE], int b[SIZE][SIZE], int i, int j)
{
	int dot=0;
	int x,y;
	for(x=0,y=0;x<SIZE;x++,y++)		
		dot+=a[i][x]*b[y][j];
	return dot;
}

int main()
{
    int a[SIZE][SIZE],b[SIZE][SIZE],c[SIZE][SIZE];
	int i,j;

	for(i=0;i<SIZE;i++)
		for(j=0;j<SIZE;j++)
			a[i][j]=b[i][j]=i;

    for(i=0;i<SIZE;i++)
		for(j=0;j<SIZE;j++)
            c[i][j]=dot(a,b,i,j);
		
}
