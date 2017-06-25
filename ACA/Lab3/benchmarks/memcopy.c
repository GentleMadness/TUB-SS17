#include <stdio.h>
#include <stdlib.h>

int main(void)
{
	int i,j,k;	
	long long *array;
	long long *array2;
	long long *array3;
	long long *array4;
	
	array= (long long*)malloc (8*65536);
	array2= (long long*)malloc (8*65536);
	array3= (long long*)malloc (8*65536);
	array4= (long long*)malloc (8*65536);

	
	memcpy(array2, array, 8*65536);
	memcpy(array4, array3, 8*65536);
	
	for (i=0; i<8; i++){
		for (j=0; j<8; j++){			
			memcpy(array2+8192*j, array+8192*i,65536);
			memcpy(array4+8192*j, array3+8192*i,65536);
			for (k=0; k<8; k++){						
				array[111]= array3[222];				
				array[430]= array3[6543];
				array[1201]= array3[13050];
				array[5612]= array3[1470];
				array[1047]= array3[7806];
				array[9872]= array3[24676];
				array[24433]= array3[22213];
				array[6870]= array3[13];

				array4[111*k]= array2[222*j];				
				array4[430*k]= array2[6543*j];
				array4[1201*k]= array2[1305*j];
				array4[5612*k]= array2[1470*j];
				array4[1047*k]= array2[7806*j];
				array4[3872*k]= array2[2467*j];
				array4[2443*k]= array2[2221*j];
				array4[6870*k]= array2[13*j];

				array4[111*i]= array2[222*k];				
				array4[430*i]= array2[6543*k];
				array4[1201*i]= array2[1305*k];
				array4[5612*i]= array2[1470*k];
				array4[1047*i]= array2[7806*k];
				array4[3872*i]= array2[2467*k];
				array4[2443*i]= array2[2221*k];
				array4[6870*i]= array2[13*k];

			}
		}
		
	}



}

