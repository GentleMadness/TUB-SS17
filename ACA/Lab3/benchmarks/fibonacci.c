#include <stdio.h>
#define DEBUG

int fibonacci(int n){
	if(n==0){
		return 0;
	}
	else{
		if(n==1)
			return 1;		
		else
			return (fibonacci(n-1)+fibonacci(n-2));		
	}
}


int main(){
	int i; 
	for (i=0; i<20; i++)
	  fibonacci(i);
    	return 0;
}

