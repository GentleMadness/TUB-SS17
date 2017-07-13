#include <unistd.h>
#include <pthread.h>
#include <sched.h>
#include <sys/syscall.h>
#include <errno.h>
#include <stdint.h>
#include <stdio.h>

//define the number here with the core assigned to you
#define CPUID_SMT0 2
#define CPUID_SMT1 (CPUID_SMT0+6)

//enable and disable your SMT thread
#define SMT1_ON 1  //enable our SMT thread
#define COUNT_NUM 100000 //iteration times of tf_smt1 when tf_smt0 stall
static int fibonacci(int n){
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

static void *tf_smt0(void *targ){
  cpu_set_t mask;
  CPU_ZERO(&mask);
  CPU_SET(CPUID_SMT0, &mask);
  pid_t tid = syscall(SYS_gettid); //glibc does not provide a wrapper for gettid
  sched_setaffinity(tid, sizeof(cpu_set_t), &mask);

  printf("Fibonacci number 40 = %d\n", fibonacci(40));

  pthread_exit(NULL);
}
int while_count;
float array[COUNT_NUM];
static void *tf_smt1(void *targ){
  cpu_set_t mask;
  CPU_ZERO(&mask);
  CPU_SET(CPUID_SMT1, &mask);
  pid_t tid = syscall(SYS_gettid); //glibc does not provide a wrapper for gettid
  sched_setaffinity(tid, sizeof(cpu_set_t), &mask);
  while (1){
    // put your "background" code here
    // no sleeps are other code that schedules the thread out is allowed
    // inline assembly is not allowed

    /*************Our code and the brief explanation***********************
          
          COUNT_NUM = 100000, which is the number of times we would do this multiplication below. 
          We find that 100000 times of float point multiplication take 
          roughly the same time as time tf_smt0 is stalled.
          The reason to use float point numbers is that the original thread tf_smt0 deals with Integers,
          therefore when tf_smt0 is stalled, the the FP Multiplier is still available to be used.

          result: our thread tf_smt1 decreases the execution time 
          from around 740ms downto 
          an average of less than 580ms.
    ***********************************************************************/  
        while_count=0;
	      while(while_count<COUNT_NUM)   
		        array[while_count++]*=0.0;
  }

  //To be sure that the previous loop is not optimized away by the compiler
  printf("This should not appear\n");
  pthread_exit(NULL);
}

static int64_t diff(struct timespec start, struct timespec end)
{
  struct timespec temp;
  int64_t d;
  if ((end.tv_nsec-start.tv_nsec)<0) {
      temp.tv_sec = end.tv_sec-start.tv_sec-1;
      temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
  } else {
      temp.tv_sec = end.tv_sec-start.tv_sec;
      temp.tv_nsec = end.tv_nsec-start.tv_nsec;
  }
  d = temp.tv_sec*1000000000+temp.tv_nsec;
  return d;
}

int main ()
{
  pthread_t smt0, smt1;
  struct timespec start, end;

  clock_gettime(CLOCK_REALTIME,&start);

  pthread_create(&smt0, NULL, tf_smt0, NULL);
#if SMT1_ON
  pthread_create(&smt1, NULL, tf_smt1, NULL);
#endif

  pthread_join(smt0, NULL);
  clock_gettime(CLOCK_REALTIME,&end);
#if SMT1_ON
  pthread_cancel(smt1);
#endif

  double time = (double) diff(start,end)/1000000;
  printf("\nExecution time thread SMT0 \t %.3f ms\n", time);
  return 0;
}
