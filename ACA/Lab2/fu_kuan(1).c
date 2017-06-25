/*
Yuting Fu 387512, Shan Kuan 387516
*/


#include <immintrin.h>
#include <malloc.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

#define IDCT_SIZE         16
#define ITERATIONS        1000000
#define MAX_NEG_CROP      1024

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))

static const short g_aiT16[16][16] =
{
  { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64},
  { 90, 87, 80, 70, 57, 43, 25,  9, -9,-25,-43,-57,-70,-80,-87,-90},
  { 89, 75, 50, 18,-18,-50,-75,-89,-89,-75,-50,-18, 18, 50, 75, 89},
  { 87, 57,  9,-43,-80,-90,-70,-25, 25, 70, 90, 80, 43, -9,-57,-87},
  { 83, 36,-36,-83,-83,-36, 36, 83, 83, 36,-36,-83,-83,-36, 36, 83},
  { 80,  9,-70,-87,-25, 57, 90, 43,-43,-90,-57, 25, 87, 70, -9,-80},
  { 75,-18,-89,-50, 50, 89, 18,-75,-75, 18, 89, 50,-50,-89,-18, 75},
  { 70,-43,-87,  9, 90, 25,-80,-57, 57, 80,-25,-90, -9, 87, 43,-70},
  { 64,-64,-64, 64, 64,-64,-64, 64, 64,-64,-64, 64, 64,-64,-64, 64},
  { 57,-80,-25, 90, -9,-87, 43, 70,-70,-43, 87,  9,-90, 25, 80,-57},
  { 50,-89, 18, 75,-75,-18, 89,-50,-50, 89,-18,-75, 75, 18,-89, 50},
  { 43,-90, 57, 25,-87, 70,  9,-80, 80, -9,-70, 87,-25,-57, 90,-43},
  { 36,-83, 83,-36,-36, 83,-83, 36, 36,-83, 83,-36,-36, 83,-83, 36},
  { 25,-70, 90,-80, 43,  9,-57, 87,-87, 57, -9,-43, 80,-90, 70,-25},
  { 18,-50, 75,-89, 89,-75, 50,-18,-18, 50,-75, 89,-89, 75,-50, 18},
  {  9,-25, 43,-57, 70,-80, 87,-90, 90,-87, 80,-70, 57,-43, 25, -9}
};

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

static void compare_results(short *ref, short *res, const char *msg)
{
    int correct =1;

    printf("Comparing %s\n",msg);
    for(int j=0; j<IDCT_SIZE; j++)  {
        for(int i=0; i<IDCT_SIZE; i++){
            if(ref[j*IDCT_SIZE+i] != res[j*IDCT_SIZE+i]){
                correct=0;
                printf("failed at %d,%d\t ref=%d, res=%d\n ", i, j, ref[j*IDCT_SIZE+i],res[j*IDCT_SIZE+i]);
            }
        }
    }
    if (correct){
        printf("correct\n\n");
    }
}

// this function is for timing, do not change anything here
static void benchmark( void (*idct16)(short *, short *), short *input, short *output, const char *version )
{
    struct timespec start, end;
    clock_gettime(CLOCK_REALTIME,&start);

    for(int i=0;i<ITERATIONS;i++)
        idct16(input, output);

    clock_gettime(CLOCK_REALTIME,&end);
    double avg = (double) diff(start,end)/ITERATIONS;
    printf("%10s:\t %.3f ns\n", version, avg);
}

//scalar code for the inverse transform
static void partialButterflyInverse16(short *src, short *dst, int shift)
{
  int E[8],O[8];
  int EE[4],EO[4];
  int EEE[2],EEO[2];
  int add = 1<<(shift-1);

  for (int j=0; j<16; j++)
  {
    /* Utilizing symmetry properties to the maximum to minimize the number of multiplications */
    for (int k=0; k<8; k++)
    {
      O[k] = g_aiT16[ 1][k]*src[ 16] + g_aiT16[ 3][k]*src[ 3*16] + g_aiT16[ 5][k]*src[ 5*16] + g_aiT16[ 7][k]*src[ 7*16] +
        g_aiT16[ 9][k]*src[ 9*16] + g_aiT16[11][k]*src[11*16] + g_aiT16[13][k]*src[13*16] + g_aiT16[15][k]*src[15*16];
    }
    for (int k=0; k<4; k++)
    {
      EO[k] = g_aiT16[ 2][k]*src[ 2*16] + g_aiT16[ 6][k]*src[ 6*16] + g_aiT16[10][k]*src[10*16] + g_aiT16[14][k]*src[14*16];
    }
    EEO[0] = g_aiT16[4][0]*src[ 4*16 ] + g_aiT16[12][0]*src[ 12*16 ];
    EEE[0] = g_aiT16[0][0]*src[ 0    ] + g_aiT16[ 8][0]*src[  8*16 ];
    EEO[1] = g_aiT16[4][1]*src[ 4*16 ] + g_aiT16[12][1]*src[ 12*16 ];
    EEE[1] = g_aiT16[0][1]*src[ 0    ] + g_aiT16[ 8][1]*src[  8*16 ];

    /* Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector */
    for (int k=0; k<2; k++)
    {
      EE[k] = EEE[k] + EEO[k];
      EE[k+2] = EEE[1-k] - EEO[1-k];
    }
    for (int k=0; k<4; k++)
    {
      E[k] = EE[k] + EO[k];
      E[k+4] = EE[3-k] - EO[3-k];
    }
    for (int k=0; k<8; k++)
    {
      dst[k]   = MAX( -32768, MIN( 32767, (E[k]   + O[k]   + add)>>shift ));
      dst[k+8] = MAX( -32768, MIN( 32767, (E[7-k] - O[7-k] + add)>>shift ));
    }
    src ++;
    dst += 16;
  }
}

static void idct16_scalar(short* pCoeff, short* pDst)
{
  short tmp[ 16*16] __attribute__((aligned(16)));
  partialButterflyInverse16(pCoeff, tmp, 7);
  partialButterflyInverse16(tmp, pDst, 12);
}

/// CURRENTLY SAME CODE AS SCALAR !!
/// REPLACE HERE WITH SSE intrinsics
//scalar code for the inverse transform
static void partialButterflyInverse16_simd(short *src, short *dst, int shift)
{
  // SIMD Declaration of variables
  int add = 1<<(shift-1); //no SIMD replacement 
  
  // used to store new matrix:
  short gaiT[8][8] __attribute__ ((aligned (16)));
  short source[16][16] __attribute__ ((aligned (16)));
  __m128i *gaiT_vec = (__m128i *)gaiT;
  __m128i *source_vec = (__m128i *)source;
  //load g_aiT16, short type in c is 16 bit, each vector is 128 bit
  //cast the variables to SIMD vector type
  const  __m128i *in_vec = (__m128i *) src;       
  __m128i *out_vec = (__m128i *) dst;
  const __m128i *g_aiT16_vec = (__m128i *)g_aiT16;

  //  Load & Transpose
  //each g_aiT16_vec contains 8 short int
  //each row of g_aiT16 needs 2 vector for 16 shorts 
  //e.g  vec[0], vec[1] are the row 0, 
  //     vec[2], vec[3] the row 1

  /*====================================
  load and transpose the LHS of 
  odd row of g_aiT matrix
  -->to be used in the first loop in 

  each row of LHS 8 elements needs one 128-bit reg
  for (int j=0; j<16; j++) loop
  -->g_aiT16[1,3,5,7,9,11,13,15][k], k=0..7
  =======================================*/

  __m128i temp0 = _mm_unpacklo_epi64(_mm_unpacklo_epi32(_mm_unpacklo_epi16(_mm_load_si128(&g_aiT16_vec[2]), _mm_load_si128(&g_aiT16_vec[6])), _mm_unpacklo_epi16(_mm_load_si128(&g_aiT16_vec[10]), _mm_load_si128(&g_aiT16_vec[14]))),
                _mm_unpacklo_epi32(_mm_unpacklo_epi16(_mm_load_si128(&g_aiT16_vec[18]), _mm_load_si128(&g_aiT16_vec[22])), _mm_unpacklo_epi16(_mm_load_si128(&g_aiT16_vec[26]), _mm_load_si128(&g_aiT16_vec[30]))));  
  __m128i temp1 = _mm_unpackhi_epi64(_mm_unpacklo_epi32(_mm_unpacklo_epi16(_mm_load_si128(&g_aiT16_vec[2]), _mm_load_si128(&g_aiT16_vec[6])), _mm_unpacklo_epi16(_mm_load_si128(&g_aiT16_vec[10]), _mm_load_si128(&g_aiT16_vec[14]))),
                _mm_unpacklo_epi32(_mm_unpacklo_epi16(_mm_load_si128(&g_aiT16_vec[18]), _mm_load_si128(&g_aiT16_vec[22])), _mm_unpacklo_epi16(_mm_load_si128(&g_aiT16_vec[26]), _mm_load_si128(&g_aiT16_vec[30]))));
  __m128i temp2 = _mm_unpacklo_epi64(_mm_unpackhi_epi32(_mm_unpacklo_epi16(_mm_load_si128(&g_aiT16_vec[2]), _mm_load_si128(&g_aiT16_vec[6])), _mm_unpacklo_epi16(_mm_load_si128(&g_aiT16_vec[10]), _mm_load_si128(&g_aiT16_vec[14]))),
                _mm_unpackhi_epi32(_mm_unpacklo_epi16(_mm_load_si128(&g_aiT16_vec[18]), _mm_load_si128(&g_aiT16_vec[22])), _mm_unpacklo_epi16(_mm_load_si128(&g_aiT16_vec[26]), _mm_load_si128(&g_aiT16_vec[30]))));
  __m128i temp3 = _mm_unpackhi_epi64(_mm_unpackhi_epi32(_mm_unpacklo_epi16(_mm_load_si128(&g_aiT16_vec[2]), _mm_load_si128(&g_aiT16_vec[6])), _mm_unpacklo_epi16(_mm_load_si128(&g_aiT16_vec[10]), _mm_load_si128(&g_aiT16_vec[14]))),
                _mm_unpackhi_epi32(_mm_unpacklo_epi16(_mm_load_si128(&g_aiT16_vec[18]), _mm_load_si128(&g_aiT16_vec[22])), _mm_unpacklo_epi16(_mm_load_si128(&g_aiT16_vec[26]), _mm_load_si128(&g_aiT16_vec[30]))));
  __m128i temp4 = _mm_unpacklo_epi64(_mm_unpacklo_epi32(_mm_unpackhi_epi16(_mm_load_si128(&g_aiT16_vec[2]), _mm_load_si128(&g_aiT16_vec[6])), _mm_unpackhi_epi16(_mm_load_si128(&g_aiT16_vec[10]), _mm_load_si128(&g_aiT16_vec[14]))),
                _mm_unpacklo_epi32(_mm_unpackhi_epi16(_mm_load_si128(&g_aiT16_vec[18]), _mm_load_si128(&g_aiT16_vec[22])), _mm_unpackhi_epi16(_mm_load_si128(&g_aiT16_vec[26]), _mm_load_si128(&g_aiT16_vec[30]))));
  __m128i temp5 = _mm_unpackhi_epi64(_mm_unpacklo_epi32(_mm_unpackhi_epi16(_mm_load_si128(&g_aiT16_vec[2]), _mm_load_si128(&g_aiT16_vec[6])), _mm_unpackhi_epi16(_mm_load_si128(&g_aiT16_vec[10]), _mm_load_si128(&g_aiT16_vec[14]))),
                _mm_unpacklo_epi32(_mm_unpackhi_epi16(_mm_load_si128(&g_aiT16_vec[18]), _mm_load_si128(&g_aiT16_vec[22])), _mm_unpackhi_epi16(_mm_load_si128(&g_aiT16_vec[26]), _mm_load_si128(&g_aiT16_vec[30]))));
  __m128i temp6 = _mm_unpacklo_epi64(_mm_unpackhi_epi32(_mm_unpackhi_epi16(_mm_load_si128(&g_aiT16_vec[2]), _mm_load_si128(&g_aiT16_vec[6])), _mm_unpackhi_epi16(_mm_load_si128(&g_aiT16_vec[10]), _mm_load_si128(&g_aiT16_vec[14]))),
                _mm_unpackhi_epi32(_mm_unpackhi_epi16(_mm_load_si128(&g_aiT16_vec[18]), _mm_load_si128(&g_aiT16_vec[22])), _mm_unpackhi_epi16(_mm_load_si128(&g_aiT16_vec[26]), _mm_load_si128(&g_aiT16_vec[30]))));
  __m128i temp7 = _mm_unpackhi_epi64(_mm_unpackhi_epi32(_mm_unpackhi_epi16(_mm_load_si128(&g_aiT16_vec[2]), _mm_load_si128(&g_aiT16_vec[6])), _mm_unpackhi_epi16(_mm_load_si128(&g_aiT16_vec[10]), _mm_load_si128(&g_aiT16_vec[14]))),
                _mm_unpackhi_epi32(_mm_unpackhi_epi16(_mm_load_si128(&g_aiT16_vec[18]), _mm_load_si128(&g_aiT16_vec[22])), _mm_unpackhi_epi16(_mm_load_si128(&g_aiT16_vec[26]), _mm_load_si128(&g_aiT16_vec[30]))));

  //store the transposed matrix in gaiT_vec
  _mm_store_si128(&gaiT_vec[0], temp0);
  _mm_store_si128(&gaiT_vec[1], temp1); 
  _mm_store_si128(&gaiT_vec[2], temp2);
  _mm_store_si128(&gaiT_vec[3], temp3);
  _mm_store_si128(&gaiT_vec[4], temp4);
  _mm_store_si128(&gaiT_vec[5], temp5);
  _mm_store_si128(&gaiT_vec[6], temp6);
  _mm_store_si128(&gaiT_vec[7], temp7);
  /*=== End of loading g_aiT odd row LHS 8*8 ===*/

  //===Start loading g_aiT even row LHS 4 ==========

  temp0 = _mm_unpacklo_epi64(_mm_unpacklo_epi32(_mm_unpacklo_epi16(_mm_load_si128(&g_aiT16_vec[0]), _mm_load_si128(&g_aiT16_vec[4])), _mm_unpacklo_epi16(_mm_load_si128(&g_aiT16_vec[8]), _mm_load_si128(&g_aiT16_vec[12]))),
                _mm_unpacklo_epi32(_mm_unpacklo_epi16(_mm_load_si128(&g_aiT16_vec[16]), _mm_load_si128(&g_aiT16_vec[20])), _mm_unpacklo_epi16(_mm_load_si128(&g_aiT16_vec[24]), _mm_load_si128(&g_aiT16_vec[28]))));  
  temp1 = _mm_unpackhi_epi64(_mm_unpacklo_epi32(_mm_unpacklo_epi16(_mm_load_si128(&g_aiT16_vec[0]), _mm_load_si128(&g_aiT16_vec[4])), _mm_unpacklo_epi16(_mm_load_si128(&g_aiT16_vec[8]), _mm_load_si128(&g_aiT16_vec[12]))),
                _mm_unpacklo_epi32(_mm_unpacklo_epi16(_mm_load_si128(&g_aiT16_vec[16]), _mm_load_si128(&g_aiT16_vec[20])), _mm_unpacklo_epi16(_mm_load_si128(&g_aiT16_vec[24]), _mm_load_si128(&g_aiT16_vec[28]))));
  temp2 = _mm_unpacklo_epi64(_mm_unpackhi_epi32(_mm_unpacklo_epi16(_mm_load_si128(&g_aiT16_vec[0]), _mm_load_si128(&g_aiT16_vec[4])), _mm_unpacklo_epi16(_mm_load_si128(&g_aiT16_vec[8]), _mm_load_si128(&g_aiT16_vec[12]))),
                _mm_unpackhi_epi32(_mm_unpacklo_epi16(_mm_load_si128(&g_aiT16_vec[16]), _mm_load_si128(&g_aiT16_vec[20])), _mm_unpacklo_epi16(_mm_load_si128(&g_aiT16_vec[24]), _mm_load_si128(&g_aiT16_vec[28]))));
  temp3 = _mm_unpackhi_epi64(_mm_unpackhi_epi32(_mm_unpacklo_epi16(_mm_load_si128(&g_aiT16_vec[0]), _mm_load_si128(&g_aiT16_vec[4])), _mm_unpacklo_epi16(_mm_load_si128(&g_aiT16_vec[8]), _mm_load_si128(&g_aiT16_vec[12]))),
                _mm_unpackhi_epi32(_mm_unpacklo_epi16(_mm_load_si128(&g_aiT16_vec[16]), _mm_load_si128(&g_aiT16_vec[20])), _mm_unpacklo_epi16(_mm_load_si128(&g_aiT16_vec[24]), _mm_load_si128(&g_aiT16_vec[28]))));

  //store the transposed matrix in gaiT_vec
  _mm_store_si128(&gaiT_vec[8],  temp0);//g[0][0],g[2][0],g[4][0],g[6][0],...g[12][0],g[14][0]
  _mm_store_si128(&gaiT_vec[9],  temp1);//g[0][1],g[2][1]...g[14][1] 
  _mm_store_si128(&gaiT_vec[10], temp2);//g[0][2],g[2][2],g[4][2]...g[12][2],g[14][2]
  _mm_store_si128(&gaiT_vec[11], temp3);//g[0][3].g[2][3],...g[12][3],g[14][3]

  //===End of storing even row, only LHS 4 ====
  //====End of loading everything needed from g_aiT matrix

  /*========================================== 
       Start loading and transposing
       odd row of src
       All these src are used for O[k] loop
    ==========================================*/

  //====== First 8 elements of the row

  temp0 = _mm_unpacklo_epi64(_mm_unpacklo_epi32(_mm_unpacklo_epi16(_mm_load_si128(&in_vec[2]), _mm_load_si128(&in_vec[6])), _mm_unpacklo_epi16(_mm_load_si128(&in_vec[10]), _mm_load_si128(&in_vec[14]))),
                _mm_unpacklo_epi32(_mm_unpacklo_epi16(_mm_load_si128(&in_vec[18]), _mm_load_si128(&in_vec[22])), _mm_unpacklo_epi16(_mm_load_si128(&in_vec[26]), _mm_load_si128(&in_vec[30]))));  //a0b0c0d0e0f0g0h0
  temp1 = _mm_unpackhi_epi64(_mm_unpacklo_epi32(_mm_unpacklo_epi16(_mm_load_si128(&in_vec[2]), _mm_load_si128(&in_vec[6])), _mm_unpacklo_epi16(_mm_load_si128(&in_vec[10]), _mm_load_si128(&in_vec[14]))),
                _mm_unpacklo_epi32(_mm_unpacklo_epi16(_mm_load_si128(&in_vec[18]), _mm_load_si128(&in_vec[22])), _mm_unpacklo_epi16(_mm_load_si128(&in_vec[26]), _mm_load_si128(&in_vec[30]))));  //a1b1c1d1e1f1g1h1
  temp2 = _mm_unpacklo_epi64(_mm_unpackhi_epi32(_mm_unpacklo_epi16(_mm_load_si128(&in_vec[2]), _mm_load_si128(&in_vec[6])), _mm_unpacklo_epi16(_mm_load_si128(&in_vec[10]), _mm_load_si128(&in_vec[14]))),
                _mm_unpackhi_epi32(_mm_unpacklo_epi16(_mm_load_si128(&in_vec[18]), _mm_load_si128(&in_vec[22])), _mm_unpacklo_epi16(_mm_load_si128(&in_vec[26]), _mm_load_si128(&in_vec[30]))));
  temp3 = _mm_unpackhi_epi64(_mm_unpackhi_epi32(_mm_unpacklo_epi16(_mm_load_si128(&in_vec[2]), _mm_load_si128(&in_vec[6])), _mm_unpacklo_epi16(_mm_load_si128(&in_vec[10]), _mm_load_si128(&in_vec[14]))),
                _mm_unpackhi_epi32(_mm_unpacklo_epi16(_mm_load_si128(&in_vec[18]), _mm_load_si128(&in_vec[22])), _mm_unpacklo_epi16(_mm_load_si128(&in_vec[26]), _mm_load_si128(&in_vec[30]))));
  temp4 = _mm_unpacklo_epi64(_mm_unpacklo_epi32(_mm_unpackhi_epi16(_mm_load_si128(&in_vec[2]), _mm_load_si128(&in_vec[6])), _mm_unpackhi_epi16(_mm_load_si128(&in_vec[10]), _mm_load_si128(&in_vec[14]))),
                _mm_unpacklo_epi32(_mm_unpackhi_epi16(_mm_load_si128(&in_vec[18]), _mm_load_si128(&in_vec[22])), _mm_unpackhi_epi16(_mm_load_si128(&in_vec[26]), _mm_load_si128(&in_vec[30]))));
  temp5 = _mm_unpackhi_epi64(_mm_unpacklo_epi32(_mm_unpackhi_epi16(_mm_load_si128(&in_vec[2]), _mm_load_si128(&in_vec[6])), _mm_unpackhi_epi16(_mm_load_si128(&in_vec[10]), _mm_load_si128(&in_vec[14]))),
                _mm_unpacklo_epi32(_mm_unpackhi_epi16(_mm_load_si128(&in_vec[18]), _mm_load_si128(&in_vec[22])), _mm_unpackhi_epi16(_mm_load_si128(&in_vec[26]), _mm_load_si128(&in_vec[30]))));
  temp6 = _mm_unpacklo_epi64(_mm_unpackhi_epi32(_mm_unpackhi_epi16(_mm_load_si128(&in_vec[2]), _mm_load_si128(&in_vec[6])), _mm_unpackhi_epi16(_mm_load_si128(&in_vec[10]), _mm_load_si128(&in_vec[14]))),
                _mm_unpackhi_epi32(_mm_unpackhi_epi16(_mm_load_si128(&in_vec[18]), _mm_load_si128(&in_vec[22])), _mm_unpackhi_epi16(_mm_load_si128(&in_vec[26]), _mm_load_si128(&in_vec[30]))));
  temp7 = _mm_unpackhi_epi64(_mm_unpackhi_epi32(_mm_unpackhi_epi16(_mm_load_si128(&in_vec[2]), _mm_load_si128(&in_vec[6])), _mm_unpackhi_epi16(_mm_load_si128(&in_vec[10]), _mm_load_si128(&in_vec[14]))),
                _mm_unpackhi_epi32(_mm_unpackhi_epi16(_mm_load_si128(&in_vec[18]), _mm_load_si128(&in_vec[22])), _mm_unpackhi_epi16(_mm_load_si128(&in_vec[26]), _mm_load_si128(&in_vec[30]))));

  //store the transposed matrix in source_vec
  _mm_store_si128(&source_vec[0], temp0);
  _mm_store_si128(&source_vec[1], temp1);
  _mm_store_si128(&source_vec[2], temp2);
  _mm_store_si128(&source_vec[3], temp3);
  _mm_store_si128(&source_vec[4], temp4);
  _mm_store_si128(&source_vec[5], temp5);
  _mm_store_si128(&source_vec[6], temp6);
  _mm_store_si128(&source_vec[7], temp7); 

  //======End of transpose the first 8 of each odd row of src

  //=====the other 8 of each src odd row ======
  //Repeat what we did just now

  temp0 = _mm_unpacklo_epi64(_mm_unpacklo_epi32(_mm_unpacklo_epi16(_mm_load_si128(&in_vec[3]), _mm_load_si128(&in_vec[7])), _mm_unpacklo_epi16(_mm_load_si128(&in_vec[11]), _mm_load_si128(&in_vec[15]))),
                _mm_unpacklo_epi32(_mm_unpacklo_epi16(_mm_load_si128(&in_vec[19]), _mm_load_si128(&in_vec[23])), _mm_unpacklo_epi16(_mm_load_si128(&in_vec[27]), _mm_load_si128(&in_vec[31]))));
  temp1 = _mm_unpackhi_epi64(_mm_unpacklo_epi32(_mm_unpacklo_epi16(_mm_load_si128(&in_vec[3]), _mm_load_si128(&in_vec[7])), _mm_unpacklo_epi16(_mm_load_si128(&in_vec[11]), _mm_load_si128(&in_vec[15]))),
                _mm_unpacklo_epi32(_mm_unpacklo_epi16(_mm_load_si128(&in_vec[19]), _mm_load_si128(&in_vec[23])), _mm_unpacklo_epi16(_mm_load_si128(&in_vec[27]), _mm_load_si128(&in_vec[31]))));
  temp2 = _mm_unpacklo_epi64(_mm_unpackhi_epi32(_mm_unpacklo_epi16(_mm_load_si128(&in_vec[3]), _mm_load_si128(&in_vec[7])), _mm_unpacklo_epi16(_mm_load_si128(&in_vec[11]), _mm_load_si128(&in_vec[15]))),
                _mm_unpackhi_epi32(_mm_unpacklo_epi16(_mm_load_si128(&in_vec[19]), _mm_load_si128(&in_vec[23])), _mm_unpacklo_epi16(_mm_load_si128(&in_vec[27]), _mm_load_si128(&in_vec[31]))));
  temp3 = _mm_unpackhi_epi64(_mm_unpackhi_epi32(_mm_unpacklo_epi16(_mm_load_si128(&in_vec[3]), _mm_load_si128(&in_vec[7])), _mm_unpacklo_epi16(_mm_load_si128(&in_vec[11]), _mm_load_si128(&in_vec[15]))),
                _mm_unpackhi_epi32(_mm_unpacklo_epi16(_mm_load_si128(&in_vec[19]), _mm_load_si128(&in_vec[23])), _mm_unpacklo_epi16(_mm_load_si128(&in_vec[27]), _mm_load_si128(&in_vec[31]))));
  temp4 = _mm_unpacklo_epi64(_mm_unpacklo_epi32(_mm_unpackhi_epi16(_mm_load_si128(&in_vec[3]), _mm_load_si128(&in_vec[7])), _mm_unpackhi_epi16(_mm_load_si128(&in_vec[11]), _mm_load_si128(&in_vec[15]))),
                _mm_unpacklo_epi32(_mm_unpackhi_epi16(_mm_load_si128(&in_vec[19]), _mm_load_si128(&in_vec[23])), _mm_unpackhi_epi16(_mm_load_si128(&in_vec[27]), _mm_load_si128(&in_vec[31]))));
  temp5 = _mm_unpackhi_epi64(_mm_unpacklo_epi32(_mm_unpackhi_epi16(_mm_load_si128(&in_vec[3]), _mm_load_si128(&in_vec[7])), _mm_unpackhi_epi16(_mm_load_si128(&in_vec[11]), _mm_load_si128(&in_vec[15]))),
                _mm_unpacklo_epi32(_mm_unpackhi_epi16(_mm_load_si128(&in_vec[19]), _mm_load_si128(&in_vec[23])), _mm_unpackhi_epi16(_mm_load_si128(&in_vec[27]), _mm_load_si128(&in_vec[31]))));
  temp6 = _mm_unpacklo_epi64(_mm_unpackhi_epi32(_mm_unpackhi_epi16(_mm_load_si128(&in_vec[3]), _mm_load_si128(&in_vec[7])), _mm_unpackhi_epi16(_mm_load_si128(&in_vec[11]), _mm_load_si128(&in_vec[15]))),
                _mm_unpackhi_epi32(_mm_unpackhi_epi16(_mm_load_si128(&in_vec[19]), _mm_load_si128(&in_vec[23])), _mm_unpackhi_epi16(_mm_load_si128(&in_vec[27]), _mm_load_si128(&in_vec[31]))));
  temp7 = _mm_unpackhi_epi64(_mm_unpackhi_epi32(_mm_unpackhi_epi16(_mm_load_si128(&in_vec[3]), _mm_load_si128(&in_vec[7])), _mm_unpackhi_epi16(_mm_load_si128(&in_vec[11]), _mm_load_si128(&in_vec[15]))),
                _mm_unpackhi_epi32(_mm_unpackhi_epi16(_mm_load_si128(&in_vec[19]), _mm_load_si128(&in_vec[23])), _mm_unpackhi_epi16(_mm_load_si128(&in_vec[27]), _mm_load_si128(&in_vec[31])))); 

  //store the LHS 8 elements of each odd row
  _mm_store_si128(&source_vec[8], temp0);  
  _mm_store_si128(&source_vec[9], temp1);
  _mm_store_si128(&source_vec[10], temp2);
  _mm_store_si128(&source_vec[11], temp3);
  _mm_store_si128(&source_vec[12], temp4);
  _mm_store_si128(&source_vec[13], temp5);
  _mm_store_si128(&source_vec[14], temp6);
  _mm_store_si128(&source_vec[15], temp7);

  /*========================================== 
       End of loading and transposing 
       odd row of src
       All these src are used for O[k]
    ==========================================*/

  /*========================================== 
       Start loading and transposing
       even row of src
       All these src are used for EO[k],EE,EEE,EEO[k]
  ==========================================*/
  //====== First 8 elements of each even row of src
  temp0 = _mm_unpacklo_epi64(_mm_unpacklo_epi32(_mm_unpacklo_epi16(_mm_load_si128(&in_vec[0]), _mm_load_si128(&in_vec[4])), _mm_unpacklo_epi16(_mm_load_si128(&in_vec[8]), _mm_load_si128(&in_vec[12]))),
                _mm_unpacklo_epi32(_mm_unpacklo_epi16(_mm_load_si128(&in_vec[16]), _mm_load_si128(&in_vec[20])), _mm_unpacklo_epi16(_mm_load_si128(&in_vec[24]), _mm_load_si128(&in_vec[28]))));  //a0b0c0d0e0f0g0h0
  temp1 = _mm_unpackhi_epi64(_mm_unpacklo_epi32(_mm_unpacklo_epi16(_mm_load_si128(&in_vec[0]), _mm_load_si128(&in_vec[4])), _mm_unpacklo_epi16(_mm_load_si128(&in_vec[8]), _mm_load_si128(&in_vec[12]))),
                _mm_unpacklo_epi32(_mm_unpacklo_epi16(_mm_load_si128(&in_vec[16]), _mm_load_si128(&in_vec[20])), _mm_unpacklo_epi16(_mm_load_si128(&in_vec[24]), _mm_load_si128(&in_vec[28]))));  //a1b1c1d1e1f1g1h1
  temp2 = _mm_unpacklo_epi64(_mm_unpackhi_epi32(_mm_unpacklo_epi16(_mm_load_si128(&in_vec[0]), _mm_load_si128(&in_vec[4])), _mm_unpacklo_epi16(_mm_load_si128(&in_vec[8]), _mm_load_si128(&in_vec[12]))),
                _mm_unpackhi_epi32(_mm_unpacklo_epi16(_mm_load_si128(&in_vec[16]), _mm_load_si128(&in_vec[20])), _mm_unpacklo_epi16(_mm_load_si128(&in_vec[24]), _mm_load_si128(&in_vec[28]))));
  temp3 = _mm_unpackhi_epi64(_mm_unpackhi_epi32(_mm_unpacklo_epi16(_mm_load_si128(&in_vec[0]), _mm_load_si128(&in_vec[4])), _mm_unpacklo_epi16(_mm_load_si128(&in_vec[8]), _mm_load_si128(&in_vec[12]))),
                _mm_unpackhi_epi32(_mm_unpacklo_epi16(_mm_load_si128(&in_vec[16]), _mm_load_si128(&in_vec[20])), _mm_unpacklo_epi16(_mm_load_si128(&in_vec[24]), _mm_load_si128(&in_vec[28]))));
  temp4 = _mm_unpacklo_epi64(_mm_unpacklo_epi32(_mm_unpackhi_epi16(_mm_load_si128(&in_vec[0]), _mm_load_si128(&in_vec[4])), _mm_unpackhi_epi16(_mm_load_si128(&in_vec[8]), _mm_load_si128(&in_vec[12]))),
                _mm_unpacklo_epi32(_mm_unpackhi_epi16(_mm_load_si128(&in_vec[16]), _mm_load_si128(&in_vec[20])), _mm_unpackhi_epi16(_mm_load_si128(&in_vec[24]), _mm_load_si128(&in_vec[28]))));
  temp5 = _mm_unpackhi_epi64(_mm_unpacklo_epi32(_mm_unpackhi_epi16(_mm_load_si128(&in_vec[0]), _mm_load_si128(&in_vec[4])), _mm_unpackhi_epi16(_mm_load_si128(&in_vec[8]), _mm_load_si128(&in_vec[12]))),
                _mm_unpacklo_epi32(_mm_unpackhi_epi16(_mm_load_si128(&in_vec[16]), _mm_load_si128(&in_vec[20])), _mm_unpackhi_epi16(_mm_load_si128(&in_vec[24]), _mm_load_si128(&in_vec[28]))));
  temp6 = _mm_unpacklo_epi64(_mm_unpackhi_epi32(_mm_unpackhi_epi16(_mm_load_si128(&in_vec[0]), _mm_load_si128(&in_vec[4])), _mm_unpackhi_epi16(_mm_load_si128(&in_vec[8]), _mm_load_si128(&in_vec[12]))),
                _mm_unpackhi_epi32(_mm_unpackhi_epi16(_mm_load_si128(&in_vec[16]), _mm_load_si128(&in_vec[20])), _mm_unpackhi_epi16(_mm_load_si128(&in_vec[24]), _mm_load_si128(&in_vec[28]))));
  temp7 = _mm_unpackhi_epi64(_mm_unpackhi_epi32(_mm_unpackhi_epi16(_mm_load_si128(&in_vec[0]), _mm_load_si128(&in_vec[4])), _mm_unpackhi_epi16(_mm_load_si128(&in_vec[8]), _mm_load_si128(&in_vec[12]))),
                _mm_unpackhi_epi32(_mm_unpackhi_epi16(_mm_load_si128(&in_vec[16]), _mm_load_si128(&in_vec[20])), _mm_unpackhi_epi16(_mm_load_si128(&in_vec[24]), _mm_load_si128(&in_vec[28]))));

  //store the transposed matrix in source_vec

  _mm_store_si128(&source_vec[16], temp0);
  _mm_store_si128(&source_vec[17], temp1);
  _mm_store_si128(&source_vec[18], temp2);
  _mm_store_si128(&source_vec[19], temp3);
  _mm_store_si128(&source_vec[20], temp4);
  _mm_store_si128(&source_vec[21], temp5);
  _mm_store_si128(&source_vec[22], temp6);
  _mm_store_si128(&source_vec[23], temp7);

  //======End of transpose the first 8 of each even row of src

  //=====the other 8 of each src even row ======

  temp0 = _mm_unpacklo_epi64(_mm_unpacklo_epi32(_mm_unpacklo_epi16(_mm_load_si128(&in_vec[1]), _mm_load_si128(&in_vec[5])), _mm_unpacklo_epi16(_mm_load_si128(&in_vec[9]), _mm_load_si128(&in_vec[13]))),
          _mm_unpacklo_epi32(_mm_unpacklo_epi16(_mm_load_si128(&in_vec[17]), _mm_load_si128(&in_vec[21])), _mm_unpacklo_epi16(_mm_load_si128(&in_vec[25]), _mm_load_si128(&in_vec[29]))));
  temp1 = _mm_unpackhi_epi64(_mm_unpacklo_epi32(_mm_unpacklo_epi16(_mm_load_si128(&in_vec[1]), _mm_load_si128(&in_vec[5])), _mm_unpacklo_epi16(_mm_load_si128(&in_vec[9]), _mm_load_si128(&in_vec[13]))),
          _mm_unpacklo_epi32(_mm_unpacklo_epi16(_mm_load_si128(&in_vec[17]), _mm_load_si128(&in_vec[21])), _mm_unpacklo_epi16(_mm_load_si128(&in_vec[25]), _mm_load_si128(&in_vec[29]))));
  temp2 = _mm_unpacklo_epi64(_mm_unpackhi_epi32(_mm_unpacklo_epi16(_mm_load_si128(&in_vec[1]), _mm_load_si128(&in_vec[5])), _mm_unpacklo_epi16(_mm_load_si128(&in_vec[9]), _mm_load_si128(&in_vec[13]))),
          _mm_unpackhi_epi32(_mm_unpacklo_epi16(_mm_load_si128(&in_vec[17]), _mm_load_si128(&in_vec[21])), _mm_unpacklo_epi16(_mm_load_si128(&in_vec[25]), _mm_load_si128(&in_vec[29]))));
  temp3 = _mm_unpackhi_epi64(_mm_unpackhi_epi32(_mm_unpacklo_epi16(_mm_load_si128(&in_vec[1]), _mm_load_si128(&in_vec[5])), _mm_unpacklo_epi16(_mm_load_si128(&in_vec[9]), _mm_load_si128(&in_vec[13]))),
          _mm_unpackhi_epi32(_mm_unpacklo_epi16(_mm_load_si128(&in_vec[17]), _mm_load_si128(&in_vec[21])), _mm_unpacklo_epi16(_mm_load_si128(&in_vec[25]), _mm_load_si128(&in_vec[29]))));
  temp4 = _mm_unpacklo_epi64(_mm_unpacklo_epi32(_mm_unpackhi_epi16(_mm_load_si128(&in_vec[1]), _mm_load_si128(&in_vec[5])), _mm_unpackhi_epi16(_mm_load_si128(&in_vec[9]), _mm_load_si128(&in_vec[13]))),
          _mm_unpacklo_epi32(_mm_unpackhi_epi16(_mm_load_si128(&in_vec[17]), _mm_load_si128(&in_vec[21])), _mm_unpackhi_epi16(_mm_load_si128(&in_vec[25]), _mm_load_si128(&in_vec[29]))));
  temp5 = _mm_unpackhi_epi64(_mm_unpacklo_epi32(_mm_unpackhi_epi16(_mm_load_si128(&in_vec[1]), _mm_load_si128(&in_vec[5])), _mm_unpackhi_epi16(_mm_load_si128(&in_vec[9]), _mm_load_si128(&in_vec[13]))),
          _mm_unpacklo_epi32(_mm_unpackhi_epi16(_mm_load_si128(&in_vec[17]), _mm_load_si128(&in_vec[21])), _mm_unpackhi_epi16(_mm_load_si128(&in_vec[25]), _mm_load_si128(&in_vec[29]))));
  temp6 = _mm_unpacklo_epi64(_mm_unpackhi_epi32(_mm_unpackhi_epi16(_mm_load_si128(&in_vec[1]), _mm_load_si128(&in_vec[5])), _mm_unpackhi_epi16(_mm_load_si128(&in_vec[9]), _mm_load_si128(&in_vec[13]))),
          _mm_unpackhi_epi32(_mm_unpackhi_epi16(_mm_load_si128(&in_vec[17]), _mm_load_si128(&in_vec[21])), _mm_unpackhi_epi16(_mm_load_si128(&in_vec[25]), _mm_load_si128(&in_vec[29]))));
  temp7 = _mm_unpackhi_epi64(_mm_unpackhi_epi32(_mm_unpackhi_epi16(_mm_load_si128(&in_vec[1]), _mm_load_si128(&in_vec[5])), _mm_unpackhi_epi16(_mm_load_si128(&in_vec[9]), _mm_load_si128(&in_vec[13]))),
          _mm_unpackhi_epi32(_mm_unpackhi_epi16(_mm_load_si128(&in_vec[17]), _mm_load_si128(&in_vec[21])), _mm_unpackhi_epi16(_mm_load_si128(&in_vec[25]), _mm_load_si128(&in_vec[29])))); 

  //store the LHS 8 elements of each even row
  _mm_store_si128(&source_vec[24], temp0);  
  _mm_store_si128(&source_vec[25], temp1);
  _mm_store_si128(&source_vec[26], temp2);
  _mm_store_si128(&source_vec[27], temp3);
  _mm_store_si128(&source_vec[28], temp4);
  _mm_store_si128(&source_vec[29], temp5);
  _mm_store_si128(&source_vec[30], temp6);
  _mm_store_si128(&source_vec[31], temp7);

  /*========================================== 
       End of loading and transposing
       even row of src
       All these src are used for O[k]
  ==========================================*/

  for (int j=0; j<16; j++)
  {
    /* Utilizing symmetry properties to the maximum to minimize the number of multiplications */
    __m128i source_O = _mm_load_si128 (&source_vec[j]);
    //odd row first 8 elements of g_aiT
    //k=0 
    __m128i O_0 = _mm_madd_epi16 (_mm_load_si128 (&gaiT_vec[0]), source_O);//g[1][0]*s[1][j]+g[3][0]*s[3][j],g[5]*s[5]+g[7]*s[7]...+g[13]*+g[15]*s[15][j]    
    //k=1   
    __m128i O_1 = _mm_madd_epi16 (_mm_load_si128 (&gaiT_vec[1]), source_O);//g[1][1]*s[1][j]+g[3][1]*s[3][j],  ...  
    //k=2
    __m128i O_2 = _mm_madd_epi16 (_mm_load_si128 (&gaiT_vec[2]), source_O);
    //k=3 
    __m128i O_3 = _mm_madd_epi16 (_mm_load_si128 (&gaiT_vec[3]), source_O);
    //k=4
    __m128i O_4 = _mm_madd_epi16 (_mm_load_si128 (&gaiT_vec[4]), source_O);
    //k=5
    __m128i O_5 = _mm_madd_epi16 (_mm_load_si128 (&gaiT_vec[5]), source_O);
    //k=6
    __m128i O_6 = _mm_madd_epi16 (_mm_load_si128 (&gaiT_vec[6]), source_O);  
    //k=7  
    __m128i O_7 = _mm_madd_epi16 (_mm_load_si128 (&gaiT_vec[7]), source_O);//g[1][7]*s[1][j]+g[3][7]*s[3][j]
    __m128i O_zero2three = _mm_hadd_epi32 (_mm_hadd_epi32 (O_0, O_1), _mm_hadd_epi32 (O_2, O_3));//O[0]~O[3]
    __m128i O_four2seven = _mm_hadd_epi32 (_mm_hadd_epi32 (O_4, O_5), _mm_hadd_epi32 (O_6, O_7));//O[4]~O[7]
    //===End of first part of O[K],  g_aiT16[1][k]*src[16],k=0-7


    //even row of src
    __m128i source_EE = _mm_load_si128 (&source_vec[j+16]); //src[0][j], src[2][j],src[4][j],...,src[14][j]
    
    //put the pair of Lo and Hi next to each other
    __m128i loEE = _mm_unpacklo_epi16(_mm_mullo_epi16 (_mm_load_si128 (&gaiT_vec[8]), source_EE), _mm_mulhi_epi16 (_mm_load_si128 (&gaiT_vec[8]), source_EE));//g[0][0]*s[0][j]lo+hi... g[6][0]*s[6][j]lo+hi
    __m128i hiEE = _mm_unpackhi_epi16(_mm_mullo_epi16 (_mm_load_si128 (&gaiT_vec[8]), source_EE), _mm_mulhi_epi16 (_mm_load_si128 (&gaiT_vec[8]), source_EE));//g[8][0]*s[8][j]lo+hi... g[14][0]*s[14][j]lo+hi
    __m128i temp0_zero2fourteen_add = _mm_add_epi32(loEE, hiEE);//g[0][0]*+g[8][0]*,g[2][0]*+g[10][0]*,g[4][0]+g[12][0],g[6][0]+g[14][0]
    __m128i temp0_zero2fourteen_sub = _mm_hsub_epi32 (loEE, hiEE);//g[0][0]-g[2][0],g[4][0]-g[6][0],g[8][0]-g[10][0],g[12][0]-g[14][0]

    //repeat same pattern as gaiT_vec[8]
    loEE = _mm_unpacklo_epi16(_mm_mullo_epi16 (_mm_load_si128 (&gaiT_vec[9]), source_EE), _mm_mulhi_epi16 (_mm_load_si128 (&gaiT_vec[9]), source_EE));
    hiEE = _mm_unpackhi_epi16(_mm_mullo_epi16 (_mm_load_si128 (&gaiT_vec[9]), source_EE), _mm_mulhi_epi16 (_mm_load_si128 (&gaiT_vec[9]), source_EE));    
    __m128i temp1_zero2fourteen_add = _mm_add_epi32(loEE, hiEE); //g[0][1]*+g[8][1]*,g[2][1]*+g[10][1]*,g[4][1]+g[12][1],g[6][1]+g[14][1] 
    __m128i temp1_zero2fourteen_sub = _mm_hsub_epi32 (loEE, hiEE);//g[0][1]-g[2][1],g[4][1]-g[6][1],g[8][1]-g[10][1],g[12][1]-g[14][1]   
             

    loEE = _mm_unpacklo_epi16(_mm_mullo_epi16 (_mm_load_si128 (&gaiT_vec[10]), source_EE), _mm_mulhi_epi16 (_mm_load_si128 (&gaiT_vec[10]), source_EE));
    hiEE = _mm_unpackhi_epi16(_mm_mullo_epi16 (_mm_load_si128 (&gaiT_vec[10]), source_EE), _mm_mulhi_epi16 (_mm_load_si128 (&gaiT_vec[10]), source_EE));    
    __m128i temp2_zero2fourteen_add = _mm_add_epi32(loEE, hiEE);  
    __m128i temp2_zero2fourteen_sub = _mm_hsub_epi32 (loEE, hiEE);//g[0][2]-g[2][2],g[4][2]-g[6][2],g[8][2]-g[10][2],g[12][2]-g[14][2]    


    loEE = _mm_unpacklo_epi16(_mm_mullo_epi16 (_mm_load_si128 (&gaiT_vec[11]), source_EE), _mm_mulhi_epi16 (_mm_load_si128 (&gaiT_vec[11]), source_EE));
    hiEE = _mm_unpackhi_epi16(_mm_mullo_epi16 (_mm_load_si128 (&gaiT_vec[11]), source_EE), _mm_mulhi_epi16 (_mm_load_si128 (&gaiT_vec[11]), source_EE));    
    __m128i temp3_zero2fourteen_add = _mm_add_epi32(loEE,hiEE);  
    __m128i temp3_zero2fourteen_sub = _mm_hsub_epi32 (loEE, hiEE);//g[0][3]-g[2][3],g[4][3]-g[6][3],g[8][3]-g[10][3],g[12][3]-g[14][3]   


    //Then add all 4 cycles together
    __m128i E_add =_mm_hadd_epi32 (_mm_hadd_epi32 (temp0_zero2fourteen_add, temp1_zero2fourteen_add), _mm_hadd_epi32 (temp2_zero2fourteen_add, temp3_zero2fourteen_add)); //g[0][0]+g[2][0]+g[4][0]+g[6][0]+g[8][0]+g[10][0]+g[12][0]+g[14][0],
                                           //g[0][1]+g[2][1]+g[4][1]+g[6][1]+g[8][1]+g[10][1]+g[12][1]+g[14][1]
                                           //g[0][2]+g[2][2]+g[4][2]+g[6][2]+g[8][2]+g[10][2]+g[12][2]+g[14][2],
                                           //g[0][3]+g[2][3]+g[4][3]+g[6][3]+g[8][3]+g[10][3]+g[12][3]+g[14][3]


      

    __m128i E_sub =_mm_hadd_epi32 (_mm_hadd_epi32 (temp3_zero2fourteen_sub, temp2_zero2fourteen_sub), _mm_hadd_epi32 (temp1_zero2fourteen_sub, temp0_zero2fourteen_sub));//g[0][3]-g[2][3]+g[4][3]-g[6][3]+g[8][3]-g[10][3]+g[12][3]-g[14][3],
                                          //g[0][2]-g[2][2]+g[4][2]-g[6][2]+g[8][2]-g[10][2]+g[12][2]-g[14][2]
                                          //g[0][1]-g[2][1]+g[4][1]-g[6][1]+g[8][1]-g[10][1]+g[12][1]-g[14][1]
                                          //g[0][0]-g[2][0]+g[4][0]-g[6][0]+g[8][0]-g[10][0]+g[12][0]-g[14][0]
    //dst[0] ~ dst[7] => E[k] + O[k] k=0~7
    __m128i add_vec = _mm_set1_epi32(add);
    __m128i sum_zero2three = _mm_add_epi32(E_add, O_zero2three); // k=0~3
    __m128i sum_four2seven = _mm_add_epi32(E_sub, O_four2seven); // k=4~7
    //Convert packed 32-bit integers from source 
    //to packed 16-bit integers using signed saturation, 
    //and store the results in dst.
    __m128i E_K_add_O_K = _mm_packs_epi32(_mm_srai_epi32(_mm_add_epi32(sum_zero2three, add_vec), shift), _mm_srai_epi32(_mm_add_epi32(sum_four2seven, add_vec), shift));
    
    //dst[8]~ dst[15] E[7-k] - O[7-k] k=0~7
    
    __m128i temp0_eight2fifteen = _mm_srai_epi32(_mm_add_epi32(_mm_sub_epi32(E_sub, O_four2seven), add_vec), shift); // shift right    

    __m128i temp1_eight2fifteen = _mm_srai_epi32(_mm_add_epi32(_mm_sub_epi32(E_add, O_zero2three), add_vec), shift); // shift right

    __m128i temp2_eight2fifteen = _mm_unpacklo_epi32(_mm_unpackhi_epi32(_mm_unpackhi_epi32(temp0_eight2fifteen, temp1_eight2fifteen), _mm_unpacklo_epi32(temp0_eight2fifteen, temp1_eight2fifteen)),
         _mm_unpacklo_epi32(_mm_unpackhi_epi32(temp0_eight2fifteen, temp1_eight2fifteen), _mm_unpacklo_epi32(temp0_eight2fifteen, temp1_eight2fifteen)));
         //g[0][0]-g[2][0]+g[4][0]-g[6][0]+g[8][0]-g[10][0]+g[12][0]-g[14][0]-g[1][7]-g[3][7]-..-g[15][7]
         //g[0][1]-g[2][1]+g[4][1]-g[6][1]+g[8][1]-g[10][1]+g[12][1]-g[14][1]-g[1][6]-g[3][6]-..-g[15][6]
         //g[0][2]-g[2][2]+g[4][2]-g[6][2]+g[8][2]-g[10][2]+g[12][2]-g[14][2]-g[1][5]-g[3][5]-..-g[15][5] 
         //g[0][3]-g[2][3]+g[4][3]-g[6][3]+g[8][3]-g[10][3]+g[12][3]-g[14][3]-g[1][4]-g[3][4]-..-g[15][4]

    __m128i temp3_eight2fifteen = _mm_unpackhi_epi32(_mm_unpackhi_epi32(_mm_unpackhi_epi32(temp0_eight2fifteen, temp1_eight2fifteen), _mm_unpacklo_epi32(temp0_eight2fifteen, temp1_eight2fifteen)),
         _mm_unpacklo_epi32(_mm_unpackhi_epi32(temp0_eight2fifteen, temp1_eight2fifteen), _mm_unpacklo_epi32(temp0_eight2fifteen, temp1_eight2fifteen)));
         //g[0][3]+g[2][3]+g[4][3]+g[6][3]+g[8][3]+g[10][3]+g[12][3]+g[14][3]-g[1][3]-g[3][3]-..-g[15][3]
         //g[0][2]+g[2][2]+g[4][2]+g[6][2]+g[8][2]+g[10][2]+g[12][2]+g[14][2]-g[1][2]-g[3][2]-..-g[15][2]
         //g[0][1]+g[2][1]+g[4][1]+g[6][1]+g[8][1]+g[10][1]+g[12][1]+g[14][1]-g[1][1]-g[3][1]-..-g[15][1]
         //g[0][0]+g[2][0]+g[4][0]+g[6][0]+g[8][0]+g[10][0]+g[12][0]+g[14][0]-g[1][0]-g[3][0]-..-g[15][0]
         
    __m128i E_K_sub_O_K = _mm_packs_epi32(temp2_eight2fifteen, temp3_eight2fifteen); 
         //g[0][0]-g[2][0]+g[4][0]-g[6][0]+g[8][0]-g[10][0]+g[12][0]-g[14][0]-g[1][7]-g[3][7]-..-g[15][7]
         //g[0][1]-g[2][1]+g[4][1]-g[6][1]+g[8][1]-g[10][1]+g[12][1]-g[14][1]-g[1][6]-g[3][6]-..-g[15][6]
         //g[0][2]-g[2][2]+g[4][2]-g[6][2]+g[8][2]-g[10][2]+g[12][2]-g[14][2]-g[1][5]-g[3][5]-..-g[15][5] 
         //g[0][3]-g[2][3]+g[4][3]-g[6][3]+g[8][3]-g[10][3]+g[12][3]-g[14][3]-g[1][4]-g[3][4]-..-g[15][4]
         //g[0][3]+g[2][3]+g[4][3]+g[6][3]+g[8][3]+g[10][3]+g[12][3]+g[14][3]-g[1][3]-g[3][3]-..-g[15][3]
         //g[0][2]+g[2][2]+g[4][2]+g[6][2]+g[8][2]+g[10][2]+g[12][2]+g[14][2]-g[1][2]-g[3][2]-..-g[15][2]
         //g[0][1]+g[2][1]+g[4][1]+g[6][1]+g[8][1]+g[10][1]+g[12][1]+g[14][1]-g[1][1]-g[3][1]-..-g[15][1]
         //g[0][0]+g[2][0]+g[4][0]+g[6][0]+g[8][0]+g[10][0]+g[12][0]+g[14][0]-g[1][0]-g[3][0]-..-g[15][0]
    _mm_store_si128 (&out_vec[2*j], E_K_add_O_K);// if j=0, out_vec[0] includes the first 8 elements of the row 0 of the output matrix
                                                //that is, j=0, out_vec[0]=dst[0],dst[1]...,dst[6],dst[7]
    _mm_store_si128 (&out_vec[2*j+1], E_K_sub_O_K);
    src ++;
    dst += 16;
  }
}
///end of simd codes!

static void idct16_simd(short* pCoeff, short* pDst)
{
  short tmp[ 16*16] __attribute__((aligned(16)));
  partialButterflyInverse16_simd(pCoeff, tmp, 7);
  partialButterflyInverse16_simd(tmp, pDst, 12);
}

int main(int argc, char **argv)
{
    //allocate memory 16-byte aligned
    short *scalar_input = (short*) memalign(16, IDCT_SIZE*IDCT_SIZE*sizeof(short));
    short *scalar_output = (short *) memalign(16, IDCT_SIZE*IDCT_SIZE*sizeof(short));

    short *simd_input = (short*) memalign(16, IDCT_SIZE*IDCT_SIZE*sizeof(short));
    short *simd_output = (short *) memalign(16, IDCT_SIZE*IDCT_SIZE*sizeof(short));

    //initialize input
    printf("input array:\n");
    for(int j=0;j<IDCT_SIZE;j++){
        for(int i=0;i<IDCT_SIZE;i++){
            short value = rand()%2 ? (rand()%32768) : -(rand()%32768) ;
            scalar_input[j*IDCT_SIZE+i] = value;
            simd_input  [j*IDCT_SIZE+i] = value;
	    printf("%d\t", value);
        }
        printf("\n");
    }

    idct16_scalar(scalar_input, scalar_output);
    idct16_simd  (simd_input  , simd_output);

    //check for correctness
    compare_results (scalar_output, simd_output, "scalar and simd");

    printf("output array:\n");
    for(int j=0;j<IDCT_SIZE;j++){
        for(int i=0;i<IDCT_SIZE;i++){
	    printf("%d\t", scalar_output[j*IDCT_SIZE+i]);
        }
        printf("\n");
    }

    //Measure the performance of each kernel
    benchmark (idct16_scalar, scalar_input, scalar_output, "scalar");
    benchmark (idct16_simd, simd_input, simd_output, "simd");

    //cleanup
    free(scalar_input);    free(scalar_output);
    free(simd_input); free(simd_output);
}
