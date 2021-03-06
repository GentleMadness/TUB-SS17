//scalar code for the inverse transform
static void simd(short *src, short *dst, int shift)
{
// SIMD Declaration of variables
  int E[8] __attribute__ ((aligned (16)));
  int O[8] __attribute__ ((aligned (16)));
  int EE[4] __attribute__ ((aligned (16)));
  int EO[4] __attribute__ ((aligned (16)));
  int EEE[4] __attribute__ ((aligned (16)));
  int EEO[4] __attribute__ ((aligned (16)));

  int add = 1<<(shift-1); //no SIMD replacement 
  
  // used to store new matrix:
  short gaiT[8][8] __attribute__ ((aligned (16)));
  short source[16][16] __attribute__ ((aligned (16)));
 
//cast the variables to SIMD vector type
  __m128i *in_vec = (__m128i *) src;       
  __m128i *out_vec = (__m128i *) dst;
  __m128i *g_aiT16_vec = (__m128i *)g_aiT16;
  __m128i *E_vec = (__m128i *)E;
  __m128i *O_vec = (__m128i *)O; 
  __m128i *EE_vec = (__m128i *)EE; 
  __m128i *EO_vec = (__m128i *)EO; 
  __m128i *EEE_vec = (__m128i *)EEE; 
  __m128i *EEO_vec = (__m128i *)EEO;
  
  __m128i *gaiT_vec = (__m128i *)gaiT;
  __m128i *source_vec = (__m128i *)source;
//load g_aiT16, short type in c is 16 bit, each vector is 128 bit

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

//step 1 load the first 8 elements of each odd row
__m128i a = _mm_load_si128(&g_aiT16_vec[2]);//first 8 num in row 1
__m128i b = _mm_load_si128(&g_aiT16_vec[6]);
__m128i c = _mm_load_si128(&g_aiT16_vec[10]);
__m128i d = _mm_load_si128(&g_aiT16_vec[14]);
__m128i e = _mm_load_si128(&g_aiT16_vec[18]);
__m128i f = _mm_load_si128(&g_aiT16_vec[22]);
__m128i g = _mm_load_si128(&g_aiT16_vec[26]);
__m128i h = _mm_load_si128(&g_aiT16_vec[30]);

//step 2 transpose every 16 bit
__m128i temp1 = _mm_unpacklo_epi16(a, b); 
__m128i temp2 = _mm_unpacklo_epi16(c, d);
__m128i temp3 = _mm_unpacklo_epi16(e, f);
__m128i temp4 = _mm_unpacklo_epi16(g, h);

__m128i temp5 = _mm_unpackhi_epi16(a, b);
__m128i temp6 = _mm_unpackhi_epi16(c, d);
__m128i temp7 = _mm_unpackhi_epi16(e, f);
__m128i temp8 = _mm_unpackhi_epi16(g, h);

//step 3 transpose every 32 bit
__m128i temp9  = _mm_unpacklo_epi32(temp1, temp2);
__m128i temp10 = _mm_unpackhi_epi32(temp1, temp2);
__m128i temp11 = _mm_unpacklo_epi32(temp3, temp4);
__m128i temp12 = _mm_unpackhi_epi32(temp3, temp4);
__m128i temp13 = _mm_unpacklo_epi32(temp5, temp6);
__m128i temp14 = _mm_unpackhi_epi32(temp5, temp6);
__m128i temp15 = _mm_unpacklo_epi32(temp7, temp8);
__m128i temp16 = _mm_unpackhi_epi32(temp7, temp8);
 
//step 4 transpose every 64 bit
a = _mm_unpacklo_epi64(temp9, temp11);  
b = _mm_unpackhi_epi64(temp9, temp11);
c = _mm_unpacklo_epi64(temp10, temp12);
d = _mm_unpackhi_epi64(temp10, temp12);
e = _mm_unpacklo_epi64(temp13, temp15);
f = _mm_unpackhi_epi64(temp13, temp15);
g = _mm_unpacklo_epi64(temp14, temp16);
h = _mm_unpackhi_epi64(temp14, temp16);

//step 5 store the transposed matrix 
_mm_store_si128(&gaiT_vec[0], a);
_mm_store_si128(&gaiT_vec[1], b); 
_mm_store_si128(&gaiT_vec[2], c);
_mm_store_si128(&gaiT_vec[3], d);
_mm_store_si128(&gaiT_vec[4], e);
_mm_store_si128(&gaiT_vec[5], f);
_mm_store_si128(&gaiT_vec[6], g);
_mm_store_si128(&gaiT_vec[7], h);
/*=== End of loading g_aiT odd row LHS 8*8 ===*/

//===Start loading g_aiT even row LHS 4 ==========
a = _mm_load_si128(&g_aiT16_vec[0]);//first 8 num in row 1
b = _mm_load_si128(&g_aiT16_vec[4]);// Actually only need 4 num of each row
c = _mm_load_si128(&g_aiT16_vec[8]);
d = _mm_load_si128(&g_aiT16_vec[12]);
e = _mm_load_si128(&g_aiT16_vec[16]);
f = _mm_load_si128(&g_aiT16_vec[20]);
g = _mm_load_si128(&g_aiT16_vec[24]);
h = _mm_load_si128(&g_aiT16_vec[38]);

//only first 4 of each line are needed
//only unpacklo is needed
temp1 = _mm_unpacklo_epi16(a, b);
temp2 = _mm_unpacklo_epi16(c, d);
temp3 = _mm_unpacklo_epi16(e, f);
temp4 = _mm_unpacklo_epi16(g, h);

temp5 = _mm_unpacklo_epi32(temp1, temp2);
temp6 = _mm_unpackhi_epi32(temp1, temp2);
temp7 = _mm_unpacklo_epi32(temp3, temp4);
temp8 = _mm_unpackhi_epi32(temp3, temp4);

T0 = _mm_unpacklo_epi64(temp5, temp7);  
T1 = _mm_unpackhi_epi64(temp5, temp7);
T2 = _mm_unpacklo_epi64(temp6, temp8);
T3 = _mm_unpackhi_epi64(temp6, temp8);

_mm_store_si128(&gaiT_vec[8],  T0);
_mm_store_si128(&gaiT_vec[9],  T1); 
_mm_store_si128(&gaiT_vec[10], T2);
_mm_store_si128(&gaiT_vec[11], T3);
//===End of storing even row, only LHS 4 ====
//====End of loading everything needed from g_aiT matrix

/*========================================== 
     Start loading and transposing
     odd row of src
     All these src are used for O[k] loop
  ==========================================*/

//====== First 8 elements of the row
a = _mm_load_si128(&in_vec[2]);
b = _mm_load_si128(&in_vec[6]);
c = _mm_load_si128(&in_vec[10]);
d = _mm_load_si128(&in_vec[14]);
e = _mm_load_si128(&in_vec[18]);
f = _mm_load_si128(&in_vec[22]);
g = _mm_load_si128(&in_vec[26]);
h = _mm_load_si128(&in_vec[30]);

temp1 = _mm_unpacklo_epi16(a, b); //a0b0a1b1a2b2a3b3
temp2 = _mm_unpacklo_epi16(c, d);
temp3 = _mm_unpacklo_epi16(e, f);
temp4 = _mm_unpacklo_epi16(g, h);
temp5 = _mm_unpackhi_epi16(a, b);
temp6 = _mm_unpackhi_epi16(c, d);
temp7 = _mm_unpackhi_epi16(e, f);
temp8 = _mm_unpackhi_epi16(g, h);

temp9  = _mm_unpacklo_epi32(temp1, temp2); //a0b0c0d0a1b1c1d1
temp10 = _mm_unpackhi_epi32(temp1, temp2);
temp11 = _mm_unpacklo_epi32(temp3, temp4); //e0f0g0h0e1f1g1h1
temp12 = _mm_unpackhi_epi32(temp3, temp4);
temp13 = _mm_unpacklo_epi32(temp5, temp6);
temp14 = _mm_unpackhi_epi32(temp5, temp6);
temp15 = _mm_unpacklo_epi32(temp7, temp8);
temp16 = _mm_unpackhi_epi32(temp7, temp8);
 
__m128i T0 = _mm_unpacklo_epi64(temp9, temp11);  //a0b0c0d0e0f0g0h0
__m128i T1 = _mm_unpackhi_epi64(temp9, temp11);  //a1b1c1d1e1f1g1h1
__m128i T2 = _mm_unpacklo_epi64(temp10, temp12);
__m128i T3 = _mm_unpackhi_epi64(temp10, temp12);
__m128i T4 = _mm_unpacklo_epi64(temp13, temp15);
__m128i T5 = _mm_unpackhi_epi64(temp13, temp15);
__m128i T6 = _mm_unpacklo_epi64(temp14, temp16);
__m128i T7 = _mm_unpackhi_epi64(temp14, temp16);

//store the transposed matrix in source_vec
_mm_store_si128(&source_vec[0], T0);//src [0] element of each odd row
_mm_store_si128(&source_vec[1], T1);//src [1] element of each odd row
_mm_store_si128(&source_vec[2], T2);
_mm_store_si128(&source_vec[3], T3);
_mm_store_si128(&source_vec[4], T4);
_mm_store_si128(&source_vec[5], T5);
_mm_store_si128(&source_vec[6], T6);
_mm_store_si128(&source_vec[7], T7);//src [7] element of each odd row  
//======End of transpose the first 8 of each odd row of src

//=====the other 8 of each src odd row ======
//Repeat what we did just now
a = _mm_load_si128(&in_vec[3]);
b = _mm_load_si128(&in_vec[7]);
c = _mm_load_si128(&in_vec[11]);
d = _mm_load_si128(&in_vec[15]);
e = _mm_load_si128(&in_vec[19]);
f = _mm_load_si128(&in_vec[23]);
g = _mm_load_si128(&in_vec[27]);
h = _mm_load_si128(&in_vec[31]);

temp1 = _mm_unpacklo_epi16(a, b);
temp2 = _mm_unpacklo_epi16(c, d);
temp3 = _mm_unpacklo_epi16(e, f);
temp4 = _mm_unpacklo_epi16(g, h);
temp5 = _mm_unpackhi_epi16(a, b);
temp6 = _mm_unpackhi_epi16(c, d);
temp7 = _mm_unpackhi_epi16(e, f);
temp8 = _mm_unpackhi_epi16(g, h);

temp9  = _mm_unpacklo_epi32(temp1, temp2);
temp10 = _mm_unpackhi_epi32(temp1, temp2);
temp11 = _mm_unpacklo_epi32(temp3, temp4);
temp12 = _mm_unpackhi_epi32(temp3, temp4);
temp13 = _mm_unpacklo_epi32(temp5, temp6);
temp14 = _mm_unpackhi_epi32(temp5, temp6);
temp15 = _mm_unpacklo_epi32(temp7, temp8);
temp16 = _mm_unpackhi_epi32(temp7, temp8);

T0 = _mm_unpacklo_epi64(temp9, temp11);
T1 = _mm_unpackhi_epi64(temp9, temp11);
T2 = _mm_unpacklo_epi64(temp10, temp12);
T3 = _mm_unpackhi_epi64(temp10, temp12);
T4 = _mm_unpacklo_epi64(temp13, temp15);
T5 = _mm_unpackhi_epi64(temp13, temp15);
T6 = _mm_unpacklo_epi64(temp14, temp16);
T7 = _mm_unpackhi_epi64(temp14, temp16); 

//store the LHS 8 elements of each odd row
_mm_store_si128(&source_vec[8], T0);  
_mm_store_si128(&source_vec[9], T1);
_mm_store_si128(&source_vec[10], T2);
_mm_store_si128(&source_vec[11], T3);
_mm_store_si128(&source_vec[12], T4);
_mm_store_si128(&source_vec[13], T5);
_mm_store_si128(&source_vec[14], T6);
_mm_store_si128(&source_vec[15], T7);
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
a = _mm_load_si128(&in_vec[0]);
b = _mm_load_si128(&in_vec[4]);
c = _mm_load_si128(&in_vec[8]);
d = _mm_load_si128(&in_vec[12]);
e = _mm_load_si128(&in_vec[16]);
f = _mm_load_si128(&in_vec[20]);
g = _mm_load_si128(&in_vec[24]);
h = _mm_load_si128(&in_vec[28]);

temp1 = _mm_unpacklo_epi16(a, b); //a0b0a1b1a2b2a3b3
temp2 = _mm_unpacklo_epi16(c, d);
temp3 = _mm_unpacklo_epi16(e, f);
temp4 = _mm_unpacklo_epi16(g, h);
temp5 = _mm_unpackhi_epi16(a, b);
temp6 = _mm_unpackhi_epi16(c, d);
temp7 = _mm_unpackhi_epi16(e, f);
temp8 = _mm_unpackhi_epi16(g, h);

temp9  = _mm_unpacklo_epi32(temp1, temp2); //a0b0c0d0a1b1c1d1
temp10 = _mm_unpackhi_epi32(temp1, temp2);
temp11 = _mm_unpacklo_epi32(temp3, temp4); //e0f0g0h0e1f1g1h1
temp12 = _mm_unpackhi_epi32(temp3, temp4);
temp13 = _mm_unpacklo_epi32(temp5, temp6);
temp14 = _mm_unpackhi_epi32(temp5, temp6);
temp15 = _mm_unpacklo_epi32(temp7, temp8);
temp16 = _mm_unpackhi_epi32(temp7, temp8);
 
T0 = _mm_unpacklo_epi64(temp9, temp11);  //a0b0c0d0e0f0g0h0
T1 = _mm_unpackhi_epi64(temp9, temp11);  //a1b1c1d1e1f1g1h1
T2 = _mm_unpacklo_epi64(temp10, temp12);
T3 = _mm_unpackhi_epi64(temp10, temp12);
T4 = _mm_unpacklo_epi64(temp13, temp15);
T5 = _mm_unpackhi_epi64(temp13, temp15);
T6 = _mm_unpacklo_epi64(temp14, temp16);
T7 = _mm_unpackhi_epi64(temp14, temp16);

//store the transposed matrix in source_vec
_mm_store_si128(&source_vec[16], T0);//src [0] element of each even row
_mm_store_si128(&source_vec[17], T1);//src [1] element of each even row
_mm_store_si128(&source_vec[18], T2);
_mm_store_si128(&source_vec[19], T3);
_mm_store_si128(&source_vec[20], T4);
_mm_store_si128(&source_vec[21], T5);
_mm_store_si128(&source_vec[22], T6);
_mm_store_si128(&source_vec[23], T7);//src [7] element of each even row  
//======End of transpose the first 8 of each even row of src

//=====the other 8 of each src even row ======
//Repeat what we did just now
a = _mm_load_si128(&in_vec[1]);
b = _mm_load_si128(&in_vec[5]);
c = _mm_load_si128(&in_vec[9]);
d = _mm_load_si128(&in_vec[13]);
e = _mm_load_si128(&in_vec[17]);
f = _mm_load_si128(&in_vec[21]);
g = _mm_load_si128(&in_vec[25]);
h = _mm_load_si128(&in_vec[29]);

temp1 = _mm_unpacklo_epi16(a, b);
temp2 = _mm_unpacklo_epi16(c, d);
temp3 = _mm_unpacklo_epi16(e, f);
temp4 = _mm_unpacklo_epi16(g, h);
temp5 = _mm_unpackhi_epi16(a, b);
temp6 = _mm_unpackhi_epi16(c, d);
temp7 = _mm_unpackhi_epi16(e, f);
temp8 = _mm_unpackhi_epi16(g, h);

temp9  = _mm_unpacklo_epi32(temp1, temp2);
temp10 = _mm_unpackhi_epi32(temp1, temp2);
temp11 = _mm_unpacklo_epi32(temp3, temp4);
temp12 = _mm_unpackhi_epi32(temp3, temp4);
temp13 = _mm_unpacklo_epi32(temp5, temp6);
temp14 = _mm_unpackhi_epi32(temp5, temp6);
temp15 = _mm_unpacklo_epi32(temp7, temp8);
temp16 = _mm_unpackhi_epi32(temp7, temp8);

T0 = _mm_unpacklo_epi64(temp9, temp11);
T1 = _mm_unpackhi_epi64(temp9, temp11);
T2 = _mm_unpacklo_epi64(temp10, temp12);
T3 = _mm_unpackhi_epi64(temp10, temp12);
T4 = _mm_unpacklo_epi64(temp13, temp15);
T5 = _mm_unpackhi_epi64(temp13, temp15);
T6 = _mm_unpacklo_epi64(temp14, temp16);
T7 = _mm_unpackhi_epi64(temp14, temp16); 

//store the LHS 8 elements of each even row
_mm_store_si128(&source_vec[24], T0);  
_mm_store_si128(&source_vec[25], T1);
_mm_store_si128(&source_vec[26], T2);
_mm_store_si128(&source_vec[27], T3);
_mm_store_si128(&source_vec[28], T4);
_mm_store_si128(&source_vec[29], T5);
_mm_store_si128(&source_vec[30], T6);
_mm_store_si128(&source_vec[31], T7);
/*========================================== 
     End of loading and transposing
     even row of src
     All these src are used for O[k]
==========================================*/

  for (int j=0; j<16; j++)
  {
    /* Utilizing symmetry properties to the maximum to minimize the number of multiplications */
    //odd row of src
    __m128i I0 = _mm_load_si128 (&source_vec[j]); 

/*  for (int k=0; k<8; k++)
    {
      O[k] = g_aiT16[ 1][k]*src[ 16] + g_aiT16[ 3][k]*src[ 3*16] + g_aiT16[ 5][k]*src[ 5*16] + g_aiT16[ 7][k]*src[ 7*16] +
        g_aiT16[ 9][k]*src[ 9*16] + g_aiT16[11][k]*src[11*16] + g_aiT16[13][k]*src[13*16] + g_aiT16[15][k]*src[15*16];
    }
*/
    //k=0
      //odd row first 8 elements of g_aiT
        __m128i I1 = _mm_load_si128 (&gaiT_vec[0]);   
        __m128i I2 = _mm_madd_epi16 (I1, I0);
      
    //k=1   
        __m128i I3 = _mm_load_si128 (&gaiT_vec[1]);   
        __m128i I4 = _mm_madd_epi16 (I3, I0);
    //k=2
        __m128i I5 = _mm_load_si128 (&gaiT_vec[2]);   
        __m128i I6 = _mm_madd_epi16 (I5, I0);
    //k=3
        __m128i I7 = _mm_load_si128 (&gaiT_vec[3]);   
        __m128i I8 = _mm_madd_epi16 (I7, I0);
    //k=4
        __m128i I9 = _mm_load_si128 (&gaiT_vec[4]);   
        __m128i I10 = _mm_madd_epi16 (I9, I0);
    //k=5
        __m128i I11 = _mm_load_si128 (&gaiT_vec[5]);   
        __m128i I12 = _mm_madd_epi16 (I11, I0);
    //k=6
        __m128i I13 = _mm_load_si128 (&gaiT_vec[6]);   
        __m128i I14 = _mm_madd_epi16 (I13, I0);  
    //k=7
        __m128i I15 = _mm_load_si128 (&gaiT_vec[7]);   
        __m128i I16 = _mm_madd_epi16 (I15, I0);
//===End of g_aiT16[1][k]*src[16],k=0-7

// Add all the products together:
       //S1= sum of k=0 and 1 
       __m128i S1 =_mm_hadd_epi32 (I2, I4);
       //S2= sum of k=2 and 3
       __m128i S2 =_mm_hadd_epi32 (I6, I8);

       //R1->sum of k=3|sum of k=2|sum of k=1|sum of k=0
       __m128i R1 =_mm_hadd_epi32 (S1, S2);

       __m128i S3 =_mm_hadd_epi32 (I10, I12);
       __m128i S4 =_mm_hadd_epi32 (I14, I16);
       //R2 = sum of k=7,6,5,4
       __m128i R2 =_mm_hadd_epi32 (S3, S4);
//======End of the 1st loop, O[k]=R2(K>3) R1(K<3) ==========

//=================Begin of the second loop=================    
/*
    for (int k=0; k<4; k++)
    { 
      EO[k] = g_aiT16[ 2][k]*src[ 2*16] + g_aiT16[ 6][k]*src[ 6*16] + 
              g_aiT16[10][k]*src[10*16] + g_aiT16[14][k]*src[14*16];
    }
*/
      //even row of src
      __m128i II0 = _mm_load_si128 (&source_vec[j+16]); 
      
       I1 = _mm_load_si128 (&gaiT_vec[8]);//No.0 of each even row 
       I2 = _mm_mullo_epi16 (I1, II0);//II0,odd row of source_vec
       I3 = _mm_mulhi_epi16 (I1, II0);
      //put the pair of Lo and Hi next to each other
      __m128i loI23 = _mm_unpacklo_epi16(I2,I3);
      __m128i hiI23 = _mm_unpackhi_epi16(I2,I3);
      //gt[0][0]*src[0][0]lo+hi, gt[2][0]*src[2][0]lo+hi,4,6   
      __m128i temp1 = _mm_add_epi32(loI23,hiI23);
      //gt[0][0]*src[0][0]-[2][0],(4-6),[8-10]
      __m128i temp2 = _mm_hsub_epi32 (loI23, hiI23);

      //repeat same pattern as gaiT_vec[8]
       I4 = _mm_load_si128 (&gaiT_vec[9]);       
       I5 = _mm_mullo_epi16 (I4, II0);
       I6 = _mm_mulhi_epi16 (I4, II0);
      __m128i loI56 = _mm_unpacklo_epi16(I5,I6);
      __m128i hiI56 = _mm_unpackhi_epi16(I5,I6);    
      __m128i temp3 = _mm_add_epi32(loI56,hiI56);  
      __m128i temp4 = _mm_hsub_epi32 (loI56, hiI56);   
             
       I7 = _mm_load_si128 (&gaiT_vec[10]);      
       I8 = _mm_mullo_epi16 (I7, II0);
       I9 = _mm_mulhi_epi16 (I7, II0);
      __m128i loI89 = _mm_unpacklo_epi16(I8,I9);
      __m128i hiI89 = _mm_unpackhi_epi16(I8,I9);    
      __m128i temp5 = _mm_add_epi32(loI89,hiI89);  
      __m128i temp6 = _mm_hsub_epi32 (loI89, hiI89);    

       I10 = _mm_load_si128 (&gaiT_vec[11]);       
       I11 = _mm_mullo_epi16 (I10, II0);
       I12 = _mm_mulhi_epi16 (I10, II0);
      __m128i loI1112 = _mm_unpacklo_epi16(I11,I12);
      __m128i hiI1112 = _mm_unpackhi_epi16(I11,I12);    
      __m128i temp7 = _mm_add_epi32(loI1112,hiI1112);  
      __m128i temp8 = _mm_hsub_epi32 (loI1112, hiI1112);   


//Then add all 4 cycles together
      __m128i A5 =_mm_hadd_epi32 (temp1, temp3);
      __m128i A6 =_mm_hadd_epi32 (temp5, temp7);
      __m128i R3 =_mm_hadd_epi32 (A5, A6);

      __m128i A7 =_mm_hadd_epi32 (temp8, temp6);
      __m128i A8 =_mm_hadd_epi32 (temp4, temp2);
      __m128i R4 =_mm_hadd_epi32 (A7, A8);
  
/*  for (int k=0; k<8; k++)
    {
      dst[k] = ...(E[k]   + O[k]   + add)>>shift ));
      dst[k+8]= ...(E[7-k] - O[7-k] + add)>>shift ));
    }
*/
      __m128i add_reg = _mm_set1_epi32(add);
      //O[k]+E[k]
      __m128i sum0_vec = _mm_add_epi32(R3,R1);        
      //O[k]+E[k]+add
      sum0_vec = _mm_add_epi32(sum0_vec,add_reg);
      //shift
      sum0_vec = _mm_srai_epi32(sum0_vec, shift);
   
      __m128i sum1_vec = _mm_add_epi32(R4,R2);
      sum1_vec = _mm_add_epi32(sum1_vec,add_reg);
      sum1_vec = _mm_srai_epi32(sum1_vec, shift); // shift right
    
    //Convert packed 32-bit integers from a and b 
    //to packed 16-bit integers using signed saturation, 
    //and store the results in dst.
      __m128i finalres0 = _mm_packs_epi32(sum0_vec, sum1_vec); // shrink packed 32bit to packed 16 bit and saturate
      _mm_store_si128 (&out_vec[2*j], finalres0);
  

      __m128i  sum2_vec = _mm_sub_epi32(R4, R2);
      sum2_vec = _mm_add_epi32(sum2_vec,add_reg);
      sum2_vec = _mm_srai_epi32(sum2_vec, shift); // shift right    

      __m128i sum3_vec = _mm_sub_epi32(R3, R1);
      sum3_vec = _mm_add_epi32(sum3_vec,add_reg);
      sum3_vec = _mm_srai_epi32(sum3_vec, shift); // shift right

         I1 = _mm_unpackhi_epi32(sum2_vec, sum3_vec);
         I2 = _mm_unpacklo_epi32(sum2_vec, sum3_vec);
         I3 = _mm_unpackhi_epi32(I1, I2);
         I4 = _mm_unpacklo_epi32(I1, I2);
         I5 = _mm_unpacklo_epi32(I3, I4);
         I6 = _mm_unpackhi_epi32(I3, I4);
         
         sum3_vec = _mm_packs_epi32(I5, I6); // shrink packed 32bit to packed 16 bit and saturate
         _mm_store_si128 (&out_vec[2*j+1], sum3_vec);

    src ++;
    dst += 16;
  }
}
