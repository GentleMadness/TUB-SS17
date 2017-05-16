;Yuting Fu 387512, Shan Kuan 387516

    .data
N_COEFFS:   .word   3                                      ;the number of coeffs is fixed
coeff:      .double 0.5, 1.0, 0.5
N_SAMPLES:  .word   5
sample:     .double 1.0, 2.0, 1.0, 2.0, 1.0
result:     .double 0.0, 0.0, 0.0, 0.0, 0.0
fp_one:     .double 1.0

    .text
main:
    lw     $t8, N_SAMPLES($zero)  ; t8 = n
    lw     $t9, N_COEFFS($zero)     
    l.d    f5, fp_one($zero)      ; f5 = 1.0  Moving this line up here decreases CPI by 0.03+0.13 
    
    sltu   $t9, $t9, $t8          ; if N_SAMPLES >= N_COEFFS, t9=1
    l.d    f1, coeff($zero)       ; f1 = coeff[0]， Move this up decrese the CPI by 0.013
    beqz   $t9, bail              ; if N_SAMPLES < N_COEFFS, stop   
	l.d    f2, coeff+8($zero)     ; f2 = coeff[1]
    
    c.lt.d f1, f0                 ; compare if f1<0
    bc1f   coeff1                 ; if coeff[0]>0, goto coeff1
    l.d    f3, coeff+16($zero)    ; f3 = coeff[2]  
    sub.d  f1, f0, f1             ; if coeff[0]<0, f1 =-coeff[0]
    
coeff1:
    c.lt.d f2, f0
    bc1f   coeff2                 ; if f2>0, goto coeff2
    l.d    f6, sample($zero)      ; sample[0], word is 32-bits, only even regnumber is allowed
    sub.d  f2, f0, f2             ; if f2<0, f2 =-f2
    
coeff2:
    c.lt.d f3, f0
    add.d  f4, f2, f1             ; f4 = f1+f2, Move this down here decrese the CPI by 0.013
    bc1f   coeff3                 ; if f3>0, goto coeff3
    daddi  $t1, $zero, 1          ; t1 = 1
    sub.d  f3, f0, f3             ; if f3<0, f3=-f3   
    
coeff3:
    add.d  f4, f4, f3             ; f4= coeff[0]+coeff[1]+coeff[2]
    dsub   $t7, $t8, $t1          ; t8 = n, t7 = n-1  
    daddi  $t2, $zero, 16
    div.d  f5, f5, f4             ; f5 = 1/norm
	
    s.d    f6, result($zero)      ; result[0] = sample[0]
    dsll   $t7, $t7, 3            ; t7 = (n-1)*8
    dsll   $t9, $t8, 3            ; t9 = n*8, a double is 8 bytes
    
    l.d    f8, sample($t7)        ; f8 = sample[n-1]
    
    daddi  $t0, $zero, 0
    daddi  $t1, $zero, 8
    
    s.d    f8, result($t7)        ; result[n-1] = sample[n-1]  
                                  ; f1 = |coeff[0]|, f2 = |coeff[1]|, f3 = |coeff[2]|, f5 = 1/norm
								  
    l.d   f6, sample($t1)         ; sample >= 3, the first 3 samples are always loaded, so do it outside the for-loop
    l.d   f8, sample($t2)         ; load/store -> 2 cycles
    l.d   f4, sample($t0)         ; finish loading the first 3 samples to f4, f6, f8, go into loop from below	
	
    mul.d  f1, f1, f5             ; f1 = coeff[0]/norm	
    mul.d  f2, f2, f5             ; f2 = coeff[1]/norm
    mul.d  f3, f3, f5             ; f3 = coeff[2]/norm

forloop: 	

    mul.d f10, f4, f1             ;  1st coeff * 1st sample
	mul.d f12, f6, f2             ;  2nd coeff * 2nd sample
	
    daddi $t2, $t2, 8
    daddi $t0, $t0, 8
	
    mul.d f16, f8, f3             ;  3rd coeff * 3rd sample, Multiple -> 4 cycles (seen from MIPS cycles window)
    
    mov.d f4, f6                  ; move f6 to f4,  prepare for the next loop 
    add.d f14, f10, f12           ; don't move this line. if moved to anywhere below, cycles increses to 955. now it's 907
    mov.d f6, f8                  ; move f8 to f6           ;  here finish calculating each element for one result (3 mul.d needed) 
    
    
    l.d   f8, sample($t2)         ; load or store needs 2 cycles, load the next sample 
    add.d f14, f14, f16           ; Do not move this line up, that will increase the cycles a lot ; final sum to be stored in result, one result finished completely
     	
    bne   $t9, $t2, forloop       ; t1 != t9 =(n-1)*8, do the loop, till store result[n-2]
    s.d   f14, result($t0)        ; Store result, from result[1], THIS LINE will be executed anyway. Moving s.d down to here decreases
                                 
bail:
	halt