    .data
N_COEFFS:   .word   3            ;the number of coeffs is fixed
coeff:      .double 0.5, 1.0, 0.5
N_SAMPLES:  .word   50
sample:     .double 2,3,4,1,1,1,1,1,2,1,1,2,1,2.0,1.0,1.0,2.0,1.0,2.0,1.0,1.0,2.0,1.0,2.0,1.0,1.0,2.0,1.0,2.0,1.0,1.0,2.0,1.0,2.0,1.0,1.0,2.0,1.0,2.0,1.0,1.0,2.0,1.0,2.0,1.0,1.0,2.0,1.0,2.0,1.0
result:     .double 0.0, 0.0, 0.0, 0.0, 0.0
fp_one:     .double 1.0

;=====================================================================================================
; Function          Latency      Initiation_Interval         total_cycles_needed
; Integer ALU         0                 1                          1
; Data Memory         1                 1                          2
; FP adder            3                 1                          4
; Multiplier          5                 1                          6
; Divider             24                25                     25, not pipelined
;======================================================================================================


; MIPS Cycles Window:
;An instruction that causes a stall is highlighted in BLUE
;instructions held up as a result of a stall are GRAY.


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
    nop    
    sub.d  f1, f0, f1             ; if coeff[0]<0, f1 =-coeff[0]
    
coeff1:
    c.lt.d f2, f0
    bc1f   coeff2                 ; if f2>0, goto coeff2
    l.d    f3, coeff+16($zero)    ; f3 = coeff[2]
    sub.d  f2, f0, f2             ; if f2<0, f2 =-f2
    
coeff2:
    c.lt.d f3, f0
    add.d  f4, f2, f1             ; f4 = f1+f2, Move this down here decrese the CPI by 0.013
    bc1f   coeff3                 ; if f3>0, goto coeff3
    l.d    f6, sample($zero)      ; sample[0], word is 32-bits, only even regnumber is allowed
    sub.d  f3, f0, f3             ; if f3<0, f3=-f3   
    
coeff3:
    add.d  f4, f4, f3             ; f4= coeff[0]+coeff[1]+coeff[2]
    div.d  f4, f5, f4             ; f4 = 1/norm
    daddi  $t1, $zero, 1          ; t1 = 1
    s.d    f6, result($zero)      ; result[0] = sample[0]
    
    dsub   $t9, $t8, $t1          ; t8 = n, t9 = n-1
    daddi  $t2, $zero, 16
    dsll   $t9, $t9, 3            ; t9 = (n-1)*8, a double is 8 bytes
    l.d    f8, sample($t9)        ; f8 = sample[n-1]
    daddi  $t0, $zero, 0
    s.d    f8, result($t9)        ; result[n-1] = sample[n-1]  

    daddi  $t1, $zero, 8
                                  ; f1 = |coeff[0]|, f2 = |coeff[1]|, 
                                  ; f3 = |coeff[2]|, f4 = 1/norm

	mul.d  f1, f1, f4             ; f1 = coeff[0]/norm	
    mul.d  f2, f2, f4             ; f2 = coeff[1]/norm
    mul.d  f3, f3, f4             ; f3 = coeff[2]/norm

forloop:
    l.d   f8, sample($t2)
    l.d   f4, sample($t0)     
	l.d   f6, sample($t1)
    
	mul.d f16, f8, f3             ;  3rd coeff * 3rd sample
    mul.d f10, f4, f1             ;  1st coeff * 1st sample
	mul.d f12, f6, f2             ;  2nd coeff * 2nd sample
	
    daddi $t0, $t0, 8        
	add.d f14, f10,f12
    add.d f4,  f14,f16           ; final sum to be stored in result
    daddi $t1, $t1, 8       
	daddi $t2, $t2, 8     
	
    s.d   f4, result($t0)         ; Store result, from result[1]  	

    bne   $t9, $t1, forloop       ; t1 != t9 =(n-1)*8, do the loop, till store result[n-2]
     
    ; load all the results to the register f18-f26 to check if they are correct.     
    ;l.d   f18, result($zero)    
    ;l.d   f20, result+8($zero)
    ;l.d   f22, result+16($zero)
    ;l.d   f24, result+24($zero)
    ;l.d   f26, result+32($zero)
    
bail:
    nop
	halt