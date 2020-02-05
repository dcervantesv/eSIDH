//*******************************************************************************************
// SIDH: an efficient supersingular isogeny cryptography library
//
// Abstract: field arithmetic in x64 assembly for P443 on Linux 
//*******************************************************************************************  
.intel_syntax noprefix 

// Registers that are used for parameter passing:
#define reg_p1  rdi
#define reg_p2  rsi
#define reg_p3  rdx

// p443 + 1
#define p443p1_2   0x0000000000000000
#define p443p1_3   0x1D236727C0000000
#define p443p1_4   0xD9EBCCB301480652
#define p443p1_5   0xEA6695D5C29E0A13
#define p443p1_6   0x048F5AB5AA19ECE9

#define p443_0     0xFFFFFFFFFFFFFFFF
#define p443_3     0x1D236727BFFFFFFF
#define p443_4     0xD9EBCCB301480652
#define p443_5     0xEA6695D5C29E0A13
#define p443_6     0x048F5AB5AA19ECE9

#define p443x2_0   0xFFFFFFFFFFFFFFFE
#define p443x2_1   0xFFFFFFFFFFFFFFFF
#define p443x2_3   0x3A46CE4F7FFFFFFF
#define p443x2_4   0xB3D7996602900CA4
#define p443x2_5   0xD4CD2BAB853C1427
#define p443x2_6   0x091EB56B5433D9D3   

.text
//***********************************************************************
//  Field addition
//  Operation: c [reg_p3] = a [reg_p1] + b [reg_p2]
//*********************************************************************** 
.global fpadd443_asm
.p2align 4
fpadd443_asm:
  push   r12
  push   r13
  push   r14
  push   r15
  
  mov   r8 , [reg_p1]
  mov   r9 , [reg_p1+8]
  mov   r10, [reg_p1+16]
  mov   r11, [reg_p1+24]
  mov   r12, [reg_p1+32]
  mov   r13, [reg_p1+40]
  mov   r14, [reg_p1+48]
  add   r8 , [reg_p2] 
  adc   r9 , [reg_p2+8] 
  adc   r10, [reg_p2+16] 
  adc   r11, [reg_p2+24] 
  adc   r12, [reg_p2+32] 
  adc   r13, [reg_p2+40] 
  adc   r14, [reg_p2+48] 

  mov   rax, p443x2_0
  sub   r8 , rax
  mov   rax, p443x2_1
  sbb   r9 , rax
  sbb   r10, rax
  mov   rax, p443x2_3
  sbb   r11, rax
  mov   rax, p443x2_4
  sbb   r12, rax
  mov   rax, p443x2_5
  sbb   r13, rax
  mov   rax, p443x2_6
  sbb   r14, rax
  
  mov   [reg_p3]   , r8
  mov   [reg_p3+8] , r9
  mov   [reg_p3+16], r10
  mov   [reg_p3+24], r11
  mov   [reg_p3+32], r12
  mov   [reg_p3+40], r13
  mov   [reg_p3+48], r14
  mov   rax, 0
  sbb   rax, 0
  
  mov   rsi, p443x2_0
  mov   r8 , p443x2_1
  mov   r10, p443x2_3
  mov   r11, p443x2_4
  mov   r12, p443x2_5
  mov   r13, p443x2_6
  and   rsi, rax
  and   r8 , rax
  and   r10, rax
  and   r11, rax
  and   r12, rax
  and   r13, rax
  and   r14, rax
  
  mov   rax, [reg_p3]
  add   rax, rsi  
  mov   [reg_p3], rax
  mov   rax, [reg_p3+8]
  adc   rax, r8 
  mov   [reg_p3+8], rax  
  mov   rax, [reg_p3+16]
  adc   rax, r8 
  mov   [reg_p3+16], rax  
  mov   rax, [reg_p3+24]  
  adc   rax, r10 
  mov   [reg_p3+24], rax 
  mov   rax, [reg_p3+32]  
  adc   rax, r11 
  mov   [reg_p3+32], rax 
  mov   rax, [reg_p3+40]    
  adc   rax, r12 
  mov   [reg_p3+40], rax 
  mov   rax, [reg_p3+48]   
  adc   rax, r13 
  mov   [reg_p3+48], rax 
  
  pop    r15
  pop    r14
  pop    r13
  pop    r12
  ret

//***********************************************************************
//  Field subtraction
//  Operation: c [reg_p3] = a [reg_p1] - b [reg_p2]
//*********************************************************************** 
.global fpsub443_asm
.p2align 4
fpsub443_asm:
  push   r12
  push   r13
  push   r14
  push   r15

  mov   r8 , [reg_p1]
  mov   r9 , [reg_p1+8]
  mov   r10, [reg_p1+16]
  mov   r11, [reg_p1+24]
  sub   r8 , [reg_p2] 
  sbb   r9 , [reg_p2+8] 
  sbb   r10, [reg_p2+16] 
  sbb   r11, [reg_p2+24] 
  mov   [reg_p3]  ,  r8
  mov   [reg_p3+8],  r9
  mov   [reg_p3+16], r10
  mov   [reg_p3+24], r11
  mov   r12, [reg_p1+32]
  mov   r13, [reg_p1+40]
  mov   r14, [reg_p1+48]
  sbb   r12, [reg_p2+32] 
  sbb   r13, [reg_p2+40] 
  sbb   r14, [reg_p2+48] 
  mov   [reg_p3+32], r12
  mov   [reg_p3+40], r13
  mov   [reg_p3+48], r14
  mov   r12, 0
  sbb   r12, 0

  mov   rsi, p443x2_0
  and   rsi, r12
  mov   r8 , p443x2_1
  and   r8 , r12
  mov   r9 , p443x2_3
  and   r9 , r12
  mov   r11, p443x2_4
  and   r11, r12
  mov   r10, p443x2_5
  and   r10, r12
  mov   r13, p443x2_6
  and   r13, r12
  
  mov   rax, [reg_p3]
  add   rax, rsi  
  mov   [reg_p3], rax
  mov   rax, [reg_p3+8]
  adc   rax, r8 
  mov   [reg_p3+8], rax  
  mov   rax, [reg_p3+16]
  adc   rax, r8 
  mov   [reg_p3+16], rax  
  mov   rax, [reg_p3+24]  
  adc   rax, r9 
  mov   [reg_p3+24], rax 
  mov   rax, [reg_p3+32]  
  adc   rax, r11 
  mov   [reg_p3+32], rax 
  mov   rax, [reg_p3+40]    
  adc   rax, r10 
  mov   [reg_p3+40], rax 
  mov   rax, [reg_p3+48]   
  adc   rax, r13
  mov   [reg_p3+48], rax 

  pop    r15
  pop    r14
  pop    r13
  pop    r12
  ret

//***********************************************************************
//  Integer multiplication
//  Based on Schoolbook method
//  Operation: c [reg_p3] = a [reg_p1] * b [reg_p2]
//  NOTE: a=c or b=c are not allowed
//***********************************************************************
.global mul443_asm
mul443_asm:    
    push   r12
    push   r13 
    push   r14 
    push   r15
	push   rbx
    mov   rcx, reg_p3 

	mov  rdx, [reg_p1]
	mulx  r9 , r8 , [reg_p2]
	mov  [rcx], r8
    mulx  r10, rax, [reg_p2+8]
	add  r9 , rax
	mulx  r11, rax, [reg_p2+16]
	adc  r10, rax
	mulx  r12, rax, [reg_p2+24]
	adc  r11, rax
	mulx  r13, rax, [reg_p2+32]
	adc  r12, rax
	mulx  r14, rax, [reg_p2+40]
	adc  r13, rax
	mulx  r15, rax, [reg_p2+48]
	adc  r14, rax
	adc  r15, 0

    xor  r8 , r8
    mov  rdx, [reg_p1+8]
    mulx  rbx, rax, [reg_p2]
    adcx  r9 , rax
    adox  r10, rbx
	mov  [rcx+8], r9
    mulx  rbx, rax, [reg_p2+8]
    adcx  r10, rax
    adox  r11, rbx
    mulx  rbx, rax, [reg_p2+16]
    adcx  r11, rax
    adox  r12, rbx
    mulx  rbx, rax, [reg_p2+24]
    adcx  r12, rax
    adox  r13, rbx
    mulx  rbx, rax, [reg_p2+32]
    adcx  r13, rax
    adox  r14, rbx
    mulx  rbx, rax, [reg_p2+40]
    adcx  r14, rax
    adox  r15, rbx
    mulx  rbx, rax, [reg_p2+48]
    adcx  r15, rax
    adox  r8 , rbx
    adc  r8 , 0

    xor  r9 , r9
    mov  rdx, [reg_p1+16]
    mulx  rbx, rax, [reg_p2]
    adcx  r10, rax
    adox  r11, rbx
	mov  [rcx+16], r10
    mulx  rbx, rax, [reg_p2+8]
    adcx  r11, rax
    adox  r12, rbx
    mulx  rbx, rax, [reg_p2+16]
    adcx  r12, rax
    adox  r13, rbx
    mulx  rbx, rax, [reg_p2+24]
    adcx  r13, rax
    adox  r14, rbx
    mulx  rbx, rax, [reg_p2+32]
    adcx  r14, rax
    adox  r15, rbx
    mulx  rbx, rax, [reg_p2+40]
    adcx  r15, rax
    adox  r8 , rbx
    mulx  rbx, rax, [reg_p2+48]
    adcx  r8 , rax
    adox  r9 , rbx
    adc  r9 , 0

    xor  r10, r10
    mov  rdx, [reg_p1+24]
    mulx  rbx, rax, [reg_p2]
    adcx  r11, rax
    adox  r12, rbx
	mov  [rcx+24], r11
    mulx  rbx, rax, [reg_p2+8]
    adcx  r12, rax
    adox  r13, rbx
    mulx  rbx, rax, [reg_p2+16]
    adcx  r13, rax
    adox  r14, rbx
    mulx  rbx, rax, [reg_p2+24]
    adcx  r14, rax
    adox  r15, rbx
    mulx  rbx, rax, [reg_p2+32]
    adcx  r15, rax
    adox  r8 , rbx
    mulx  rbx, rax, [reg_p2+40]
    adcx  r8 , rax
    adox  r9 , rbx
    mulx  rbx, rax, [reg_p2+48]
    adcx  r9 , rax
    adox  r10, rbx
    adc  r10, 0

    xor  r11, r11
    mov  rdx, [reg_p1+32]
    mulx  rbx, rax, [reg_p2]
    adcx  r12, rax
    adox  r13, rbx
	mov  [rcx+32], r12
    mulx  rbx, rax, [reg_p2+8]
    adcx  r13, rax
    adox  r14, rbx
    mulx  rbx, rax, [reg_p2+16]
    adcx  r14, rax
    adox  r15, rbx
    mulx  rbx, rax, [reg_p2+24]
    adcx  r15, rax
    adox  r8 , rbx
    mulx  rbx, rax, [reg_p2+32]
    adcx  r8 , rax
    adox  r9 , rbx
    mulx  rbx, rax, [reg_p2+40]
    adcx  r9 , rax
    adox  r10, rbx
    mulx  rbx, rax, [reg_p2+48]
    adcx  r10, rax
    adox  r11, rbx
    adc  r11, 0

    xor  r12, r12
    mov  rdx, [reg_p1+40]
    mulx  rbx, rax, [reg_p2]
    adcx  r13, rax
    adox  r14, rbx
	mov  [rcx+40], r13
    mulx  rbx, rax, [reg_p2+8]
    adcx  r14, rax
    adox  r15, rbx
    mulx  rbx, rax, [reg_p2+16]
    adcx  r15, rax
    adox  r8 , rbx
    mulx  rbx, rax, [reg_p2+24]
    adcx  r8 , rax
    adox  r9 , rbx
    mulx  rbx, rax, [reg_p2+32]
    adcx  r9 , rax
    adox  r10, rbx
    mulx  rbx, rax, [reg_p2+40]
    adcx  r10, rax
    adox  r11, rbx
    mulx  rbx, rax, [reg_p2+48]
    adcx  r11, rax
    adox  r12, rbx
    adc  r12, 0

    xor  r13, r13
    mov  rdx, [reg_p1+48]
    mulx  rbx, rax, [reg_p2]
    adcx  r14, rax
    adox  r15, rbx
	mov  [rcx+48], r14
    mulx  rbx, rax, [reg_p2+8]
    adcx  r15, rax
    adox  r8 , rbx
    mulx  rbx, rax, [reg_p2+16]
    adcx  r8 , rax
    adox  r9 , rbx
    mulx  rbx, rax, [reg_p2+24]
    adcx  r9 , rax
    adox  r10, rbx
    mulx  rbx, rax, [reg_p2+32]
    adcx  r10, rax
    adox  r11, rbx
    mulx  rbx, rax, [reg_p2+40]
    adcx  r11, rax
    adox  r12, rbx
    mulx  rbx, rax, [reg_p2+48]
    adcx  r12, rax
    adox  r13, rbx
    adc  r13, 0

	mov  [rcx+56], r15
	mov  [rcx+64], r8
	mov  [rcx+72], r9
	mov  [rcx+80], r10
	mov  [rcx+88], r11
	mov  [rcx+96], r12
	mov  [rcx+104], r13

    pop    rbx
    pop    r15
    pop    r14
    pop    r13
    pop    r12
    ret

//***********************************************************************
//  Integer squaring
//  Based on Schoolbook method
//  Operation: c [reg_p3] = a [reg_p1] * b [reg_p2]
//  NOTE: a=c or b=c are not allowed
//***********************************************************************
.global sqr443_asm
sqr443_asm:    
    push   r12
    push   r13 
    push   r14 
    push   r15
	push   rbx
    mov   rcx, reg_p2 

	mov  rdx, [reg_p1]
	mulx  r9 , r8 , [reg_p1]
	mov  [rcx], r8
    mulx  r10, rax, [reg_p1+8]
	add  r9 , rax
	mulx  r11, rax, [reg_p1+16]
	adc  r10, rax
	mulx  r12, rax, [reg_p1+24]
	adc  r11, rax
	mulx  r13, rax, [reg_p1+32]
	adc  r12, rax
	mulx  r14, rax, [reg_p1+40]
	adc  r13, rax
	mulx  r15, rax, [reg_p1+48]
	adc  r14, rax
	adc  r15, 0

    xor  r8 , r8
    mov  rdx, [reg_p1+8]
    mulx  rbx, rax, [reg_p1]
    adcx  r9 , rax
    adox  r10, rbx
	mov  [rcx+8], r9
    mulx  rbx, rax, [reg_p1+8]
    adcx  r10, rax
    adox  r11, rbx
    mulx  rbx, rax, [reg_p1+16]
    adcx  r11, rax
    adox  r12, rbx
    mulx  rbx, rax, [reg_p1+24]
    adcx  r12, rax
    adox  r13, rbx
    mulx  rbx, rax, [reg_p1+32]
    adcx  r13, rax
    adox  r14, rbx
    mulx  rbx, rax, [reg_p1+40]
    adcx  r14, rax
    adox  r15, rbx
    mulx  rbx, rax, [reg_p1+48]
    adcx  r15, rax
    adox  r8 , rbx
    adc  r8 , 0

    xor  r9 , r9
    mov  rdx, [reg_p1+16]
    mulx  rbx, rax, [reg_p1]
    adcx  r10, rax
    adox  r11, rbx
	mov  [rcx+16], r10
    mulx  rbx, rax, [reg_p1+8]
    adcx  r11, rax
    adox  r12, rbx
    mulx  rbx, rax, [reg_p1+16]
    adcx  r12, rax
    adox  r13, rbx
    mulx  rbx, rax, [reg_p1+24]
    adcx  r13, rax
    adox  r14, rbx
    mulx  rbx, rax, [reg_p1+32]
    adcx  r14, rax
    adox  r15, rbx
    mulx  rbx, rax, [reg_p1+40]
    adcx  r15, rax
    adox  r8 , rbx
    mulx  rbx, rax, [reg_p1+48]
    adcx  r8 , rax
    adox  r9 , rbx
    adc  r9 , 0

    xor  r10, r10
    mov  rdx, [reg_p1+24]
    mulx  rbx, rax, [reg_p1]
    adcx  r11, rax
    adox  r12, rbx
	mov  [rcx+24], r11
    mulx  rbx, rax, [reg_p1+8]
    adcx  r12, rax
    adox  r13, rbx
    mulx  rbx, rax, [reg_p1+16]
    adcx  r13, rax
    adox  r14, rbx
    mulx  rbx, rax, [reg_p1+24]
    adcx  r14, rax
    adox  r15, rbx
    mulx  rbx, rax, [reg_p1+32]
    adcx  r15, rax
    adox  r8 , rbx
    mulx  rbx, rax, [reg_p1+40]
    adcx  r8 , rax
    adox  r9 , rbx
    mulx  rbx, rax, [reg_p1+48]
    adcx  r9 , rax
    adox  r10, rbx
    adc  r10, 0

    xor  r11, r11
    mov  rdx, [reg_p1+32]
    mulx  rbx, rax, [reg_p1]
    adcx  r12, rax
    adox  r13, rbx
	mov  [rcx+32], r12
    mulx  rbx, rax, [reg_p1+8]
    adcx  r13, rax
    adox  r14, rbx
    mulx  rbx, rax, [reg_p1+16]
    adcx  r14, rax
    adox  r15, rbx
    mulx  rbx, rax, [reg_p1+24]
    adcx  r15, rax
    adox  r8 , rbx
    mulx  rbx, rax, [reg_p1+32]
    adcx  r8 , rax
    adox  r9 , rbx
    mulx  rbx, rax, [reg_p1+40]
    adcx  r9 , rax
    adox  r10, rbx
    mulx  rbx, rax, [reg_p1+48]
    adcx  r10, rax
    adox  r11, rbx
    adc  r11, 0

    xor  r12, r12
    mov  rdx, [reg_p1+40]
    mulx  rbx, rax, [reg_p1]
    adcx  r13, rax
    adox  r14, rbx
	mov  [rcx+40], r13
    mulx  rbx, rax, [reg_p1+8]
    adcx  r14, rax
    adox  r15, rbx
    mulx  rbx, rax, [reg_p1+16]
    adcx  r15, rax
    adox  r8 , rbx
    mulx  rbx, rax, [reg_p1+24]
    adcx  r8 , rax
    adox  r9 , rbx
    mulx  rbx, rax, [reg_p1+32]
    adcx  r9 , rax
    adox  r10, rbx
    mulx  rbx, rax, [reg_p1+40]
    adcx  r10, rax
    adox  r11, rbx
    mulx  rbx, rax, [reg_p1+48]
    adcx  r11, rax
    adox  r12, rbx
    adc  r12, 0

    xor  r13, r13
    mov  rdx, [reg_p1+48]
    mulx  rbx, rax, [reg_p1]
    adcx  r14, rax
    adox  r15, rbx
	mov  [rcx+48], r14
    mulx  rbx, rax, [reg_p1+8]
    adcx  r15, rax
    adox  r8 , rbx
    mulx  rbx, rax, [reg_p1+16]
    adcx  r8 , rax
    adox  r9 , rbx
    mulx  rbx, rax, [reg_p1+24]
    adcx  r9 , rax
    adox  r10, rbx
    mulx  rbx, rax, [reg_p1+32]
    adcx  r10, rax
    adox  r11, rbx
    mulx  rbx, rax, [reg_p1+40]
    adcx  r11, rax
    adox  r12, rbx
    mulx  rbx, rax, [reg_p1+48]
    adcx  r12, rax
    adox  r13, rbx
    adc  r13, 0

	mov  [rcx+56], r15
	mov  [rcx+64], r8
	mov  [rcx+72], r9
	mov  [rcx+80], r10
	mov  [rcx+88], r11
	mov  [rcx+96], r12
	mov  [rcx+104], r13

    pop    rbx
    pop    r15
    pop    r14
    pop    r13
    pop    r12
    ret

//***********************************************************************
//  Montgomery reduction
//  Based on comba method
//  Operation: c [reg_p2] = a [reg_p1]
//  NOTE: a=c is not allowed
//*********************************************************************** 
.global rdc443_asm
.p2align 4
rdc443_asm:
  push   r12
  push   r13 
  push   r14 
  push   r15

  mov  rdx, p443p1_3
  mulx  r9 , r8 , [reg_p1]
  mulx  r10, rax, [reg_p1+8]
  add  r9 , rax
  mulx  r11, rax, [reg_p1+16]
  adc  r10, rax
  adc  r11, 0

  xor  r12, r12
  mov  rdx, p443p1_4
  mulx  rcx, rax, [reg_p1]
  adcx  r9 , rax
  adox  r10, rcx
  mulx  rcx, rax, [reg_p1+8]
  adcx  r10, rax
  adox  r11, rcx
  mulx  rcx, rax, [reg_p1+16]
  adcx  r11, rax
  adox  r12, rcx
  adc  r12, 0

  xor  r13, r13
  mov  rdx, p443p1_5
  mulx  rcx, rax, [reg_p1]
  adcx  r10, rax
  adox  r11, rcx
  mulx  rcx, rax, [reg_p1+8]
  adcx  r11, rax
  adox  r12, rcx
  mulx  rcx, rax, [reg_p1+16]
  adcx  r12, rax
  adox  r13, rcx
  adc  r13, 0

  xor  r14, r14
  mov  rdx, p443p1_6
  mulx  rcx, rax, [reg_p1]
  adcx  r11, rax
  adox  r12, rcx
  mulx  rcx, rax, [reg_p1+8]
  adcx  r12, rax
  adox  r13, rcx
  mulx  rcx, rax, [reg_p1+16]
  adcx  r13, rax
  adox  r14, rcx
  adc  r14, 0

  mov  rax, 0
  add  [reg_p1+24], r8
  adc  [reg_p1+32], r9
  adc  [reg_p1+40], r10
  adc  [reg_p1+48], r11
  adc  [reg_p1+56], r12
  adc  [reg_p1+64], r13
  adc  [reg_p1+72], r14
  adc  [reg_p1+80], rax

  mov   rdx, p443p1_3
  mulx  r9 , r8 , [reg_p1+24]
  mulx  r10, rax, [reg_p1+32]
  add  r9 , rax
  mulx  r11, rax, [reg_p1+40]
  adc  r10, rax
  adc  r11, 0

  xor  r12, r12
  mov  rdx, p443p1_4
  mulx  rcx, rax, [reg_p1+24]
  adcx  r9 , rax
  adox  r10, rcx
  mulx  rcx, rax, [reg_p1+32]
  adcx  r10, rax
  adox  r11, rcx
  mulx  rcx, rax, [reg_p1+40]
  adcx  r11, rax
  adox  r12, rcx
  adc  r12, 0

  xor  r13, r13
  mov  rdx, p443p1_5
  mulx  rcx, rax, [reg_p1+24]
  adcx  r10, rax
  adox  r11, rcx
  mulx  rcx, rax, [reg_p1+32]
  adcx  r11, rax
  adox  r12, rcx
  mulx  rcx, rax, [reg_p1+40]
  adcx  r12, rax
  adox  r13, rcx
  adc  r13, 0

  xor  r14, r14
  mov  rdx, p443p1_6
  mulx  rcx, rax, [reg_p1+24]
  adcx  r11, rax
  adox  r12, rcx
  mulx  rcx, rax, [reg_p1+32]
  adcx  r12, rax
  adox  r13, rcx
  mulx  rcx, rax, [reg_p1+40]
  adcx  r13, rax
  adox  r14, rcx
  adc  r14, 0

  mov  rax, 0
  add  [reg_p1+48], r8
  adc  [reg_p1+56], r9
  adc  [reg_p1+64], r10
  adc  [reg_p1+72], r11
  adc  [reg_p1+80], r12
  adc  [reg_p1+88], r13
  adc  [reg_p1+96], r14
  adc  [reg_p1+104], rax

  mov  rdx, [reg_p1+48]
  mov  r8 , p443p1_3
  mov  r9 , p443p1_4
  mov  r10, p443p1_5
  mov  rcx, p443p1_6

  mulx  r12, r11, r8
  mulx  r13, rax, r9
  add  r12, rax
  mulx  r14, rax, r10
  adc  r13, rax
  mulx  r15, rax, rcx
  adc  r14, rax
  adc  r15, 0

  mov  r8 , [reg_p1+56]
  mov  r9 , [reg_p1+64]
  add  r11, [reg_p1+72]
  adc  r12, [reg_p1+80]
  adc  r13, [reg_p1+88]
  adc  r14, [reg_p1+96]
  adc  r15, [reg_p1+104]

  mov  [reg_p2]   , r8
  mov  [reg_p2+8] , r9
  mov  [reg_p2+16], r11
  mov  [reg_p2+24], r12
  mov  [reg_p2+32], r13
  mov  [reg_p2+40], r14
  mov  [reg_p2+48], r15

  pop    r15
  pop    r14
  pop    r13
  pop    r12
  ret


  
//***********************************************************************
//  443-bit multiprecision addition
//  Operation: c [reg_p3] = a [reg_p1] + b [reg_p2]
//*********************************************************************** 
.global mp_add443_asm
mp_add443_asm:
  push   r12
  push   r13
  push   r14
  push   r15
  
  mov    r8, [reg_p1]
  mov    r9, [reg_p1+8]
  mov    r10, [reg_p1+16]
  mov    r11, [reg_p1+24]
  mov    r12, [reg_p1+32]
  mov    r13, [reg_p1+40]
  mov    r14, [reg_p1+48]

  add    r8, [reg_p2] 
  adc    r9, [reg_p2+8] 
  adc    r10, [reg_p2+16] 
  adc    r11, [reg_p2+24] 
  adc    r12, [reg_p2+32] 
  adc    r13, [reg_p2+40] 
  adc    r14, [reg_p2+48] 

  mov    [reg_p3], r8
  mov    [reg_p3+8], r9
  mov    [reg_p3+16], r10
  mov    [reg_p3+24], r11
  mov    [reg_p3+32], r12
  mov    [reg_p3+40], r13
  mov    [reg_p3+48], r14
  
  pop    r15
  pop    r14
  pop    r13
  pop    r12
  ret


//***********************************************************************
//  2x443-bit multiprecision addition
//  Operation: c [reg_p3] = a [reg_p1] + b [reg_p2]
//*********************************************************************** 
.global mp_add443x2_asm
mp_add443x2_asm:
  push   r12
  push   r13
  push   r14
  push   r15
  
  mov    r8, [reg_p1]
  mov    r9, [reg_p1+8]
  mov    r10, [reg_p1+16]
  mov    r11, [reg_p1+24]
  mov    r12, [reg_p1+32]
  mov    r13, [reg_p1+40]
  mov    r14, [reg_p1+48]

  add    r8, [reg_p2] 
  adc    r9, [reg_p2+8] 
  adc    r10, [reg_p2+16] 
  adc    r11, [reg_p2+24] 
  adc    r12, [reg_p2+32] 
  adc    r13, [reg_p2+40] 
  adc    r14, [reg_p2+48] 

  mov    [reg_p3], r8
  mov    [reg_p3+8], r9
  mov    [reg_p3+16], r10
  mov    [reg_p3+24], r11
  mov    [reg_p3+32], r12
  mov    [reg_p3+40], r13
  mov    [reg_p3+48], r14
  
  mov    r8, [reg_p1+56]
  mov    r9, [reg_p1+64]
  mov    r10, [reg_p1+72]
  mov    r11, [reg_p1+80]
  mov    r12, [reg_p1+88]
  mov    r13, [reg_p1+96]
  mov    r14, [reg_p1+104]

  adc    r8, [reg_p2+56] 
  adc    r9, [reg_p2+64] 
  adc    r10, [reg_p2+72] 
  adc    r11, [reg_p2+80] 
  adc    r12, [reg_p2+88] 
  adc    r13, [reg_p2+96] 
  adc    r14, [reg_p2+104] 

  mov    [reg_p3+56], r8
  mov    [reg_p3+64], r9
  mov    [reg_p3+72], r10
  mov    [reg_p3+80], r11
  mov    [reg_p3+88], r12
  mov    [reg_p3+96], r13
  mov    [reg_p3+104], r14
  
  pop    r15
  pop    r14
  pop    r13
  pop    r12
  ret
  
//***********************************************************************
//  2x443-bit multiprecision subtraction
//  Operation: c [reg_p3] = a [reg_p1] - b [reg_p2]. Returns borrow mask
//*********************************************************************** 
.global mp_sub443x2_asm
mp_sub443x2_asm:
  xor    rax, rax
  mov    r8, [reg_p1]
  mov    r9, [reg_p1+8]
  mov    r10, [reg_p1+16]
  mov    r11, [reg_p1+24]
  mov    rcx, [reg_p1+32]
  sub    r8, [reg_p2] 
  sbb    r9, [reg_p2+8] 
  sbb    r10, [reg_p2+16] 
  sbb    r11, [reg_p2+24] 
  sbb    rcx, [reg_p2+32] 
  mov    [reg_p3], r8
  mov    [reg_p3+8], r9
  mov    [reg_p3+16], r10
  mov    [reg_p3+24], r11
  mov    [reg_p3+32], rcx

  mov    r8, [reg_p1+40]
  mov    r9, [reg_p1+48]
  mov    r10, [reg_p1+56] 
  mov    r11, [reg_p1+64]
  mov    rcx, [reg_p1+72] 
  sbb    r8, [reg_p2+40] 
  sbb    r9, [reg_p2+48] 
  sbb    r10, [reg_p2+56]
  sbb    r11, [reg_p2+64] 
  sbb    rcx, [reg_p2+72]
  mov    [reg_p3+40], r8
  mov    [reg_p3+48], r9
  mov    [reg_p3+56], r10
  mov    [reg_p3+64], r11
  mov    [reg_p3+72], rcx
  
  mov    r8, [reg_p1+80]
  mov    r9, [reg_p1+88] 
  mov    r10, [reg_p1+96]
  mov    r11, [reg_p1+104]
  sbb    r8, [reg_p2+80]
  sbb    r9, [reg_p2+88]
  sbb    r10, [reg_p2+96] 
  sbb    r11, [reg_p2+104] 
  sbb    rax, 0
  mov    [reg_p3+80], r8 
  mov    [reg_p3+88], r9
  mov    [reg_p3+96], r10
  mov    [reg_p3+104], r11
  
  ret  
  
//***********************************************************************
//  Double 2x443-bit multiprecision subtraction
//  Operation: c [reg_p3] = c [reg_p3] - a [reg_p1] - b [reg_p2]
//*********************************************************************** 
.global mp_dblsub443x2_asm
mp_dblsub443x2_asm:
  push   r12
  push   r13
  push   r14
  push   r15
  
  xor    rax, rax
  mov    r8, [reg_p3]
  mov    r9, [reg_p3+8]
  mov    r10, [reg_p3+16]
  mov    r11, [reg_p3+24]
  mov    r12, [reg_p3+32]
  mov    r13, [reg_p3+40]
  mov    r14, [reg_p3+48]
  mov    r15, [reg_p3+56]
  sub    r8, [reg_p1]
  sbb    r9, [reg_p1+8] 
  sbb    r10, [reg_p1+16] 
  sbb    r11, [reg_p1+24] 
  sbb    r12, [reg_p1+32] 
  sbb    r13, [reg_p1+40] 
  sbb    r14, [reg_p1+48] 
  sbb    r15, [reg_p1+56]
  adc    rax, 0
  sub    r8, [reg_p2]
  sbb    r9, [reg_p2+8] 
  sbb    r10, [reg_p2+16] 
  sbb    r11, [reg_p2+24] 
  sbb    r12, [reg_p2+32] 
  sbb    r13, [reg_p2+40] 
  sbb    r14, [reg_p2+48] 
  sbb    r15, [reg_p2+56]
  adc    rax, 0  
  mov    [reg_p3], r8
  mov    [reg_p3+8], r9
  mov    [reg_p3+16], r10
  mov    [reg_p3+24], r11
  mov    [reg_p3+32], r12
  mov    [reg_p3+40], r13
  mov    [reg_p3+48], r14
  mov    [reg_p3+56], r15
    
  mov    r8, [reg_p3+64]
  mov    r9, [reg_p3+72]
  mov    r10, [reg_p3+80]
  mov    r11, [reg_p3+88]
  mov    r12, [reg_p3+96]
  mov    r13, [reg_p3+104]
  mov    r14, [reg_p3+112]
  mov    r15, [reg_p3+120]
  sub    r8, rax 
  sbb    r8, [reg_p1+64] 
  sbb    r9, [reg_p1+72] 
  sbb    r10, [reg_p1+80] 
  sbb    r11, [reg_p1+88] 
  sbb    r12, [reg_p1+96] 
  sbb    r13, [reg_p1+104] 
  sub    r8, [reg_p2+64] 
  sbb    r9, [reg_p2+72] 
  sbb    r10, [reg_p2+80] 
  sbb    r11, [reg_p2+88] 
  sbb    r12, [reg_p2+96] 
  sbb    r13, [reg_p2+104] 
  mov    [reg_p3+64], r8
  mov    [reg_p3+72], r9
  mov    [reg_p3+80], r10
  mov    [reg_p3+88], r11
  mov    [reg_p3+96], r12
  mov    [reg_p3+104], r13
  
  pop    r15
  pop    r14
  pop    r13
  pop    r12
  ret
  
  
