//*******************************************************************************************
// SIDH: an efficient supersingular isogeny cryptography library
//
// Abstract: field arithmetic in x64 assembly for P765 on Linux 
//*******************************************************************************************  

.intel_syntax noprefix 

// Registers that are used for parameter passing:
#define reg_p1  rdi
#define reg_p2  rsi
#define reg_p3  rdx

// p765 + 1
#define p765p1_5   0x0000000000000000
#define p765p1_6   0xAFA37E4ADE8AD370
#define p765p1_7   0x295AB216AE6E7F61
#define p765p1_8   0x4CFB7AF5840E37DB
#define p765p1_9   0xC62A44C95DC67031
#define p765p1_10  0xBCBF3BE3A0B668BC
#define p765p1_11  0x19C0F8DEC9ABD1D5

#define p765_0     0xFFFFFFFFFFFFFFFF
#define p765_6     0xAFA37E4ADE8AD36F
#define p765_7     0x295AB216AE6E7F61
#define p765_8     0x4CFB7AF5840E37DB
#define p765_9     0xC62A44C95DC67031
#define p765_10    0xBCBF3BE3A0B668BC
#define p765_11    0x19C0F8DEC9ABD1D5

#define p765x2_0   0xFFFFFFFFFFFFFFFE
#define p765x2_1   0xFFFFFFFFFFFFFFFF
#define p765x2_6   0x5F46FC95BD15A6DF
#define p765x2_7   0x52B5642D5CDCFEC3
#define p765x2_8   0x99F6F5EB081C6FB6
#define p765x2_9   0x8C548992BB8CE062
#define p765x2_10  0x797E77C7416CD179
#define p765x2_11  0x3381F1BD9357A3AB

.text
//***********************************************************************
//  Field addition
//  Operation: c [reg_p3] = a [reg_p1] + b [reg_p2]
//*********************************************************************** 
.global fpadd765_asm
.p2align 4
fpadd765_asm:
  push   r12
  push   r13
  push   r14
  push   r15
  push   rbx
  
  mov   r8 , [reg_p1]
  mov   r9 , [reg_p1+8]
  mov   r10, [reg_p1+16]
  mov   r11, [reg_p1+24]
  add   r8 , [reg_p2] 
  adc   r9 , [reg_p2+8] 
  adc   r10, [reg_p2+16] 
  adc   r11, [reg_p2+24] 
  mov   r12, [reg_p1+32]
  mov   r13, [reg_p1+40]
  mov   r14, [reg_p1+48]
  mov   r15, [reg_p1+56] 
  adc   r12, [reg_p2+32] 
  adc   r13, [reg_p2+40] 
  adc   r14, [reg_p2+48] 
  adc   r15, [reg_p2+56]
  mov   rcx, [reg_p1+64]
  mov   rax, [reg_p1+72]
  adc   rcx, [reg_p2+64] 
  adc   rax, [reg_p2+72] 
  mov   [reg_p3+72], rax
  mov   rax, [reg_p1+80]
  adc   rax, [reg_p2+80] 
  mov   [reg_p3+80], rax
  mov   rax, [reg_p1+88]
  adc   rax, [reg_p2+88] 
  mov   [reg_p3+88], rax

  mov   rax, p765x2_0
  sub   r8, rax
  mov   rax, p765x2_1
  sbb   r9, rax
  sbb   r10, rax
  sbb   r11, rax
  sbb   r12, rax
  sbb   r13, rax
  mov   rax, p765x2_6
  sbb   r14, rax
  mov   rax, p765x2_7
  sbb   r15, rax
  mov   rax, p765x2_8
  sbb   rcx, rax
  mov   [reg_p3], r8
  mov   [reg_p3+8], r9
  mov   [reg_p3+16], r10
  mov   [reg_p3+24], r11
  mov   [reg_p3+32], r12
  mov   [reg_p3+40], r13
  mov   [reg_p3+48], r14
  mov   [reg_p3+56], r15
  mov   [reg_p3+64], rcx
  mov   r8, [reg_p3+72]
  mov   r9, [reg_p3+80]
  mov   r10, [reg_p3+88]
  mov   rax, p765x2_9
  sbb   r8, rax
  mov   rax, p765x2_10
  sbb   r9, rax
  mov   rax, p765x2_11
  sbb   r10, rax
  mov   [reg_p3+72], r8
  mov   [reg_p3+80], r9
  mov   [reg_p3+88], r10
  mov   rax, 0
  sbb   rax, 0
  
  mov   rbx, p765x2_0
  mov   r8, p765x2_1
  mov   r10, p765x2_6
  mov   r11, p765x2_7
  mov   r12, p765x2_8
  mov   r13, p765x2_9
  mov   r14, p765x2_10
  mov   r15, p765x2_11
  and   rbx, rax
  and   r8 , rax
  and   r10, rax
  and   r11, rax
  and   r12, rax
  and   r13, rax
  and   r14, rax
  and   r15, rax
  
  mov   rax, [reg_p3]
  add   rax, rbx  
  mov   [reg_p3], rax
  mov   rax, [reg_p3+8]
  adc   rax, r8 
  mov   [reg_p3+8], rax  
  mov   rax, [reg_p3+16]
  adc   rax, r8 
  mov   [reg_p3+16], rax  
  mov   rax, [reg_p3+24]  
  adc   rax, r8 
  mov   [reg_p3+24], rax 
  mov   rax, [reg_p3+32]  
  adc   rax, r8 
  mov   [reg_p3+32], rax 
  mov   rax, [reg_p3+40]    
  adc   rax, r8 
  mov   [reg_p3+40], rax 
  mov   rax, [reg_p3+48]   
  adc   rax, r10 
  mov   [reg_p3+48], rax 
  mov   rax, [reg_p3+56]   
  adc   rax, r11  
  mov   [reg_p3+56], rax 
  mov   rax, [reg_p3+64]  
  adc   rax, r12 
  mov   [reg_p3+64], rax 
  mov   rax, [reg_p3+72]   
  adc   rax, r13 
  mov   [reg_p3+72], rax 
  mov   rax, [reg_p3+80]   
  adc   rax, r14 
  mov   [reg_p3+80], rax 
  mov   rax, [reg_p3+88]   
  adc   rax, r15
  mov   [reg_p3+88], rax 
  
  pop    rbx
  pop    r15
  pop    r14
  pop    r13
  pop    r12
  ret


//***********************************************************************
//  Field subtraction
//  Operation: c [reg_p3] = a [reg_p1] - b [reg_p2]
//*********************************************************************** 
.global fpsub765_asm
.p2align 4
fpsub765_asm:
  push   r12
  push   r13
  push   r14
  push   r15
  push   rbx

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
  mov   r15, [reg_p1+56] 
  sbb   r12, [reg_p2+32] 
  sbb   r13, [reg_p2+40] 
  sbb   r14, [reg_p2+48] 
  sbb   r15, [reg_p2+56]
  mov   [reg_p3+32], r12
  mov   [reg_p3+40], r13
  mov   [reg_p3+48], r14
  mov   [reg_p3+56], r15
  mov   r8 , [reg_p1+64]
  mov   r9 , [reg_p1+72]
  mov   r10, [reg_p1+80]
  mov   r11, [reg_p1+88]
  mov   r12, 0
  sbb   r8 , [reg_p2+64]
  sbb   r9 , [reg_p2+72]
  sbb   r10, [reg_p2+80] 
  sbb   r11, [reg_p2+88]
  sbb   r12, 0
  mov   [reg_p3+64], r8
  mov   [reg_p3+72], r9
  mov   [reg_p3+80], r10
  mov   [reg_p3+88], r11

  mov   rbx, p765x2_0
  and   rbx, r12
  mov   r8 , p765x2_1
  and   r8 , r12
  mov   r9 , p765x2_6
  and   r9 , r12
  mov   r11, p765x2_7
  and   r11, r12
  mov   r10, p765x2_8
  and   r10, r12
  mov   r13, p765x2_9
  and   r13, r12
  mov   r14, p765x2_10
  and   r14, r12
  mov   r15, p765x2_11
  and   r15, r12
  
  mov   rax, [reg_p3]
  add   rax, rbx  
  mov   [reg_p3], rax
  mov   rax, [reg_p3+8]
  adc   rax, r8 
  mov   [reg_p3+8], rax  
  mov   rax, [reg_p3+16]
  adc   rax, r8 
  mov   [reg_p3+16], rax  
  mov   rax, [reg_p3+24]  
  adc   rax, r8 
  mov   [reg_p3+24], rax 
  mov   rax, [reg_p3+32]  
  adc   rax, r8 
  mov   [reg_p3+32], rax 
  mov   rax, [reg_p3+40]    
  adc   rax, r8 
  mov   [reg_p3+40], rax 
  mov   rax, [reg_p3+48]   
  adc   rax, r9 
  mov   [reg_p3+48], rax 
  mov   rax, [reg_p3+56]   
  adc   rax, r11  
  mov   [reg_p3+56], rax 
  mov   rax, [reg_p3+64]  
  adc   rax, r10 
  mov   [reg_p3+64], rax 
  mov   rax, [reg_p3+72]   
  adc   rax, r13 
  mov   [reg_p3+72], rax 
  mov   rax, [reg_p3+80]   
  adc   rax, r14 
  mov   [reg_p3+80], rax 
  mov   rax, [reg_p3+88]   
  adc   rax, r15
  mov   [reg_p3+88], rax 

  pop    rbx
  pop    r15
  pop    r14
  pop    r13
  pop    r12
  ret

///////////////////////////////////////////////////////////////// MACRO
// Schoolbook integer multiplication
// Inputs:  memory pointers M0 and M1
// Outputs: memory pointer C
// Temps:   stack pointer for 15 64-bit values, regs T0:T7
/////////////////////////////////////////////////////////////////
.macro MUL384_SCHOOL M0, M1, C, T0, T1, T2, T3, T4, T5, T6, T7
    mov    rdx, \M0
    mulx   \T1, \T0, \M1
    mulx   \T2, rax, 8\M1
    add    \T1, rax
    mulx   \T3, rax, 16\M1
    adc    \T2, rax
    mulx   \T4, rax, 24\M1
    adc    \T3, rax
    mulx   \T5, rax, 32\M1
    adc    \T4, rax
    mulx   \T6, rax, 40\M1
    adc    \T5, rax
    adc    \T6, 0
    mov    \C, \T0
    
    xor    \T0, \T0
    mov    rdx, 8\M0
    mulx   \T7, rax, \M1
    adcx   \T1, rax
    adox   \T2, \T7
    mulx   \T7, rax, 8\M1
    adcx   \T2, rax
    adox   \T3, \T7
    mulx   \T7, rax, 16\M1
    adcx   \T3, rax
    adox   \T4, \T7
    mulx   \T7, rax, 24\M1
    adcx   \T4, rax
    adox   \T5, \T7
    mulx   \T7, rax, 32\M1
    adcx   \T5, rax
    adox   \T6, \T7
    mulx   \T7, rax, 40\M1
    adcx   \T6, rax
    adox   \T0, \T7
    adc    \T0, 0
    mov   8\C, \T1
    
    xor    \T1, \T1
    mov    rdx, 16\M0
    mulx   \T7, rax, \M1
    adcx   \T2, rax
    adox   \T3, \T7
    mulx   \T7, rax, 8\M1
    adcx   \T3, rax
    adox   \T4, \T7
    mulx   \T7, rax, 16\M1
    adcx   \T4, rax
    adox   \T5, \T7
    mulx   \T7, rax, 24\M1
    adcx   \T5, rax
    adox   \T6, \T7
    mulx   \T7, rax, 32\M1
    adcx   \T6, rax
    adox   \T0, \T7
    mulx   \T7, rax, 40\M1
    adcx   \T0, rax
    adox   \T1, \T7
    adc    \T1, 0
    mov   16\C, \T2

    xor    \T2, \T2
    mov    rdx, 24\M0
    mulx   \T7, rax, \M1
    adcx   \T3, rax
    adox   \T4, \T7
    mulx   \T7, rax, 8\M1
    adcx   \T4, rax
    adox   \T5, \T7
    mulx   \T7, rax, 16\M1
    adcx   \T5, rax
    adox   \T6, \T7
    mulx   \T7, rax, 24\M1
    adcx   \T6, rax
    adox   \T0, \T7
    mulx   \T7, rax, 32\M1
    adcx   \T0, rax
    adox   \T1, \T7
    mulx   \T7, rax, 40\M1
    adcx   \T1, rax
    adox   \T2, \T7
    adc    \T2, 0
    mov   24\C, \T3
    
    xor    \T3, \T3
    mov    rdx, 32\M0
    mulx   \T7, rax, \M1
    adcx   \T4, rax
    adox   \T5, \T7
    mulx   \T7, rax, 8\M1
    adcx   \T5, rax
    adox   \T6, \T7
    mulx   \T7, rax, 16\M1
    adcx   \T6, rax
    adox   \T0, \T7
    mulx   \T7, rax, 24\M1
    adcx   \T0, rax
    adox   \T1, \T7
    mulx   \T7, rax, 32\M1
    adcx   \T1, rax
    adox   \T2, \T7
    mulx   \T7, rax, 40\M1
    adcx   \T2, rax
    adox   \T3, \T7
    adc    \T3, 0
    mov   32\C, \T4

    xor    \T4, \T4
    mov    rdx, 40\M0
    mulx   \T7, rax, \M1
    adcx   \T5, rax
    adox   \T6, \T7
    mulx   \T7, rax, 8\M1
    adcx   \T6, rax
    adox   \T0, \T7
    mulx   \T7, rax, 16\M1
    adcx   \T0, rax
    adox   \T1, \T7
    mulx   \T7, rax, 24\M1
    adcx   \T1, rax
    adox   \T2, \T7
    mulx   \T7, rax, 32\M1
    adcx   \T2, rax
    adox   \T3, \T7
    mulx   \T7, rax, 40\M1
    adcx   \T3, rax
    adox   \T4, \T7
    adc    \T4, 0
    
    mov  40\C, \T5
    mov  48\C, \T6
    mov  56\C, \T0
    mov  64\C, \T1
    mov  72\C, \T2
    mov  80\C, \T3
    mov  88\C, \T4
.endm

///////////////////////////////////////////////////////////////// MACRO
// Schoolbook integer multiplication
// Inputs:  memory pointers M0 and M1
// Outputs: memory pointer C
// Temps:   stack pointer for 15 64-bit values, regs T0:T7
/////////////////////////////////////////////////////////////////
.macro MUL384_SCHOOL_V2 M0, M1, C, T0, T1, T2, T3, T4, T5, T6, T7
    mov    rdx, \M0
    mulx   \T1, \T0, \M1
    mulx   \T2, rax, 8\M1
    add    \T1, rax
    mulx   \T3, rax, 16\M1
    adc    \T2, rax
    mulx   \T4, rax, 24\M1
    adc    \T3, rax
    mulx   \T5, rax, 32\M1
    adc    \T4, rax
    mulx   \T6, rax, 40\M1
    adc    \T5, rax
    adc    \T6, 0
    mov    \C, \T0
    
    xor    \T0, \T0
    mov    rdx, 8\M0
    mulx   \T7, rax, \M1
    adcx   \T1, rax
    adox   \T2, \T7
    mulx   \T7, rax, 8\M1
    adcx   \T2, rax
    adox   \T3, \T7
    mulx   \T7, rax, 16\M1
    adcx   \T3, rax
    adox   \T4, \T7
    mulx   \T7, rax, 24\M1
    adcx   \T4, rax
    adox   \T5, \T7
    mulx   \T7, rax, 32\M1
    adcx   \T5, rax
    adox   \T6, \T7
    mulx   \T7, rax, 40\M1
    adcx   \T6, rax
    adox   \T0, \T7
    adc    \T0, 0
    mov   8\C, \T1
    
    xor    \T1, \T1
    mov    rdx, 16\M0
    mulx   \T7, rax, \M1
    adcx   \T2, rax
    adox   \T3, \T7
    mulx   \T7, rax, 8\M1
    adcx   \T3, rax
    adox   \T4, \T7
    mulx   \T7, rax, 16\M1
    adcx   \T4, rax
    adox   \T5, \T7
    mulx   \T7, rax, 24\M1
    adcx   \T5, rax
    adox   \T6, \T7
    mulx   \T7, rax, 32\M1
    adcx   \T6, rax
    adox   \T0, \T7
    mulx   \T7, rax, 40\M1
    adcx   \T0, rax
    adox   \T1, \T7
    adc    \T1, 0
    mov   16\C, \T2

    xor    \T2, \T2
    mov    rdx, 24\M0
    mulx   \T7, rax, \M1
    adcx   \T3, rax
    adox   \T4, \T7
    mulx   \T7, rax, 8\M1
    adcx   \T4, rax
    adox   \T5, \T7
    mulx   \T7, rax, 16\M1
    adcx   \T5, rax
    adox   \T6, \T7
    mulx   \T7, rax, 24\M1
    adcx   \T6, rax
    adox   \T0, \T7
    mulx   \T7, rax, 32\M1
    adcx   \T0, rax
    adox   \T1, \T7
    mulx   \T7, rax, 40\M1
    adcx   \T1, rax
    adox   \T2, \T7
    adc    \T2, 0
    mov   24\C, \T3
    
    xor    \T3, \T3
    mov    rdx, 32\M0
    mulx   \T7, rax, \M1
    adcx   \T4, rax
    adox   \T5, \T7
    mulx   \T7, rax, 8\M1
    adcx   \T5, rax
    adox   \T6, \T7
    mulx   \T7, rax, 16\M1
    adcx   \T6, rax
    adox   \T0, \T7
    mulx   \T7, rax, 24\M1
    adcx   \T0, rax
    adox   \T1, \T7
    mulx   \T7, rax, 32\M1
    adcx   \T1, rax
    adox   \T2, \T7
    mulx   \T7, rax, 40\M1
    adcx   \T2, rax
    adox   \T3, \T7
    adc    \T3, 0
    mov   32\C, \T4

    xor    \T4, \T4
    mov    rdx, 40\M0
    mulx   \T7, rax, \M1
    adcx   \T5, rax
    adox   \T6, \T7
    mulx   \T7, rax, 8\M1
    adcx   \T6, rax
    adox   \T0, \T7
    mulx   \T7, rax, 16\M1
    adcx   \T0, rax
    adox   \T1, \T7
    mulx   \T7, rax, 24\M1
    adcx   \T1, rax
    adox   \T2, \T7
    mulx   \T7, rax, 32\M1
    adcx   \T2, rax
    adox   \T3, \T7
    mulx   \T7, rax, 40\M1
    adcx   \T3, rax
    adox   \T4, \T7
    adc    \T4, 0
    /*
    mov  40\C, \T5
    mov  48\C, \T6
    mov  56\C, \T0
    mov  64\C, \T1
    mov  72\C, \T2
    mov  80\C, \T3
    mov  88\C, \T4*/  
.endm

///////////////////////////////////////////////////////////////// MACRO
// Karatsuba integer squaring
// Inputs:  memory pointers M0 and M1
// Outputs: memory pointer C
// Temps:   stack pointer for 6 64-bit values, regs T0:T11
/////////////////////////////////////////////////////////////////
.macro SQR384_SCHOOL M0, C, T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11
    // M0_L*M0_L
    mov    rdx, \M0
    mulx   \T1, \T0,   \M0
    mulx   \T3, \T2,  8\M0
    mulx   \T5, \T4, 16\M0
    mov    rdx, 8\M0
    mulx   \T7, \T6,   rdx
    mulx   \T9, \T8, 16\M0
    mov    rdx, 16\M0
    mulx   \T11, \T10,   rdx
    
    add     \T1, \T2
    adc     \T6, \T3
    adc     \T7, \T5
    adc    \T10, \T9
    adc    \T11, 0
    add     \T6, \T4
    adc     \T7, \T8
    adc    \T10, 0
    add     \T1, \T2
    adc     \T6, \T3
    adc     \T7, \T5
    adc    \T10, \T9
    adc    \T11, 0
    add     \T6, \T4
    adc     \T7, \T8
    adc    \T10, 0
        
    mov   \C, \T0
    mov  8\C, \T1
    mov 16\C, \T6
    mov 24\C, \T7
    mov 32\C, \T10
    mov 40\C, \T11

    // M0_H*M0_H
    mov    rdx, 24\M0
    mulx   \T1, \T0, 24\M0
    mulx   \T3, \T2, 32\M0
    mulx   \T5, \T4, 40\M0
    mov    rdx,32\M0
    mulx   \T7, \T6,   rdx
    mulx   \T9, \T8, 40\M0
    mov    rdx, 40\M0
    mulx   \T11, \T10,   rdx
    
    add     \T1, \T2
    adc     \T6, \T3
    adc     \T7, \T5
    adc    \T10, \T9
    adc    \T11, 0
    add     \T6, \T4
    adc     \T7, \T8
    adc    \T10, 0
    add     \T1, \T2
    adc     \T6, \T3
    adc     \T7, \T5
    adc    \T10, \T9
    adc    \T11, 0
    add     \T6, \T4
    adc     \T7, \T8
    adc    \T10, 0
        
    mov 80\C, \T10
    mov 88\C, \T11

    // M0_L*M0_H
    mov    rdx, \M0
    mulx   \T3, \T2, 24\M0
    mulx   \T4, \T8, 32\M0
    add    \T3, \T8
    mulx   \T5, \T8, 40\M0
    adc    \T4, \T8
    adc    \T5, 0
  
    xor    \T9, \T9
    mov    rdx, 8\M0
    mulx   \T10, \T8, 24\M0
    adcx   \T3, \T8
    adox   \T4, \T10
    mulx   \T10, \T8, 32\M0
    adcx   \T4, \T8
    adox   \T5, \T10
    mulx   \T10, \T8, 40\M0
    adcx   \T5, \T8
    adox   \T9, \T10
    adc    \T9, 0

    xor    \T11, \T11
    mov    rdx, 16\M0
    mulx   \T10, \T8, 24\M0
    adcx   \T4, \T8
    adox   \T5, \T10
    mulx   \T10, \T8, 32\M0
    adcx   \T5, \T8
    adox   \T9, \T10
    mulx   \T10, \T8, 40\M0
    adcx   \T9, \T8
    adox   \T11, \T10
    adc    \T11, 0
    
    xor  \T8, \T8
    add  \T2, \T2
    adc  \T3, \T3
    adc  \T4, \T4
    adc  \T5, \T5
    adc  \T9, \T9
    adc  \T11, \T11
    adc  \T8, 0
    add  \T2, 24\C
    adc  \T3, 32\C
    adc  \T4, 40\C
    adc  \T0, \T5
    adc  \T1, \T9
    adc  \T6, \T11
    adc  \T7, \T8
    
    mov  24\C, \T2
    mov  32\C, \T3
    mov  40\C, \T4
    mov  48\C, \T0
    mov  56\C, \T1
    mov  64\C, \T6
    mov  72\C, \T7

.endm

//*****************************************************************************
//  765-bit multiplication using Karatsuba (one level), schoolbook (two levels)
//***************************************************************************** 
.global sqr765_asm
sqr765_asm:    
    push   r12
    push   r13 
    push   r14 
    push   r15
    push   rbx 
    push   rbp
    sub    rsp, 40

    SQR384_SCHOOL [reg_p1], [reg_p2], r8, r9, r10, r11, r12, r13, r14, r15, rax, rbx, rcx, rbp
    SQR384_SCHOOL [reg_p1+48], [reg_p2+96], r8, r9, r10, r11, r12, r13, r14, r15, rax, rbx, rcx, rbp
    
    MUL384_SCHOOL_V2  [reg_p1], [reg_p1+48], [rsp], r8, r9, r10, r11, rbx, rbp, r14, r15

    xor   rdx, rdx
    mov   r13, [rsp]
    mov   r15, [rsp+8]
    mov   rax, [rsp+16]
    mov   r12, [rsp+24]
    mov   rcx, [rsp+32]
    
    add  r13, r13
    adc  r15, r15
    adc  rax, rax
    adc  r12, r12
    adc  rcx, rcx
    adc  rbp, rbp
    adc  r14, r14
    adc  r8 , r8
    adc  r9 , r9
    adc  r10, r10
    adc  r11, r11
    adc  rbx, rbx
    adc  rdx, 0            
    add  r13, [reg_p2+48]
    adc  r15, [reg_p2+56]
    adc  rax, [reg_p2+64]
    adc  r12, [reg_p2+72]
    adc  rcx, [reg_p2+80]
    adc  rbp, [reg_p2+88]
    adc  r14, [reg_p2+96]
    adc  r8 , [reg_p2+104]
    adc  r9 , [reg_p2+112]
    adc  r10, [reg_p2+120]
    adc  r11, [reg_p2+128]
    adc  rbx, [reg_p2+136]
    adc  rdx, [reg_p2+144]

    mov [reg_p2+48], r13
    mov [reg_p2+56], r15
    mov [reg_p2+64], rax
    mov [reg_p2+72], r12
    mov [reg_p2+80], rcx
    mov [reg_p2+88], rbp
    mov [reg_p2+96], r14
    mov [reg_p2+104], r8
    mov [reg_p2+112], r9
    mov [reg_p2+120], r10
    mov [reg_p2+128], r11
    mov [reg_p2+136], rbx
    mov [reg_p2+144], rdx
    
    add    rsp, 40
    pop    rbp  
    pop    rbx
    pop    r15
    pop    r14
    pop    r13
    pop    r12
    ret

//*****************************************************************************
//  765-bit multiplication using Karatsuba (one level), schoolbook (two levels)
//***************************************************************************** 
.global mul765_asm
mul765_asm:    
    push   r12
    push   r13 
    push   r14 
    push   r15
    mov    rcx, reg_p3 

    // [rsp] <- AH + AL, rax <- mask
    xor    rax, rax
    mov    r8, [reg_p1]
    mov    r9, [reg_p1+8]
    mov    r10, [reg_p1+16]
    mov    r11, [reg_p1+24] 
    mov    r12, [reg_p1+32] 
    mov    r13, [reg_p1+40] 
    push   rbx 
    push   rbp
    sub    rsp, 168
    add    r8, [reg_p1+48]
    adc    r9, [reg_p1+56]
    adc    r10, [reg_p1+64]
    adc    r11, [reg_p1+72]
    adc    r12, [reg_p1+80]
    adc    r13, [reg_p1+88]
    sbb    rax, 0
    mov    [rsp], r8
    mov    [rsp+8], r9
    mov    [rsp+16], r10
    mov    [rsp+24], r11
    mov    [rsp+32], r12
    mov    [rsp+40], r13

    // [rsp+48] <- BH + BL, rdx <- mask
    xor    rdx, rdx
    mov    r8, [reg_p2]
    mov    r9, [reg_p2+8]
    mov    rbx, [reg_p2+16]
    mov    rbp, [reg_p2+24] 
    mov    r14, [reg_p2+32]     
    mov    r15, [reg_p2+40]     
    add    r8, [reg_p2+48]
    adc    r9, [reg_p2+56]
    adc    rbx, [reg_p2+64]
    adc    rbp, [reg_p2+72]
    adc    r14, [reg_p2+80]
    adc    r15, [reg_p2+88]
    sbb    rdx, 0
    mov    [rsp+48], r8
    mov    [rsp+56], r9
    mov    [rsp+64], rbx
    mov    [rsp+72], rbp
    mov    [rsp+80], r14     
    mov    [rsp+88], r15     
    
    // [rcx] <- masked (BH + BL)
    and    r8, rax
    and    r9, rax
    and    rbx, rax
    and    rbp, rax
    and    r14, rax     
    and    r15, rax     
    mov    [rcx], r8
    mov    [rcx+8], r9
    mov    [rcx+16], rbx    /////
    mov    [rcx+24], rbp    /////

    // r8-r13 <- masked (AH + AL)
    mov    r8, [rsp]
    mov    r9, [rsp+8]
    and    r8, rdx
    and    r9, rdx
    and    r10, rdx
    and    r11, rdx
    and    r12, rdx
    and    r13, rdx

    // [rsp+96] <- masked (AH + AL) + masked (AH + AL)
    mov    rax, [rcx]
    mov    rdx, [rcx+8]
    add    r8, rax
    adc    r9, rdx
    adc    r10, rbx
    adc    r11, rbp
    adc    r12, r14         
    adc    r13, r15      
    mov    [rsp+96], r8
    mov    [rsp+104], r9
    mov    [rsp+112], r10
    mov    [rsp+120], r11

    // [rcx] <- AL x BL
    MUL384_SCHOOL  [reg_p1], [reg_p2], [rcx], r8, r9, r10, r11, rbx, rbp, r14, r15     // Result C0-C5 

    // [rsp] <- AH x BH 
    MUL384_SCHOOL  [reg_p1+48], [reg_p2+48], [rcx+96], r8, r9, r10, r11, rbx, rbp, r14, r15

    // [rcx+96] <- (AH+AL) x (BH+BL), low part 
    MUL384_SCHOOL_V2  [rsp], [rsp+48], [rsp+128], r8, r9, r10, r11, rbx, rbp, r14, r15
    
    // r8-r13 <- (AH+AL) x (BH+BL), final step
    add    r14, [rsp+96]
    adc    r8 , [rsp+104]
    adc    r9 , [rsp+112]
    adc    r10, [rsp+120]
    adc    r11, r12
    adc    rbx, r13

    mov  [rsp], rdi
    // rdi,rdx,rbx,rbp,r14,r15,r8-r13 <- (AH+AL) x (BH+BL) - ALxBL
    mov    rdi, [rsp+128]
    sub    rdi, [rcx]
    mov    rdx, [rsp+136]
    sbb    rdx, [rcx+8]
    mov    r12, [rsp+144]
    sbb    r12, [rcx+16]
    mov    r13, [rsp+152]
    sbb    r13, [rcx+24]
    mov    r15, [rsp+160]     
    sbb    r15, [rcx+32]   
    sbb    rbp, [rcx+40]     
    sbb    r14, [rcx+48]
    sbb    r8 , [rcx+56]
    sbb    r9 , [rcx+64]
    sbb    r10, [rcx+72]
    sbb    r11, [rcx+80]
    sbb    rbx, [rcx+88]
      
    // rdi,rdx,rbx,rbp,r14,r15,r8-r13 <- (AH+AL) x (BH+BL) - ALxBL - AHxBH
    sub    rdi, [rcx+96]
    sbb    rdx, [rcx+104]
    sbb    r12, [rcx+112]
    sbb    r13, [rcx+120]
    sbb    r15, [rcx+128]     
    sbb    rbp, [rcx+136]   
    sbb    r14, [rcx+144]
    sbb    r8 , [rcx+152]
    sbb    r9 , [rcx+160]
    sbb    r10, [rcx+168]
    sbb    r11, [rcx+176]
    sbb    rbx, [rcx+184]
    
    mov   rax, [rcx+144]
    add   rdi, [rcx+48]
    adc   rdx, [rcx+56]
    adc   r12, [rcx+64]
    adc   r13, [rcx+72]
    adc   r15, [rcx+80]
    adc   rbp, [rcx+88]
    adc   r14, [rcx+96]
    adc   r8 , [rcx+104]
    adc   r9 , [rcx+112]
    adc   r10, [rcx+120]
    adc   r11, [rcx+128]
    adc   rbx, [rcx+136]
    adc   rax, 0
 
    mov [rcx+48], rdi
    mov [rcx+56], rdx
    mov [rcx+64], r12
    mov [rcx+72], r13
    mov [rcx+80], r15
    mov [rcx+88], rbp
    mov [rcx+96], r14
    mov [rcx+104], r8
    mov [rcx+112], r9
    mov [rcx+120], r10
    mov [rcx+128], r11
    mov [rcx+136], rbx
    mov [rcx+144], rax
    mov rdi,  [rsp]
    
    add    rsp, 168   
    pop    rbp  
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
.global rdc765_asm
.p2align 4
rdc765_asm:
  push   r12
  push   r13 
  push   r14 
  push   r15
  sub rsp, 24

  mov   rdx, p765p1_6
  mulx  r9 , r8 , [reg_p1]
  mulx  r10, rax, [reg_p1+8]
  add  [reg_p1+48], r8
  adc  r9 , rax
  mulx  r11, rax, [reg_p1+16]
  adc  r10, rax
  mulx  r12, rax, [reg_p1+24]
  adc  r11, rax
  mulx  r13, rax, [reg_p1+32]
  adc  r12, rax
  mulx  r14, rax, [reg_p1+40]
  adc  r13, rax
  adc  r14, 0

  xor  r15, r15
  mov  rdx, p765p1_7
  mulx  rcx, rax, [reg_p1]
  adcx  r9 , rax
  adox  r10, rcx
  mulx  rcx, rax, [reg_p1+8]
  adcx  r10, rax
  adox  r11, rcx
  mulx  rcx, rax, [reg_p1+16]
  adcx  r11, rax
  adox  r12, rcx
  mulx  rcx, rax, [reg_p1+24]
  adcx  r12, rax
  adox  r13, rcx
  mulx  rcx, rax, [reg_p1+32]
  adcx  r13, rax
  adox  r14, rcx
  mulx  rcx, rax, [reg_p1+40]
  adcx  r14, rax
  adox  r15, rcx
  adc  r15, 0

  xor  r8 , r8
  mov  [rsp], r9
  mov  rdx, p765p1_8
  mulx  rcx, rax, [reg_p1]
  adcx  r10, rax
  adox  r11, rcx
  mulx  rcx, rax, [reg_p1+8]
  adcx  r11, rax
  adox  r12, rcx
  mulx  rcx, rax, [reg_p1+16]
  adcx  r12, rax
  adox  r13, rcx
  mulx  rcx, rax, [reg_p1+24]
  adcx  r13, rax
  adox  r14, rcx
  mulx  rcx, rax, [reg_p1+32]
  adcx  r14, rax
  adox  r15, rcx
  mulx  rcx, rax, [reg_p1+40]
  adcx  r15, rax
  adox  r8 , rcx
  adc  r8 , 0

  xor  r9 , r9
  mov  [rsp+8], r10
  mov  rdx, p765p1_9
  mulx  rcx, rax, [reg_p1]
  adcx  r11, rax
  adox  r12, rcx
  mulx  rcx, rax, [reg_p1+8]
  adcx  r12, rax
  adox  r13, rcx
  mulx  rcx, rax, [reg_p1+16]
  adcx  r13, rax
  adox  r14, rcx
  mulx  rcx, rax, [reg_p1+24]
  adcx  r14, rax
  adox  r15, rcx
  mulx  rcx, rax, [reg_p1+32]
  adcx  r15, rax
  adox  r8 , rcx
  mulx  rcx, rax, [reg_p1+40]
  adcx  r8 , rax
  adox  r9 , rcx
  adc  r9 , 0

  xor  r10, r10
  mov  [rsp+16], r11
  mov  rdx, p765p1_10
  mulx  rcx, rax, [reg_p1]
  adcx  r12, rax
  adox  r13, rcx
  mulx  rcx, rax, [reg_p1+8]
  adcx  r13, rax
  adox  r14, rcx
  mulx  rcx, rax, [reg_p1+16]
  adcx  r14, rax
  adox  r15, rcx
  mulx  rcx, rax, [reg_p1+24]
  adcx  r15, rax
  adox  r8 , rcx
  mulx  rcx, rax, [reg_p1+32]
  adcx  r8 , rax
  adox  r9 , rcx
  mulx  rcx, rax, [reg_p1+40]
  adcx  r9 , rax
  adox  r10, rcx
  adc  r10, 0

  xor  r11, r11
  mov  rdx, p765p1_11
  mulx  rcx, rax, [reg_p1]
  adcx  r13, rax
  adox  r14, rcx
  mulx  rcx, rax, [reg_p1+8]
  adcx  r14, rax
  adox  r15, rcx
  mulx  rcx, rax, [reg_p1+16]
  adcx  r15, rax
  adox  r8 , rcx
  mulx  rcx, rax, [reg_p1+24]
  adcx  r8 , rax
  adox  r9 , rcx
  mulx  rcx, rax, [reg_p1+32]
  adcx  r9 , rax
  adox  r10, rcx
  mulx  rcx, rax, [reg_p1+40]
  adcx  r10, rax
  adox  r11, rcx
  adc  r11, 0

  mov  rax, [rsp]
  mov  rcx, [rsp+8]
  mov  rdx, [rsp+16]

  add  [reg_p1+56] , rax
  adc  [reg_p1+64] , rcx
  adc  [reg_p1+72] , rdx
  adc  [reg_p1+80] , r12
  adc  [reg_p1+88] , r13
  adc  [reg_p1+96] , r14
  adc  [reg_p1+104], r15
  mov  rax, 0
  adc  [reg_p1+112], r8
  adc  [reg_p1+120], r9
  adc  [reg_p1+128], r10
  adc  [reg_p1+136], r11
  adc  [reg_p1+144], rax

  mov  rcx, [reg_p1+96] 
  mov  rdx, p765p1_6
  mulx  r9 , r8 , [reg_p1+48]
  mulx  r10, rax, [reg_p1+56]
  add  rcx, r8
  adc  r9 , rax
  mulx  r11, rax, [reg_p1+64]
  adc  r10, rax
  mulx  r12, rax, [reg_p1+72]
  adc  r11, rax
  mulx  r13, rax, [reg_p1+80]
  adc  r12, rax
  mulx  r14, rax, [reg_p1+88]
  adc  r13, rax
  adc  r14, 0
  mov  [reg_p2], rcx

  xor  r15, r15
  mov  rdx, p765p1_7
  mulx  rcx, rax, [reg_p1+48]
  adcx  r9 , rax
  adox  r10, rcx
  mulx  rcx, rax, [reg_p1+56]
  adcx  r10, rax
  adox  r11, rcx
  mulx  rcx, rax, [reg_p1+64]
  adcx  r11, rax
  adox  r12, rcx
  mulx  rcx, rax, [reg_p1+72]
  adcx  r12, rax
  adox  r13, rcx
  mulx  rcx, rax, [reg_p1+80]
  adcx  r13, rax
  adox  r14, rcx
  mulx  rcx, rax, [reg_p1+88]
  adcx  r14, rax
  adox  r15, rcx
  adc  r15, 0

  xor  r8 , r8
  mov  [rsp], r9
  mov  rdx, p765p1_8
  mulx  rcx, rax, [reg_p1+48]
  adcx  r10, rax
  adox  r11, rcx
  mulx  rcx, rax, [reg_p1+56]
  adcx  r11, rax
  adox  r12, rcx
  mulx  rcx, rax, [reg_p1+64]
  adcx  r12, rax
  adox  r13, rcx
  mulx  rcx, rax, [reg_p1+72]
  adcx  r13, rax
  adox  r14, rcx
  mulx  rcx, rax, [reg_p1+80]
  adcx  r14, rax
  adox  r15, rcx
  mulx  rcx, rax, [reg_p1+88]
  adcx  r15, rax
  adox  r8 , rcx
  adc  r8 , 0

  xor  r9 , r9
  mov  [rsp+8], r10
  mov  rdx, p765p1_9
  mulx  rcx, rax, [reg_p1+48]
  adcx  r11, rax
  adox  r12, rcx
  mulx  rcx, rax, [reg_p1+56]
  adcx  r12, rax
  adox  r13, rcx
  mulx  rcx, rax, [reg_p1+64]
  adcx  r13, rax
  adox  r14, rcx
  mulx  rcx, rax, [reg_p1+72]
  adcx  r14, rax
  adox  r15, rcx
  mulx  rcx, rax, [reg_p1+80]
  adcx  r15, rax
  adox  r8 , rcx
  mulx  rcx, rax, [reg_p1+88]
  adcx  r8 , rax
  adox  r9 , rcx
  adc  r9 , 0

  xor  r10, r10
  mov  [rsp+16], r11
  mov  rdx, p765p1_10
  mulx  rcx, rax, [reg_p1+48]
  adcx  r12, rax
  adox  r13, rcx
  mulx  rcx, rax, [reg_p1+56]
  adcx  r13, rax
  adox  r14, rcx
  mulx  rcx, rax, [reg_p1+64]
  adcx  r14, rax
  adox  r15, rcx
  mulx  rcx, rax, [reg_p1+72]
  adcx  r15, rax
  adox  r8 , rcx
  mulx  rcx, rax, [reg_p1+80]
  adcx  r8 , rax
  adox  r9 , rcx
  mulx  rcx, rax, [reg_p1+88]
  adcx  r9 , rax
  adox  r10, rcx
  adc  r10, 0

  xor  r11, r11
  mov  rdx, p765p1_11
  mulx  rcx, rax, [reg_p1+48]
  adcx  r13, rax
  adox  r14, rcx
  mulx  rcx, rax, [reg_p1+56]
  adcx  r14, rax
  adox  r15, rcx
  mulx  rcx, rax, [reg_p1+64]
  adcx  r15, rax
  adox  r8 , rcx
  mulx  rcx, rax, [reg_p1+72]
  adcx  r8 , rax
  adox  r9 , rcx
  mulx  rcx, rax, [reg_p1+80]
  adcx  r9 , rax
  adox  r10, rcx
  mulx  rcx, rax, [reg_p1+88]
  adcx  r10, rax
  adox  r11, rcx
  adc  r11, 0

  mov  rax, [rsp]
  mov  rcx, [rsp+8]
  mov  rdx, [rsp+16]

  add  rax, [reg_p1+104]
  adc  rcx, [reg_p1+112]
  adc  rdx, [reg_p1+120]
  adc  r12, [reg_p1+128]
  adc  r13, [reg_p1+136]
  adc  r14, [reg_p1+144]
  adc  r15, [reg_p1+152]
  adc  r8 , [reg_p1+160]
  adc  r9 , [reg_p1+168]
  adc  r10, [reg_p1+176]
  adc  r11, [reg_p1+184]

  mov  [reg_p2+8] , rax
  mov  [reg_p2+16], rcx
  mov  [reg_p2+24], rdx
  mov  [reg_p2+32], r12
  mov  [reg_p2+40], r13
  mov  [reg_p2+48], r14
  mov  [reg_p2+56], r15
  mov  [reg_p2+64], r8
  mov  [reg_p2+72], r9
  mov  [reg_p2+80], r10
  mov  [reg_p2+88], r11

  add rsp, 24
  pop    r15
  pop    r14
  pop    r13
  pop    r12
  ret


//***********************************************************************
//  765-bit multiprecision addition
//  Operation: c [reg_p3] = a [reg_p1] + b [reg_p2]
//*********************************************************************** 
.global mp_add765_asm
mp_add765_asm:
  push   r12
  push   r13
  push   r14
  push   r15
  push   rbx
  push   rbp
  
  mov    r8, [reg_p1]
  mov    r9, [reg_p1+8]
  mov    r10, [reg_p1+16]
  mov    r11, [reg_p1+24]
  mov    r12, [reg_p1+32]
  mov    r13, [reg_p1+40]
  mov    r14, [reg_p1+48]
  mov    r15, [reg_p1+56] 
  mov    rax, [reg_p1+64]
  mov    rbx, [reg_p1+72] 
  mov    rcx, [reg_p1+80]  
  mov    rbp, [reg_p1+88] 

  add    r8, [reg_p2] 
  adc    r9, [reg_p2+8] 
  adc    r10, [reg_p2+16] 
  adc    r11, [reg_p2+24] 
  adc    r12, [reg_p2+32] 
  adc    r13, [reg_p2+40] 
  adc    r14, [reg_p2+48] 
  adc    r15, [reg_p2+56]
  adc    rax, [reg_p2+64] 
  adc    rbx, [reg_p2+72]
  adc    rcx, [reg_p2+80]
  adc    rbp, [reg_p2+88]

  mov    [reg_p3], r8
  mov    [reg_p3+8], r9
  mov    [reg_p3+16], r10
  mov    [reg_p3+24], r11
  mov    [reg_p3+32], r12
  mov    [reg_p3+40], r13
  mov    [reg_p3+48], r14
  mov    [reg_p3+56], r15
  mov    [reg_p3+64], rax
  mov    [reg_p3+72], rbx
  mov    [reg_p3+80], rcx
  mov    [reg_p3+88], rbp
  
  pop    rbp
  pop    rbx
  pop    r15
  pop    r14
  pop    r13
  pop    r12
  ret
  
//***********************************************************************
//  2x751-bit multiprecision subtraction
//  Operation: c [reg_p3] = a [reg_p1] - b [reg_p2]. Returns borrow mask
//*********************************************************************** 
.global mp_sub765x2_asm
mp_sub765x2_asm:
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
  mov    rcx, [reg_p1+112]
  sbb    r8, [reg_p2+80]
  sbb    r9, [reg_p2+88]
  sbb    r10, [reg_p2+96] 
  sbb    r11, [reg_p2+104] 
  sbb    rcx, [reg_p2+112]
  mov    [reg_p3+80], r8 
  mov    [reg_p3+88], r9
  mov    [reg_p3+96], r10
  mov    [reg_p3+104], r11
  mov    [reg_p3+112], rcx
  
  mov    r8, [reg_p1+120]
  mov    r9, [reg_p1+128]
  mov    r10, [reg_p1+136]
  mov    r11, [reg_p1+144]
  mov    rcx, [reg_p1+152]
  sbb    r8, [reg_p2+120] 
  sbb    r9, [reg_p2+128] 
  sbb    r10, [reg_p2+136] 
  sbb    r11, [reg_p2+144] 
  sbb    rcx, [reg_p2+152]
  mov    [reg_p3+120], r8
  mov    [reg_p3+128], r9
  mov    [reg_p3+136], r10
  mov    [reg_p3+144], r11
  mov    [reg_p3+152], rcx  
   
  mov    r8, [reg_p1+160]
  mov    r9, [reg_p1+168] 
  mov    r10, [reg_p1+176]  
  mov    r11, [reg_p1+184]
  sbb    r8, [reg_p2+160]
  sbb    r9, [reg_p2+168]
  sbb    r10, [reg_p2+176]
  sbb    r11, [reg_p2+184]
  sbb    rax, 0
  mov    [reg_p3+160], r8
  mov    [reg_p3+168], r9
  mov    [reg_p3+176], r10
  mov    [reg_p3+184], r11
  ret  
  
//***********************************************************************
//  Double 2x765-bit multiprecision subtraction
//  Operation: c [reg_p3] = c [reg_p3] - a [reg_p1] - b [reg_p2]
//*********************************************************************** 
.global mp_dblsub765x2_asm
mp_dblsub765x2_asm:
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
    
  xor    rcx, rcx
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
  sbb    r14, [reg_p1+112] 
  sbb    r15, [reg_p1+120]
  adc    rcx, 0
  sub    r8, [reg_p2+64] 
  sbb    r9, [reg_p2+72] 
  sbb    r10, [reg_p2+80] 
  sbb    r11, [reg_p2+88] 
  sbb    r12, [reg_p2+96] 
  sbb    r13, [reg_p2+104] 
  sbb    r14, [reg_p2+112] 
  sbb    r15, [reg_p2+120]
  adc    rcx, 0
  mov    [reg_p3+64], r8
  mov    [reg_p3+72], r9
  mov    [reg_p3+80], r10
  mov    [reg_p3+88], r11
  mov    [reg_p3+96], r12
  mov    [reg_p3+104], r13
  mov    [reg_p3+112], r14
  mov    [reg_p3+120], r15
  
  mov    r8, [reg_p3+128]
  mov    r9, [reg_p3+136]
  mov    r10, [reg_p3+144]
  mov    r11, [reg_p3+152]
  mov    r12, [reg_p3+160]
  mov    r13, [reg_p3+168]
  mov    r14, [reg_p3+176]
  mov    r15, [reg_p3+184]  
  sub    r8, rcx 
  sbb    r8, [reg_p1+128] 
  sbb    r9, [reg_p1+136] 
  sbb    r10, [reg_p1+144] 
  sbb    r11, [reg_p1+152] 
  sbb    r12, [reg_p1+160] 
  sbb    r13, [reg_p1+168] 
  sbb    r14, [reg_p1+176] 
  sbb    r15, [reg_p1+184]
  sub    r8, [reg_p2+128] 
  sbb    r9, [reg_p2+136] 
  sbb    r10, [reg_p2+144] 
  sbb    r11, [reg_p2+152] 
  sbb    r12, [reg_p2+160] 
  sbb    r13, [reg_p2+168] 
  sbb    r14, [reg_p2+176] 
  sbb    r15, [reg_p2+184]
  mov    [reg_p3+128], r8
  mov    [reg_p3+136], r9
  mov    [reg_p3+144], r10
  mov    [reg_p3+152], r11
  mov    [reg_p3+160], r12
  mov    [reg_p3+168], r13
  mov    [reg_p3+176], r14
  mov    [reg_p3+184], r15
  
  pop    r15
  pop    r14
  pop    r13
  pop    r12
  ret
