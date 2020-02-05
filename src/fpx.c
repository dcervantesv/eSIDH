/********************************************************************************************
* SIDH: an efficient supersingular isogeny cryptography library
*
* Abstract: core functions over GF(p) and GF(p^2)
*********************************************************************************************/
#include <stdio.h>

void print_fp(felm_t a, char *name)
{
    int k;
    felm_t aa;

    from_mont(a, aa);
    
    printf("%s := Fp!0x", name);
    for(k=NWORDS_FIELD-1; k>=0; k--)
        printf("%.16jX", aa[k]);
    printf(";\n");
}

void print_fp2(f2elm_t a, char *name)
{
    int k;
    f2elm_t aa;

    from_fp2mont(a, aa);
    
    printf("%s := Fp2![0x", name);
    for(k=NWORDS_FIELD-1; k>=0; k--)
        printf("%.16jX", aa[0][k]);
    printf(", 0x");
    for(k=NWORDS_FIELD-1; k>=0; k--)
        printf("%.16jX", aa[1][k]);
    printf("];\n");
}


void print_point(point_proj_t a, char *name)
{
    int k;
    f2elm_t aa, bb;

    from_fp2mont(a->X, aa);
    from_fp2mont(a->Z, bb);
    
    printf("%sX := Fp2![0x", name);
    for(k=NWORDS_FIELD-1; k>=0; k--)
        printf("%.16jX", aa[0][k]);
    printf(", 0x");
    for(k=NWORDS_FIELD-1; k>=0; k--)
        printf("%.16jX", aa[1][k]);
    printf("];\n%sZ := Fp2![0x", name);
    for(k=NWORDS_FIELD-1; k>=0; k--)
        printf("%.16jX", bb[0][k]);
    printf(", 0x");
    for(k=NWORDS_FIELD-1; k>=0; k--)
        printf("%.16jX", bb[1][k]);
    printf("];\n");
}


void print_fp_2(felm_t a, char *name)
{
    int k;
    
    printf("%s := 0x", name);
    for(k=NWORDS_FIELD-1; k>=0; k--)
        printf("%.16jX", a[k]);
    printf(";\n");
}

void print_fp2_2(f2elm_t a, char *name)
{
    int k;
    
    printf("%s := [0x", name);
    for(k=NWORDS_FIELD-1; k>=0; k--)
        printf("%.16jX", a[0][k]);
    printf(", 0x");
    for(k=NWORDS_FIELD-1; k>=0; k--)
        printf("%.16jX", a[1][k]);
    printf("];\n");
}


void print_point_2(point_proj_t a, char *name)
{
    int k;
    
    printf("%sX := [0x", name);
    for(k=NWORDS_FIELD-1; k>=0; k--)
        printf("%.16jX", a->X[0][k]);
    printf(", 0x");
    for(k=NWORDS_FIELD-1; k>=0; k--)
        printf("%.16jX", a->X[1][k]);
    printf("];\n%sZ := [0x", name);
    for(k=NWORDS_FIELD-1; k>=0; k--)
        printf("%.16jX", a->Z[0][k]);
    printf(", 0x");
    for(k=NWORDS_FIELD-1; k>=0; k--)
        printf("%.16jX", a->Z[1][k]);
    printf("];\n");
}

__inline void fpcopy(const felm_t a, felm_t c)
{ // Copy a field element, c = a.
    unsigned int i;

    for (i = 0; i < NWORDS_FIELD; i++)
        c[i] = a[i];
}

__inline void fpzero(felm_t a)
{ // Zero a field element, a = 0.
    unsigned int i;

    for (i = 0; i < NWORDS_FIELD; i++)
        a[i] = 0;
}


void to_mont(const felm_t a, felm_t mc)
{ // Conversion to Montgomery representation,
  // mc = a*R^2*R^(-1) mod p = a*R mod p, where a in [0, p-1].
  // The Montgomery constant R^2 mod p is the global value "Montgomery_R2". 

    fpmul_mont(a, (digit_t*)&Montgomery_R2, mc);
}


void from_mont(const felm_t ma, felm_t c)
{ // Conversion from Montgomery representation to standard representation,
  // c = ma*R^(-1) mod p = a mod p, where ma in [0, p-1].
    digit_t one[NWORDS_FIELD] = {0};    

    one[0] = 1;
    fpmul_mont(ma, one, c);
    fpcorrection(c);
}


void copy_words(const digit_t* a, digit_t* c, const unsigned int nwords)
{ // Copy wordsize digits, c = a, where lng(a) = nwords.
    unsigned int i;
        
    for (i = 0; i < nwords; i++) {                      
        c[i] = a[i];
    }
}


void fpmul_mont(const felm_t ma, const felm_t mb, felm_t mc)
{ // Multiprecision multiplication, c = a*b mod p.
    dfelm_t temp = {0};

    mp_mul(ma, mb, temp, NWORDS_FIELD);
    rdc_mont(temp, mc);
}


void fpsqr_mont(const felm_t ma, felm_t mc)
{ // Multiprecision squaring, c = a^2 mod p.
    dfelm_t temp = {0};
#if ((NBITS_FIELD != 751) && (NBITS_FIELD != 503))
    sqr_asm(ma, temp);
#else
    mp_mul(ma, ma, temp, NWORDS_FIELD);
#endif    
    rdc_mont(temp, mc);
    
}


void fpinv_mont(felm_t a)
{ // Field inversion using Montgomery arithmetic, a = a^(-1)*R mod p.
    felm_t tt;

    fpcopy(a, tt);
    fpinv_chain_mont(tt);
    fpsqr_mont(tt, tt);
    fpsqr_mont(tt, tt);
    fpmul_mont(a, tt, a);
}


void fp2copy(const f2elm_t a, f2elm_t c)
{ // Copy a GF(p^2) element, c = a.
    fpcopy(a[0], c[0]);
    fpcopy(a[1], c[1]);
}


void fp2zero(f2elm_t a)
{ // Zero a GF(p^2) element, a = 0.
    fpzero(a[0]);
    fpzero(a[1]);
}


void fp2neg(f2elm_t a)
{ // GF(p^2) negation, a = -a in GF(p^2).
    fpneg(a[0]);
    fpneg(a[1]);
}


__inline void fp2add(const f2elm_t a, const f2elm_t b, f2elm_t c)           
{ // GF(p^2) addition, c = a+b in GF(p^2).

    fpadd(a[0], b[0], c[0]);
    fpadd(a[1], b[1], c[1]);

}


__inline void fp2sub(const f2elm_t a, const f2elm_t b, f2elm_t c)          
{ // GF(p^2) subtraction, c = a-b in GF(p^2).
    fpsub(a[0], b[0], c[0]);
    fpsub(a[1], b[1], c[1]);
}


void fp2div2(const f2elm_t a, f2elm_t c)          
{ // GF(p^2) division by two, c = a/2  in GF(p^2).
    fpdiv2(a[0], c[0]);
    fpdiv2(a[1], c[1]);
}


void fp2correction(f2elm_t a)
{ // Modular correction, a = a in GF(p^2).
    fpcorrection(a[0]);
    fpcorrection(a[1]);
}


__inline static void mp_addfast(const digit_t* a, const digit_t* b, digit_t* c)
{ // Multiprecision addition, c = a+b.    
#if (OS_TARGET == OS_WIN) || defined(GENERIC_IMPLEMENTATION)

    mp_add(a, b, c, NWORDS_FIELD);
    
#elif (OS_TARGET == OS_LINUX)                 
    
    mp_add_asm(a, b, c);    

#endif
}


void fp2sqr_mont(const f2elm_t a, f2elm_t c)
{ // GF(p^2) squaring using Montgomery arithmetic, c = a^2 in GF(p^2).
  // Inputs: a = a0+a1*i, where a0, a1 are in [0, 2*p-1] 
  // Output: c = c0+c1*i, where c0, c1 are in [0, 2*p-1] 
    felm_t t1, t2, t3;
    
    mp_addfast(a[0], a[1], t1);                      // t1 = a0+a1 
    fpsub(a[0], a[1], t2);                           // t2 = a0-a1
    mp_addfast(a[0], a[0], t3);                      // t3 = 2a0
    fpmul_mont(t1, t2, c[0]);                        // c0 = (a0+a1)(a0-a1)
    fpmul_mont(t3, a[1], c[1]);                      // c1 = 2a0*a1
}


__inline unsigned int mp_sub(const digit_t* a, const digit_t* b, digit_t* c, const unsigned int nwords)
{ // Multiprecision subtraction, c = a-b, where lng(a) = lng(b) = nwords. Returns the borrow bit.
    unsigned int i, borrow = 0;

    for (i = 0; i < nwords; i++) {
        SUBC(borrow, a[i], b[i], borrow, c[i]);
    }

    return borrow;
}


__inline static digit_t mp_subfast(const digit_t* a, const digit_t* b, digit_t* c)
{ // Multiprecision subtraction, c = a-b, where lng(a) = lng(b) = 2*NWORDS_FIELD. 
  // If c < 0 then returns mask = 0xFF..F, else mask = 0x00..0   
#if (OS_TARGET == OS_WIN) || defined(GENERIC_IMPLEMENTATION)

    return (0 - (digit_t)mp_sub(a, b, c, 2*NWORDS_FIELD));

#elif (OS_TARGET == OS_LINUX)                 

    return mp_subx2_asm(a, b, c);

#endif
}


__inline static void mp_dblsubfast(const digit_t* a, const digit_t* b, digit_t* c)
{ // Multiprecision subtraction, c = c-a-b, where lng(a) = lng(b) = 2*NWORDS_FIELD. 
  // Inputs should be s.t. c > a and c > b  
#if (OS_TARGET == OS_WIN) || defined(GENERIC_IMPLEMENTATION)

    mp_sub(c, a, c, 2*NWORDS_FIELD);
    mp_sub(c, b, c, 2*NWORDS_FIELD);

#elif (OS_TARGET == OS_LINUX)                 

    mp_dblsubx2_asm(a, b, c);

#endif
}


void fp2mul_mont(const f2elm_t a, const f2elm_t b, f2elm_t c)
{ // GF(p^2) multiplication using Montgomery arithmetic, c = a*b in GF(p^2).
  // Inputs: a = a0+a1*i and b = b0+b1*i, where a0, a1, b0, b1 are in [0, 2*p-1] 
  // Output: c = c0+c1*i, where c0, c1 are in [0, 2*p-1] 

    felm_t t1, t2;
    dfelm_t tt1, tt2, tt3; 
    digit_t mask;
    unsigned int i;
    
    mp_addfast(a[0], a[1], t1);                      // t1 = a0+a1
    mp_addfast(b[0], b[1], t2);                      // t2 = b0+b1
    mp_mul(a[0], b[0], tt1, NWORDS_FIELD);           // tt1 = a0*b0
    mp_mul(a[1], b[1], tt2, NWORDS_FIELD);           // tt2 = a1*b1
    mp_mul(t1, t2, tt3, NWORDS_FIELD);               // tt3 = (a0+a1)*(b0+b1)
    mp_dblsubfast(tt1, tt2, tt3);                    // tt3 = (a0+a1)*(b0+b1) - a0*b0 - a1*b1 
    
    mask = mp_subfast(tt1, tt2, tt1);                // tt1 = a0*b0 - a1*b1. If tt1 < 0 then mask = 0xFF..F, else if tt1 >= 0 then mask = 0x00..0

    for (i = 0; i < NWORDS_FIELD; i++) {
        t1[i] = ((digit_t*)PRIME)[i] & mask;
    }

    rdc_mont(tt3, c[1]);                             // c[1] = (a0+a1)*(b0+b1) - a0*b0 - a1*b1 
    mp_addfast((digit_t*)&tt1[NWORDS_FIELD], t1, (digit_t*)&tt1[NWORDS_FIELD]);
    rdc_mont(tt1, c[0]);                             // c[0] = a0*b0 - a1*b1
}


void fpinv_chain_mont(felm_t a)
{ // Chain to compute a^(p-3)/4 using Montgomery arithmetic.
    unsigned int i, j;
    
#ifdef P503_3_1
    felm_t t[15], tt;

    // Precomputed table
    fpsqr_mont(a, tt);
    fpmul_mont(a, tt, t[0]);
    for (i = 0; i <= 13; i++) fpmul_mont(t[i], tt, t[i+1]);

    fpcopy(a, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(a, tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[8], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[6], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[9], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[0], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(a, tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[6], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[8], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(a, tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[10], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[0], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[10], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[10], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[6], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[3], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i = 0; i < 12; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[12], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[8], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[6], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[12], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[11], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[6], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[14], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[14], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[6], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[8], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(a, tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[4], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[6], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[7], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(a, tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[0], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[11], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[13], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[1], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[10], tt, tt);
    for (j = 0; j < 49; j++) {
        for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
        fpmul_mont(t[14], tt, tt);
    }
    fpcopy(tt, a);  

#endif


#ifdef P509_3_5_1_old
    felm_t t[14], tt;

    // Precomputed table
    fpsqr_mont(a, tt);
    fpmul_mont(a, tt, t[0]);
    fpmul_mont(t[0], tt, t[1]);
    fpmul_mont(t[1], tt, t[2]);
    fpmul_mont(t[2], tt, t[3]);
    fpmul_mont(t[3], tt, t[3]);    
    for (i = 3; i <= 12; i++) fpmul_mont(t[i], tt, t[i+1]);

    fpcopy(t[5], tt);
    for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[12], tt, tt);
    for (i=0; i < 3; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[13], tt, tt);
    for (i=0; i < 3; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[1], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[9], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[13], tt, tt);
    for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[10], tt, tt);
    for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[3], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[12], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[13], tt, tt);
    for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[10], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[11], tt, tt);
    fpsqr_mont(tt, tt);
    fpmul_mont(a, tt, tt);
    for (i=0; i < 12; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[8], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[6], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[4], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[13], tt, tt);
    for (i=0; i < 3; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[1], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[3], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[9], tt, tt);
    for (i=0; i < 3; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[0], tt, tt);
    for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[10], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i=0; i < 2; i++) fpsqr_mont(tt, tt);
    fpmul_mont(a, tt, tt);
    for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[7], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[12], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[8], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[6], tt, tt);
    for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[4], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[4], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[1], tt, tt);
    for (i=0; i < 15; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[7], tt, tt);
    for (i=0; i < 3; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[12], tt, tt);
    for (j = 0; j < 50; j++) {
        for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
        fpmul_mont(t[13], tt, tt);
    }
    for (i=0; i < 3; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, a);

#endif

#ifdef P509_3_5_1
    felm_t t[14], tt;

    // Precomputed table
    fpsqr_mont(a, tt);
    fpmul_mont(a, tt, t[0]);
    for (i = 0; i <= 9; i++) fpmul_mont(t[i], tt, t[i+1]);
    fpmul_mont(t[10], tt, t[11]);
    fpmul_mont(t[11], tt, t[11]);
    fpmul_mont(t[11], tt, t[12]);
    fpmul_mont(t[12], tt, t[13]);

fpcopy(t[13], tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[4], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[1], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[0], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[10], tt, tt);
for (i=0; i < 2; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[0], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[0], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[0], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[1], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[1], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[4], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[5], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[10], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 2; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[0], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[6], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[8], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[3], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[5], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[1], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[2], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[2], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[2], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[4], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[11], tt, tt);
for (i=0; i < 3; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[2], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[3], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[6], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[9], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[8], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[4], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[4], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[12], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[2], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[6], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[2], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[7], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[9], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[11], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, a);

#endif

#ifdef P509_3_7_5_0
    felm_t t[15], tt;

    // Precomputed table
    fpsqr_mont(a, tt);
    fpmul_mont(tt, tt, t[0]);
    fpmul_mont(t[0], a, t[0]);
    for (i = 0; i <= 12; i++) fpmul_mont(t[i], tt, t[i+1]);

fpcopy(t[7], tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[11], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[0], tt, tt);
for (i=0; i < 10; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[7], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[4], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[0], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[2], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[10], tt, tt);
for (i=0; i < 3; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[0], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[0], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[4], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[0], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[9], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[11], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[12], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[12], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[0], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[11], tt, tt);
for (i=0; i < 3; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[0], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[0], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[9], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[6], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[8], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[4], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[5], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[10], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[1], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[5], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[1], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[9], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[10], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[4], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[0], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[1], tt, tt);
for (i=0; i < 15; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[9], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[0], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[8], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[3], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[11], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[7], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[9], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[2], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[2], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[12], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[9], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 1; i++) fpsqr_mont(tt, tt);
fpmul_mont(a, tt, a);

#endif

#ifdef P434_3_1

    felm_t t[22], tt, t4;

    // Precomputed table
    fpsqr_mont(a, tt);
    fpsqr_mont(tt, t4);

    fpmul_mont(t4, a, t[0]); //a^5
    fpmul_mont(t[0], tt, t[1]); //a^7

    fpmul_mont(t[1], tt, t[2]);
    fpmul_mont(t[2], tt, t[2]); //a^11
    fpmul_mont(t[2], tt, t[3]); //a^13
    fpmul_mont(t[3], tt, t[4]); //a^15
    
    fpmul_mont(t[4], t4, t[5]);
    fpmul_mont(t[5], tt, t[5]); //a^21
    fpmul_mont(t[5], tt, t[6]); //a^23
    fpmul_mont(t[6], tt, t[7]); //a^25
    
    fpmul_mont(t[7], t4, t[8]); //a^29
    fpmul_mont(t[8], tt, t[9]); //a^31

    fpmul_mont(t[9], t4, t[10]); //a^35
    
    fpmul_mont(t[10], t4, t[11]); //a^39
    fpmul_mont(t[11], tt, t[12]); //a^41
    fpmul_mont(t[12], tt, t[13]); //a^43
    
    fpmul_mont(t[13], t4, t[14]); //a^47
    fpmul_mont(t[14], tt, t[15]); //a^49
    fpmul_mont(t[15], tt, t[16]); //a^51
    fpmul_mont(t[16], tt, t[17]); //a^53
    fpmul_mont(t[17], tt, t[18]); //a^55
    fpmul_mont(t[18], tt, t[19]); //a^57

    fpmul_mont(t[19], t4, t[20]); //a^61
    fpmul_mont(t[20], tt, t[21]); //a^63
    

fpcopy(t[10], tt);
for (i=0; i < 2; i++) fpsqr_mont(tt, tt);
fpmul_mont(a, tt, tt);
for (i=0; i < 10; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[9], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[11], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[6], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[19], tt, tt);
for (i=0; i < 3; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[0], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[10], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[16], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[4], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[8], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[3], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[1], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[12], tt, tt);
for (i=0; i < 12; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[20], tt, tt);
for (i=0; i < 3; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[1], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[7], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[1], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[4], tt, tt);
for (i=0; i < 11; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[15], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[2], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[5], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[18], tt, tt);
for (i=0; i < 10; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[6], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[16], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[17], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[15], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[4], tt, a);


#endif

#ifdef P509_3_7_5_1
    felm_t t[15], tt;

    // Precomputed table
    fpsqr_mont(a, tt);
    fpmul_mont(tt, a, t[0]);
    for (i = 0; i <= 10; i++) fpmul_mont(t[i], tt, t[i+1]);
    fpmul_mont(t[11], tt, t[12]);
    fpmul_mont(t[12], tt, t[12]);
    fpmul_mont(t[12], tt, t[13]);

fpcopy(t[0], tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[7], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[6], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[8], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[8], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[5], tt, tt);
for (i=0; i < 3; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[0], tt, tt);
for (i=0; i < 10; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[11], tt, tt);
for (i=0; i < 3; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[1], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(a, tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[12], tt, tt);
for (i=0; i < 1; i++) fpsqr_mont(tt, tt);
fpmul_mont(a, tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[5], tt, tt);
for (i=0; i < 3; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[0], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[10], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[4], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[9], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[6], tt, tt);
for (i=0; i < 3; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[0], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[3], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[1], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[5], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[2], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[9], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[10], tt, tt);
for (i=0; i < 3; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[1], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[9], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[10], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[12], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[4], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[7], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[8], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[2], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[6], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[3], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[8], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[4], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[1], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[11], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[10], tt, tt);
for (i=0; i < 3; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[2], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[2], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[8], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 2; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[0], tt, a);

#endif

#ifdef P507_3_5_1
    felm_t t[14], tt;

    // Precomputed table
    fpsqr_mont(a, tt);
    fpmul_mont(tt, a, t[0]);
    fpmul_mont(t[0], tt, t[1]);
    fpmul_mont(t[1], tt, t[2]);
    fpmul_mont(t[2], tt, t[2]);
    for (i = 2; i <= 12; i++) fpmul_mont(t[i], tt, t[i+1]);

fpcopy(t[0], tt);
for (i=0; i < 10; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[9], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[12], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[2], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[12], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[6], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[8], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[6], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[6], tt, tt);
for (i=0; i < 2; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[0], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[4], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[9], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[11], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[6], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[10], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[9], tt, tt);
fpsqr_mont(tt, tt);
fpmul_mont(a, tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[2], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[10], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[1], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[1], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[1], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[1], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[3], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[12], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[4], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[11], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[5], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[4], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[3], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[8], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[7], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[6], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[12], tt, tt);
for (i=0; i < 2; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[0], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[9], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[10], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[8], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[0], tt, tt);
for (i=0; i < 10; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[9], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
fpsqr_mont(tt, tt);
fpmul_mont(a, tt, a);

#endif

#ifdef P508_3_5_1
    felm_t t[15], tt;

    // Precomputed table
    fpsqr_mont(a, tt);
    fpmul_mont(tt, a, t[0]);
    for (i = 0; i <= 13; i++) fpmul_mont(t[i], tt, t[i+1]);
    
/*
t[0] := a^3;
t[1] := a^5;
t[2] := a^7;
t[3] := a^9;
t[4] := a^11;
t[5] := a^13;
t[6] := a^15;
t[7] := a^17;
t[8] := a^19;
t[9] := a^21;
t[10] := a^23;
t[11] := a^25;
t[12] := a^27;
t[13] := a^29;
t[14] := a^31;
*/
fpcopy(t[0], tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[2], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[7], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[10], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[9], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[8], tt, tt);
for (i=0; i < 2; i++) fpsqr_mont(tt, tt);
fpmul_mont(a, tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[6], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[4], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[8], tt, tt);
for (i=0; i < 3; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[2], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[3], tt, tt);
for (i=0; i < 10; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[10], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[11], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[7], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[4], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[2], tt, tt);
for (i=0; i < 3; i++) fpsqr_mont(tt, tt);
fpmul_mont(a, tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[11], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[9], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[7], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[8], tt, tt);
for (i=0; i < 3; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[1], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[0], tt, tt);
for (i=0; i < 12; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[12], tt, tt);
for (i=0; i < 2; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[0], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[1], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[12], tt, tt);
for (i=0; i < 2; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[0], tt, tt);
for (i=0; i < 10; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[3], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[3], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[6], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[1], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[5], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[3], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[1], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[6], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[5], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[11], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[12], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
fpsqr_mont(tt, tt);
fpmul_mont(a, tt, a);

#endif


#ifdef P751_3_1
    felm_t t[27], tt;
    
    // Precomputed table
    fpsqr_mont(a, tt);
    fpmul_mont(a, tt, t[0]);
    fpmul_mont(t[0], tt, t[1]);
    fpmul_mont(t[1], tt, t[2]);
    fpmul_mont(t[2], tt, t[3]); 
    fpmul_mont(t[3], tt, t[3]);
    for (i = 3; i <= 8; i++) fpmul_mont(t[i], tt, t[i+1]);
    fpmul_mont(t[9], tt, t[9]);
    for (i = 9; i <= 20; i++) fpmul_mont(t[i], tt, t[i+1]);
    fpmul_mont(t[21], tt, t[21]); 
    for (i = 21; i <= 24; i++) fpmul_mont(t[i], tt, t[i+1]); 
    fpmul_mont(t[25], tt, t[25]);
    fpmul_mont(t[25], tt, t[26]);

    fpcopy(a, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[20], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[24], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[11], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[8], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[23], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i = 0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i = 0; i < 10; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[15], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[13], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[26], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[20], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[11], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[10], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[14], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[4], tt, tt);
    for (i = 0; i < 10; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[18], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[1], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[22], tt, tt);
    for (i = 0; i < 10; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[6], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[24], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[9], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[18], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[17], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(a, tt, tt);
    for (i = 0; i < 10; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[16], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[7], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[0], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[12], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[19], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[22], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[25], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[10], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[22], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[18], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[4], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[14], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[13], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[23], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[21], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[23], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[12], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[9], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[3], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[13], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[17], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[26], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[8], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[11], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[20], tt, tt);
    for (j = 0; j < 61; j++) {
        for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
        fpmul_mont(t[26], tt, tt);
    }
    fpcopy(tt, a);

#endif

#ifdef P765_3_5_1

    felm_t t[26], tt, t4;
    
    // Precomputed table
    fpsqr_mont(a, tt);
	fpsqr_mont(tt, t4);
    fpmul_mont(tt, t4, t[0]);
	fpmul_mont(t[0], a, t[0]);
    fpmul_mont(t[0], t4, t[1]);
    for (i = 1; i <= 4; i++) fpmul_mont(t[i], tt, t[i+1]);
    fpmul_mont(t[5], t4, t[6]);
	for (i = 6; i <= 12; i++) fpmul_mont(t[i], tt, t[i+1]);
    fpmul_mont(t[13], t4, t[14]);
	fpmul_mont(t[14], tt, t[15]);
	fpmul_mont(t[15], tt, t[16]);
    fpmul_mont(t[16], t4, t[17]);
    for (i = 17; i <= 23; i++) fpmul_mont(t[i], tt, t[i+1]);

    fpsqr_mont(t[18], tt);
	fpmul_mont(a, tt, tt);
    for (i = 0; i < 11; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[10], tt, tt);
    for (i = 0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[20], tt, tt);
    for (i = 0; i < 4; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[1], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[15], tt, tt);
    for (i = 0; i < 4; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[9], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[16], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[21], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[24], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[22], tt, tt);
    for (i = 0; i < 3; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[0], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[9], tt, tt);
    for (i = 0; i < 11; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[16], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[4], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[3], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[17], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[4], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[14], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[13], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[22], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[12], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[0], tt, tt);
    for (i = 0; i < 12; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[17], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[10], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[20], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[15], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[19], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[11], tt, tt);
    for (i = 0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[0], tt, tt);
    for (i = 0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[20], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[8], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[7], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[13], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[19], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[7], tt, tt);
    for (i = 0; i < 10; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[16], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[6], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[20], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[24], tt, tt);
    for (i = 0; i < 4; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[1], tt, tt);
    for (i = 0; i < 10; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[19], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[23], tt, tt);
    for (i = 0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[20], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[21], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[15], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[23], tt, tt);
    for (i = 0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[15], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[20], tt, tt);
    for (j = 0; j < 63; j++) {
        for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
        fpmul_mont(t[24], tt, tt);
    }
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[10], tt, a);

#endif


#ifdef P765_3_7_5

    felm_t t[23], tt, t4;
    
    //print_fp(a, "a");

    // Precomputed table
    fpsqr_mont(a, tt);
    fpsqr_mont(tt, t4);
    fpmul_mont(tt, a, t[0]);
    fpmul_mont(t[0], tt, t[1]);
    fpmul_mont(t[1], t4, t[2]);
    fpmul_mont(t[2], tt, t[3]);
    fpmul_mont(t[3], t4, t[4]);
    for (i = 4; i <= 9; i++) fpmul_mont(t[i], tt, t[i+1]);
    fpmul_mont(t[10], t4, t[11]);
    fpmul_mont(t[11], tt, t[12]);
    fpmul_mont(t[12], tt, t[13]);
    fpmul_mont(t[13], t4, t[14]);
    fpmul_mont(t[14], tt, t[15]);
    fpmul_mont(t[15], tt, t[16]);
    fpmul_mont(t[16], t4, t[17]);
    fpmul_mont(t[17], t4, t[18]);
    fpmul_mont(t[18], tt, t[19]);
    fpmul_mont(t[19], t4, t[20]);
    fpmul_mont(t[20], tt, t[21]);
    fpmul_mont(t[21], t4, t[22]);	

    //print_fp(t[22], "b");    
/*
t[0] := a^3;
t[1] := a^5;
t[2] := a^9;
t[3] := a^11;
t[4] := a^15;
t[5] := a^17;
t[6] := a^19;
t[7] := a^21;
t[8] := a^23;
t[9] := a^25;
t[10] := a^27;
t[11] := a^31;
t[12] := a^33;
t[13] := a^35;
t[14] := a^39;
t[15] := a^41;
t[16] := a^43;
t[17] := a^47;
t[18] := a^51;
t[19] := a^53;
t[20] := a^57;
t[21] := a^59;
t[22] := a^63;
*/
fpcopy(t[7], tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[5], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[1], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[3], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[10], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[20], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[12], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[18], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[18], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[20], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[12], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[4], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[9], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(a, tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[2], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[5], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[6], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[5], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[1], tt, tt);
for (i=0; i < 10; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[20], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[9], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[4], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[2], tt, tt);
for (i=0; i < 3; i++) fpsqr_mont(tt, tt);
fpmul_mont(a, tt, tt);
for (i=0; i < 11; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[1], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(a, tt, tt);
for (i=0; i < 11; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[19], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[8], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[15], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[7], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[7], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[12], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[16], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[3], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(a, tt, tt);
for (i=0; i < 14; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[17], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 3; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[1], tt, tt);
for (i=0; i < 12; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[2], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[0], tt, tt);
for (i=0; i < 11; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[20], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[20], tt, tt);
for (i=0; i < 3; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[0], tt, tt);
for (i=0; i < 10; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[9], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
fpsqr_mont(tt, tt);
fpmul_mont(a, tt, tt);
for (i=0; i < 12; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[11], tt, a);


#endif

#ifdef P764_3_5_1

    felm_t t[24], tt, t4;

    // Precomputed table
    fpsqr_mont(a, tt);
    fpsqr_mont(tt, t4);
    fpmul_mont(t4, a, t[0]);
    fpmul_mont(t[0], t4, t[1]);
    for (i = 1; i <= 6; i++) fpmul_mont(t[i], tt, t[i+1]);
    fpmul_mont(t[7], t4, t[8]);
    fpmul_mont(t[8], t4, t[9]);
    fpmul_mont(t[9], tt, t[10]);
    fpmul_mont(t[10], tt, t[11]);
    fpmul_mont(t[11], tt, t[12]);
    fpmul_mont(t[12], t4, t[13]);
    fpmul_mont(t[13], t4, t[14]);
    fpmul_mont(t[14], tt, t[15]);
    fpmul_mont(t[15], tt, t[16]);
    fpmul_mont(t[16], tt, t[17]);
    fpmul_mont(t[17], t4, t[18]);
    fpmul_mont(t[18], tt, t[19]);
    fpmul_mont(t[19], tt, t[20]);
    fpmul_mont(t[20], tt, t[21]);
    fpmul_mont(t[21], tt, t[22]);
    fpmul_mont(t[22], tt, t[23]);

fpcopy(t[21], tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[18], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[3], tt, tt);
for (i=0; i < 3; i++) fpsqr_mont(tt, tt);
fpmul_mont(a, tt, tt);
for (i=0; i < 13; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[17], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[3], tt, tt);
for (i=0; i < 10; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[20], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[8], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[18], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[15], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[17], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[10], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[11], tt, tt);
fpsqr_mont(tt, tt);
fpmul_mont(a, tt, tt);
for (i=0; i < 10; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[9], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[3], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[4], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 3; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[0], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(a, tt, tt);
for (i=0; i < 13; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[19], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[9], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[6], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[12], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[12], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[7], tt, tt);
for (i=0; i < 11; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[17], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[5], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[1], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[4], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[6], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[7], tt, tt);
for (i=0; i < 11; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[15], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[16], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[17], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[8], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[5], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[2], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[17], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 3; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[0], tt, tt);
for (i=0; i < 14; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[20], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[15], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[10], tt, a);

#endif

#ifdef P764_3_5_23

    felm_t t[27], tt, t4;

    // Precomputed table
    fpsqr_mont(a, tt);
    fpsqr_mont(tt, t4);
    fpmul_mont(tt, a, t[0]);
    for (i = 0; i <= 5; i++) fpmul_mont(t[i], tt, t[i+1]);
    fpmul_mont(t[6], t4, t[7]);
    fpmul_mont(t[7], tt, t[8]);
    fpmul_mont(t[8], tt, t[9]);
    fpmul_mont(t[9], t4, t[10]);
    for (i = 10; i <=18 ; i++) fpmul_mont(t[i], tt, t[i+1]);
    fpmul_mont(t[19], t4, t[20]);
    fpmul_mont(t[20], t4, t[21]);
    for (i = 21; i <= 25; i++) fpmul_mont(t[i], tt, t[i+1]);
    
fpcopy(t[10], tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[10], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[15], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[18], tt, tt);
for (i=0; i < 3; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[2], tt, tt);
for (i=0; i < 13; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[19], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[9], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[12], tt, tt);
for (i=0; i < 10; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[16], tt, tt);
for (i=0; i < 10; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[20], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[8], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[1], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[19], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[24], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[3], tt, tt);
for (i=0; i < 3; i++) fpsqr_mont(tt, tt);
fpmul_mont(a, tt, tt);
for (i=0; i < 13; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[25], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[9], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[12], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[3], tt, tt);
for (i=0; i < 12; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[11], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[6], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[5], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[3], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[5], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[11], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[16], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[17], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[19], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[15], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[10], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[7], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(a, tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[0], tt, tt);
for (i=0; i < 10; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[9], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[19], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[19], tt, tt);
fpsqr_mont(tt, tt);
fpmul_mont(a, tt, tt);
for (i=0; i < 10; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[11], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[6], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[8], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[25], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[4], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[20], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[9], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[12], tt, a);

#endif    
    

#ifdef P764_3_7_1

    felm_t t[27], tt, t4;

    // Precomputed table
    fpsqr_mont(a, tt);
    fpsqr_mont(tt, t4);
    fpmul_mont(tt, a, t[0]);
    for (i = 0; i <= 9; i++) fpmul_mont(t[i], tt, t[i+1]);
    fpmul_mont(t[10], t4, t[11]);
    for (i = 11; i <= 14; i++) fpmul_mont(t[i], tt, t[i+1]);
    fpmul_mont(t[15], t4, t[16]);
    for (i = 16; i <= 22; i++) fpmul_mont(t[i], tt, t[i+1]);
    fpmul_mont(t[23], t4, t[24]);
    fpmul_mont(t[24], tt, t[25]);
    fpmul_mont(t[25], t4, t[26]);
    

fpcopy(t[11], tt);
for (i=0; i < 10; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[7], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[24], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[4], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[12], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[18], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[25], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[8], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[4], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[15], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[6], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[17], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[15], tt, tt);
for (i=0; i < 10; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[11], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 10; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[16], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[4], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[11], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[16], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[7], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[17], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[20], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[3], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[7], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[10], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[5], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[25], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[15], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[5], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[4], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[1], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[7], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[8], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[8], tt, tt);
for (i=0; i < 3; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[0], tt, tt);
for (i=0; i < 10; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[19], tt, tt);
for (i=0; i < 2; i++) fpsqr_mont(tt, tt);
fpmul_mont(a, tt, tt);
for (i=0; i < 12; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[25], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[5], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[9], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[24], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[2], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[6], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[15], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 2; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[0], tt, a);

#endif
    
    

/*   764 con 3^e3*5^e5    
    felm_t t[27], tt, t4;

    // Precomputed table
    fpsqr_mont(a, tt);
    fpsqr_mont(tt, t4);
    fpmul_mont(tt, a, t[0]);
    for (i = 0; i <= 8; i++) fpmul_mont(t[i], tt, t[i+1]);
    fpmul_mont(t[9], t4, t[10]);
    for (i = 10; i <= 17; i++) fpmul_mont(t[i], tt, t[i+1]);
    fpmul_mont(t[18], t4, t[19]);
    for (i = 19; i <= 23; i++) fpmul_mont(t[i], tt, t[i+1]);
    fpmul_mont(t[24], t[1], t[25]);
    fpmul_mont(t[25], a, t[25]);
    fpmul_mont(t[25], tt, t[26]);


    fpcopy(t[7], tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[3], tt, tt);
    for (i = 0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[13], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i = 0; i < 11; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[17], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[13], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[10], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[15], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[26], tt, tt);
    for (i = 0; i < 2; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[0], tt, tt);
    for (i = 0; i < 11; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[22], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[16], tt, tt);
    for (i = 0; i < 4; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[3], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[1], tt, tt);
    for (i = 0; i < 10; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[10], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[13], tt, tt);
    for (i = 0; i < 2; i++) fpsqr_mont(tt, tt);
    fpmul_mont(a, tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(a, tt, tt);
    for (i = 0; i < 13; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[11], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[23], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[14], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[19], tt, tt);
    for (i = 0; i < 4; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[14], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[22], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[11], tt, tt);
    for (i = 0; i < 10; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[9], tt, tt);
    for (i = 0; i < 3; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[0], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[0], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[6], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[10], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[8], tt, tt);
    for (i = 0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[21], tt, tt);
    for (i = 0; i < 3; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i = 0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[12], tt, tt);
    for (i = 0; i < 12; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[18], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[24], tt, tt);
    for (i = 0; i < 4; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[4], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[11], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[10], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[18], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[23], tt, tt);
    for (i = 0; i < 3; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[6], tt, tt);
    for (i = 0; i < 10; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[25], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[8], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[14], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[10], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[24], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[20], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[8], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[7], tt, tt);
    for (i = 0; i < 11; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[20], tt, tt);

    for (j = 0; j < 63; j++) {
        for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
        fpmul_mont(t[26], tt, tt);
    }
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[26], tt, a);
*/    

/*
#elif (NBITS_FIELD == 546)

    felm_t t[13], tt, t4;

    // Precomputed table
    fpsqr_mont(a, tt);
	fpsqr_mont(tt, t4);
    fpmul_mont(tt, a, t[0]);    //a^3
    fpmul_mont(t[0], tt, t[1]); //a^5
    fpmul_mont(t[1], tt, t[2]); //a^7
    fpmul_mont(t[2], t4, t[3]); //a^11
    fpmul_mont(t[3], tt, t[3]); //a^13
    for (i = 3; i <= 11; i++) fpmul_mont(t[i], tt, t[i+1]);

    fpcopy(t[0], tt);
    for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[4], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[8], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[7], tt, tt);
    for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[3], tt, tt);
    for (i=0; i < 3; i++) fpsqr_mont(tt, tt);
    fpmul_mont(a, tt, tt);
    for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[4], tt, tt);
    for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[7], tt, tt);
    for (i=0; i < 3; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[1], tt, tt);
    for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[9], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[9], tt, tt);
    for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[1], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[3], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[8], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[10], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[8], tt, tt);
    fpsqr_mont(tt, tt);
    fpmul_mont(a, tt, tt);
    for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[11], tt, tt);
    for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[11], tt, tt);
    for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[11], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[9], tt, tt);
    for (i=0; i < 2; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[0], tt, tt);
    for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[9], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[11], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[7], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[7], tt, tt);
    for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[4], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[9], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[1], tt, tt);
    for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[3], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[9], tt, tt);
    for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[3], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(a, tt, tt);
    for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(a, tt, tt);
    for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[4], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[6], tt, tt);
    for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[4], tt, tt);
    for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[9], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[11], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[9], tt, tt);
    for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[6], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[12], tt, tt);
    for(j=0; j < 53; j++){
        for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
        fpmul_mont(t[12], tt, tt);
    }
    fpsqr_mont(tt, tt);
    fpmul_mont(a, tt, a);
*/    

#ifdef P557_3_5_1

    felm_t t[15], tt;

    // Precomputed table
    fpsqr_mont(a, tt);
    fpmul_mont(tt, a, t[0]); // a^3
    for (i = 0; i <= 13; i++) fpmul_mont(t[i], tt, t[i+1]);
    
    fpcopy(t[6], tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[11], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[7], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[12], tt, tt);
    for (i=0; i < 3; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[8], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[6], tt, tt);
    for (i=0; i < 2; i++) fpsqr_mont(tt, tt);
    fpmul_mont(a, tt, tt);
    for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[8], tt, tt);
    for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[8], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[12], tt, tt);
    for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[8], tt, tt);
    for (i=0; i < 2; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[0], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(a, tt, tt);
    for (i=0; i < 10; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[10], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[11], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[13], tt, tt);
    for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[13], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[9], tt, tt);
    for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[1], tt, tt);
    for (i=0; i < 11; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[12], tt, tt);
    for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[3], tt, tt);
    for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[7], tt, tt);
    for (i=0; i < 2; i++) fpsqr_mont(tt, tt);
    fpmul_mont(a, tt, tt);
    for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[1], tt, tt);
    for (i=0; i < 10; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[12], tt, tt);
    for (i=0; i < 3; i++) fpsqr_mont(tt, tt);
    fpmul_mont(a, tt, tt);
    for (i=0; i < 10; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[4], tt, tt);
    for (i=0; i < 11; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[13], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[14], tt, tt);
    fpsqr_mont(tt, tt);
    fpmul_mont(a, tt, tt);
    for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[10], tt, tt);
    for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[11], tt, tt);
    for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[7], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[7], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[7], tt, tt);
    for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[0], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[1], tt, tt);
    for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[11], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[8], tt, tt);

    for(j=0; j < 55; j++){
        for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
        fpmul_mont(t[14], tt, tt);
    }
    fpsqr_mont(tt, tt);
    fpmul_mont(a, tt, a);
    
#endif    

#ifdef P443_3_5_1

    felm_t t[15], tt;

    // Precomputed table
    fpsqr_mont(a, tt);
    fpmul_mont(tt, a, t[0]);    //a^3
    for (i = 0; i <= 13; i++) fpmul_mont(t[i], tt, t[i+1]);

    fpcopy(t[3], tt);
    for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[6], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[4], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[9], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[9], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[9], tt, tt);
    for (i=0; i < 2; i++) fpsqr_mont(tt, tt);
    fpmul_mont(a, tt, tt);
    for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[11], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[13], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[8], tt, tt);
    for (i=0; i < 3; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[1], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[6], tt, tt);
    for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[1], tt, tt);
    for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[11], tt, tt);
    for (i=0; i < 3; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[1], tt, tt);
    for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[9], tt, tt);
    for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[10], tt, tt);
    for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[1], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[6], tt, tt);
    for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[1], tt, tt);
    for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[8], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[12], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[6], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[10], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[8], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[4], tt, tt);
    for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[0], tt, tt);
    for (i=0; i < 10; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[1], tt, tt);
    for (i=0; i < 3; i++) fpsqr_mont(tt, tt);
    fpmul_mont(a, tt, tt);
    for (i=0; i < 13; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[11], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[3], tt, tt);
    for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[13], tt, tt);
    for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[7], tt, tt);
    for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[4], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[8], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[12], tt, tt);
    for(j=0; j < 43; j++){
        for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
        fpmul_mont(t[14], tt, tt);
    }
    for (i=0; i < 3; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, a);
    
#endif    
/*
#elif (NBITS_FIELD == 434)

    felm_t t[14], tt;

    // Precomputed table
    fpsqr_mont(a, tt);
    fpmul_mont(tt, a, t[0]);    //a^3
    fpmul_mont(t[0], tt, t[1]);    //a^5
    fpmul_mont(t[1], tt, t[2]);    //a^7
    fpmul_mont(t[2], tt, t[3]);    //a^9
    fpmul_mont(t[3], tt, t[3]);    //a^11
    for (i = 3; i <= 12; i++) fpmul_mont(t[i], tt, t[i+1]);
        
    fpcopy(t[6], tt);
    for (i=0; i < 3; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[1], tt, tt);
    for (i=0; i < 10; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[13], tt, tt);
    for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[7], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[6], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[12], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[10], tt, tt);
    for (i=0; i < 3; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[1], tt, tt);
    for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[6], tt, tt);
    for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[3], tt, tt);
    for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[13], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[6], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[13], tt, tt);
    for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[4], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[4], tt, tt);
    for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[1], tt, tt);
    for (i=0; i < 3; i++) fpsqr_mont(tt, tt);
    fpmul_mont(a, tt, tt);
    for (i=0; i < 11; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[8], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[7], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[11], tt, tt);
    for (i=0; i < 2; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[0], tt, tt);
    for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[10], tt, tt);
    for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[0], tt, tt);
    for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[8], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[6], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[9], tt, tt);
    for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[1], tt, tt);
    for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[13], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[12], tt, tt);
    for (i=0; i < 2; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[0], tt, tt);
    for (i=0; i < 10; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[9], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[10], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[12], tt, tt);
    for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[9], tt, tt);
    for(j=0; j < 42; j++){
        for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
        fpmul_mont(t[13], tt, tt);
    }
    fpsqr_mont(tt, tt);
    fpmul_mont(a, tt, a);
*/

#ifdef P1013_3_5_1

    felm_t t[29], tt;

    // Precomputed table
    fpsqr_mont(a, tt);
    fpmul_mont(tt, a, t[0]);    //a^3
    for (i = 0; i <= 11; i++) fpmul_mont(t[i], tt, t[i+1]);
    fpmul_mont(t[12], tt, t[12]);    //a^29
    for (i = 12; i <= 25; i++) fpmul_mont(t[i], tt, t[i+1]);
    fpmul_mont(t[26], tt, t[26]);    //a^59
    fpmul_mont(t[26], tt, t[27]);    //a^61
    fpmul_mont(t[27], tt, t[28]);    //a^63

    fpcopy(t[0], tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[0], tt, tt);
    for (i=0; i < 12; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[3], tt, tt);
    for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[18], tt, tt);
    for (i=0; i < 2; i++) fpsqr_mont(tt, tt);
    fpmul_mont(a, tt, tt);
    for (i=0; i < 11; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[12], tt, tt);
    for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[16], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[4], tt, tt);
    for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[15], tt, tt);
    for (i=0; i < 3; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[1], tt, tt);
    for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[9], tt, tt);
    for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[16], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[27], tt, tt);
    for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[17], tt, tt);
    for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[23], tt, tt);
    for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[10], tt, tt);
    for (i=0; i < 2; i++) fpsqr_mont(tt, tt);
    fpmul_mont(a, tt, tt);
    for (i=0; i < 11; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[19], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[11], tt, tt);
    for (i=0; i < 10; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[15], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[27], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[8], tt, tt);
    for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[16], tt, tt);
    for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i=0; i < 11; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[19], tt, tt);
    for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[8], tt, tt);
    for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[22], tt, tt);
    for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[20], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[13], tt, tt);
    for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i=0; i < 10; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[19], tt, tt);
    for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[4], tt, tt);
    for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[22], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[13], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[24], tt, tt);
    for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[22], tt, tt);
    for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[19], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[23], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[22], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[7], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[9], tt, tt);
    for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[26], tt, tt);
    for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[19], tt, tt);
    for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i=0; i < 11; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[12], tt, tt);
    for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[23], tt, tt);
    for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[11], tt, tt);
    for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[25], tt, tt);
    for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[20], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[13], tt, tt);
    for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[16], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[17], tt, tt);
    for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[12], tt, tt);
    for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[12], tt, tt);
    for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[14], tt, tt);
    for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[4], tt, tt);
    for (i=0; i < 10; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[15], tt, tt);
    for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[18], tt, tt);
    for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[14], tt, tt);
    for (i=0; i < 10; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[18], tt, tt);
    for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[1], tt, tt);
    for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[24], tt, tt);
    for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[6], tt, tt);
    for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[8], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[4], tt, tt);
    for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i=0; i < 10; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[8], tt, tt);
    for (i=0; i < 12; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[21], tt, tt);
    
    for(j=0; j < 84; j++){
        for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
        fpmul_mont(t[28], tt, tt);
    }
    
    for (i=0; i < 3; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, a);

    
#endif

#ifdef P1005_3_7_17

    felm_t t[30], tt;

    // Precomputed table
    fpsqr_mont(a, tt);
    fpmul_mont(tt, a, t[0]);    //a^3
    for (i = 0; i <= 20; i++) fpmul_mont(t[i], tt, t[i+1]);
    fpmul_mont(t[21], tt, t[22]);
    fpmul_mont(t[22], tt, t[22]);
    for (i = 22; i <= 28; i++) fpmul_mont(t[i], tt, t[i+1]);

fpcopy(t[6], tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[18], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[18], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[8], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[3], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[7], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[19], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[0], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[11], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[10], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[16], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[12], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[1], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[2], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 3; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[1], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[25], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[11], tt, tt);
for (i=0; i < 10; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[7], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[18], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[4], tt, tt);
for (i=0; i < 10; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[12], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[4], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[24], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[2], tt, tt);
for (i=0; i < 10; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 2; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[0], tt, tt);
for (i=0; i < 10; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[24], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[4], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[24], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[17], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[15], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[25], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[20], tt, tt);
for (i=0; i < 3; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[2], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[9], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[8], tt, tt);
for (i=0; i < 2; i++) fpsqr_mont(tt, tt);
fpmul_mont(a, tt, tt);
for (i=0; i < 12; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[25], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[5], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[10], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[18], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[19], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[18], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[9], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[11], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[2], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[3], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[17], tt, tt);
for (i=0; i < 10; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[25], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[29], tt, tt);
fpsqr_mont(tt, tt);
fpmul_mont(a, tt, a);

#endif

#ifdef P1010_3_7_5

    felm_t t[28], tt, t4;

    // Precomputed table
    fpsqr_mont(a, tt);
    fpsqr_mont(tt, t4);
    fpmul_mont(tt, a, t[0]);
    for (i = 0; i <= 5; i++) fpmul_mont(t[i], tt, t[i+1]);
    fpmul_mont(t[6], t4, t[7]);
    for (i = 7; i <= 10; i++) fpmul_mont(t[i], tt, t[i+1]);
    fpmul_mont(t[11], t4, t[12]);
    fpmul_mont(t[12], tt, t[13]);
    fpmul_mont(t[13], t4, t[14]);
    for (i = 14; i <= 26; i++) fpmul_mont(t[i], tt, t[i+1]);

fpcopy(t[19], tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[12], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[8], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[6], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[3], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[18], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[1], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 3; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[2], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[4], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[0], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[0], tt, tt);
for (i=0; i < 11; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[3], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[4], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[8], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[8], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[12], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[10], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[0], tt, tt);
for (i=0; i < 10; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[17], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[5], tt, tt);
for (i=0; i < 12; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[24], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[5], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[2], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[24], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[20], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[5], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[2], tt, tt);
for (i=0; i < 10; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[16], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[10], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[20], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[11], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[2], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[24], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[7], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[7], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[25], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[18], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[15], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[25], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[20], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[24], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[9], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[12], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[8], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[1], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[11], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[7], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[20], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[17], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[11], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[8], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[4], tt, tt);
for (i=0; i < 13; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[2], tt, tt);
for (i=0; i < 11; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[15], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[12], tt, a);

#endif

#ifdef P1012_3_7_1

    felm_t t[29], tt;

    // Precomputed table
    fpsqr_mont(a, tt);
    fpmul_mont(tt, a, t[0]);    //a^3
    for (i = 0; i <= 3; i++) fpmul_mont(t[i], tt, t[i+1]);
    fpmul_mont(t[4], tt, t[5]);
    fpmul_mont(t[5], tt, t[5]);
    for (i = 5; i <= 8; i++) fpmul_mont(t[i], tt, t[i+1]);
    fpmul_mont(t[9], tt, t[10]);
    fpmul_mont(t[10], tt, t[10]);
    for (i = 10; i <= 27; i++) fpmul_mont(t[i], tt, t[i+1]);

fpcopy(t[2], tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[0], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[0], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[0], tt, tt);
for (i=0; i < 11; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[25], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[18], tt, tt);
fpsqr_mont(tt, tt);
fpmul_mont(a, tt, tt);
for (i=0; i < 12; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[17], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[3], tt, tt);
for (i=0; i < 11; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[10], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[5], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[7], tt, tt);
for (i=0; i < 3; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[0], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[4], tt, tt);
for (i=0; i < 11; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[20], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[12], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[19], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[7], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[3], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[3], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[18], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[6], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[1], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[7], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[13], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[3], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[9], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[6], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[27], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[20], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[15], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[3], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[4], tt, tt);
for (i=0; i < 5; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[1], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[2], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[8], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[8], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[16], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[18], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[19], tt, tt);
for (i=0; i < 10; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[11], tt, tt);
for (i=0; i < 2; i++) fpsqr_mont(tt, tt);
fpmul_mont(a, tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[1], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[19], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[26], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[2], tt, tt);
for (i=0; i < 10; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[21], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[11], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[9], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[18], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[24], tt, tt);
for (i=0; i < 10; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[16], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[12], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[6], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[10], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[0], tt, tt);
for (i=0; i < 13; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[14], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[5], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[22], tt, tt);
for (i=0; i < 4; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[5], tt, tt);
for (i=0; i < 8; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[23], tt, tt);
for (i=0; i < 9; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[25], tt, tt);
for (i=0; i < 7; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[20], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 6; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[28], tt, tt);
for (i=0; i < 2; i++) fpsqr_mont(tt, tt);
fpmul_mont(t[0], tt, a);
    
#endif    

}


void fp2inv_mont(f2elm_t a)
{// GF(p^2) inversion using Montgomery arithmetic, a = (a0-i*a1)/(a0^2+a1^2).
    f2elm_t t1;

    fpsqr_mont(a[0], t1[0]);                         // t10 = a0^2
    fpsqr_mont(a[1], t1[1]);                         // t11 = a1^2
    fpadd(t1[0], t1[1], t1[0]);                      // t10 = a0^2+a1^2
    fpinv_mont(t1[0]);                               // t10 = (a0^2+a1^2)^-1
    fpneg(a[1]);                                     // a = a0-i*a1
    fpmul_mont(a[0], t1[0], a[0]);
    fpmul_mont(a[1], t1[0], a[1]);                   // a = (a0-i*a1)*(a0^2+a1^2)^-1
}


void to_fp2mont(const f2elm_t a, f2elm_t mc)
{ // Conversion of a GF(p^2) element to Montgomery representation,
  // mc_i = a_i*R^2*R^(-1) = a_i*R in GF(p^2). 

    to_mont(a[0], mc[0]);
    to_mont(a[1], mc[1]);
}


void from_fp2mont(const f2elm_t ma, f2elm_t c)
{ // Conversion of a GF(p^2) element from Montgomery representation to standard representation,
  // c_i = ma_i*R^(-1) = a_i in GF(p^2).

    from_mont(ma[0], c[0]);
    from_mont(ma[1], c[1]);
}


__inline unsigned int mp_add(const digit_t* a, const digit_t* b, digit_t* c, const unsigned int nwords)
{ // Multiprecision addition, c = a+b, where lng(a) = lng(b) = nwords. Returns the carry bit.
    unsigned int i, carry = 0;
        
    for (i = 0; i < nwords; i++) {                      
        ADDC(carry, a[i], b[i], carry, c[i]);
    }

    return carry;
}


void mp_shiftleft(digit_t* x, unsigned int shift, const unsigned int nwords)
{
    unsigned int i, j = 0;

    while (shift > RADIX) {
        j += 1;
        shift -= RADIX;
    }

    for (i = 0; i < nwords-j; i++) 
        x[nwords-1-i] = x[nwords-1-i-j];
    for (i = nwords-j; i < nwords; i++) 
        x[nwords-1-i] = 0;
    if (shift != 0) {
        for (j = nwords-1; j > 0; j--) 
            SHIFTL(x[j], x[j-1], shift, x[j], RADIX);
        x[0] <<= shift;
    }
}


void mp_shiftr1(digit_t* x, const unsigned int nwords)
{ // Multiprecision right shift by one.
    unsigned int i;

    for (i = 0; i < nwords-1; i++) {
        SHIFTR(x[i+1], x[i], 1, x[i], RADIX);
    }
    x[nwords-1] >>= 1;
}


void mp_shiftl1(digit_t* x, const unsigned int nwords)
{ // Multiprecision left shift by one.
    int i;

    for (i = nwords-1; i > 0; i--) {
        SHIFTL(x[i], x[i-1], 1, x[i], RADIX);
    }
    x[0] <<= 1;
}
