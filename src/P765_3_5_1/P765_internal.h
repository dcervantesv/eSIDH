/********************************************************************************************
* SIDH: an efficient supersingular isogeny cryptography library
*
* Abstract: internal header file for P765
*********************************************************************************************/  

#ifndef __P765_INTERNAL_H__
#define __P765_INTERNAL_H__

#include "../config.h"
 

#if (TARGET == TARGET_AMD64)
    #define NWORDS_FIELD    12              // Number of words of a 765-bit field element
    #define p765_ZERO_WORDS 6               // Number of "0" digits in the least significant part of p765 + 1     
#elif (TARGET == TARGET_x86)
    #define NWORDS_FIELD    24 
    #define p765_ZERO_WORDS 12
#elif (TARGET == TARGET_ARM)
    #define NWORDS_FIELD    24
    #define p765_ZERO_WORDS 12
#elif (TARGET == TARGET_ARM64)
    #define NWORDS_FIELD    12
    #define p765_ZERO_WORDS 6
#endif 
    

// Basic constants

#define NBITS_FIELD             765  
#define MAXBITS_FIELD           768                
#define MAXWORDS_FIELD          12//((MAXBITS_FIELD+RADIX-1)/RADIX)     // Max. number of words to represent field elements
#define NWORDS64_FIELD          ((NBITS_FIELD+63)/64)               // Number of 64-bit words of a 765-bit field element 
#define NBITS_ORDER             384
#define NWORDS_ORDER            ((NBITS_ORDER+RADIX-1)/RADIX)       // Number of words of oA and oB, where oA and oB are the subgroup orders of Alice and Bob, resp.
#define NWORDS64_ORDER          ((NBITS_ORDER+63)/64)               // Number of 64-bit words of a 384-bit element 
#define MAXBITS_ORDER           NBITS_ORDER                         
#define MAXWORDS_ORDER          ((MAXBITS_ORDER+RADIX-1)/RADIX)     // Max. number of words to represent elements in [1, oA-1] or [1, oB].
#define BOB_PRIMES              2  // Number of primes for Bob's part
#define ALICE                   0
#define BOB                     1     //BOB3
#define BOB1                    2
#define OALICE_BITS             372  
#define OBOB_BITS               189   //BOB3
#define OBOB1_BITS              189 
#define OBOB_EXPON              189   //no se usa
#define MASK_ALICE              0x0F  
#define MASK_BOB                0x17  
#define MASK_BOB1               0x10  
#define PRIME                   p765
#define PARAM_A                 0  
#define PARAM_C                 1
// Fixed parameters for isogeny tree computation
#define MAX_INT_POINTS_ALICE    185      
#define MAX_INT_POINTS_BOB1     30 
#define MAX_INT_POINTS_BOB2     30 
#define MAX_Alice               186
#define MAX_Bob1                119   //BOB3
#define MAX_Bob2                81
#define Index_Bob1              3    // private key index for the Bob5 ladder 
#define KERNEL_POINTS1          1
#define KERNEL_POINTS2          2
//extras
#define MSG_BYTES               32
#define SECRETKEY_A_BYTES       47
#define SECRETKEY_B_BYTES       48
#define FP2_ENCODED_BYTES       2*((NBITS_FIELD + 7) / 8)



// SIDH's basic element definitions and point representations

typedef digit_t felm_t[NWORDS_FIELD];                                 // Datatype for representing 765-bit field elements (768-bit max.)
typedef digit_t dfelm_t[2*NWORDS_FIELD];                              // Datatype for representing double-precision 2x765-bit field elements (2x768-bit max.) 
typedef felm_t  f2elm_t[2];                                           // Datatype for representing quadratic extension field elements GF(p765^2)
        
typedef struct { f2elm_t X; f2elm_t Z; } point_proj;                  // Point representation in projective XZ Montgomery coordinates.
typedef point_proj point_proj_t[1]; 



/**************** Function prototypes ****************/
/************* Multiprecision functions **************/ 

// Copy wordsize digits, c = a, where lng(a) = nwords
void copy_words(const digit_t* a, digit_t* c, const unsigned int nwords);

// Multiprecision addition, c = a+b, where lng(a) = lng(b) = nwords. Returns the carry bit 
unsigned int mp_add(const digit_t* a, const digit_t* b, digit_t* c, const unsigned int nwords);

// 765-bit multiprecision addition, c = a+b
void mp_add765(const digit_t* a, const digit_t* b, digit_t* c);
void mp_add765_asm(const digit_t* a, const digit_t* b, digit_t* c);

// Multiprecision subtraction, c = a-b, where lng(a) = lng(b) = nwords. Returns the borrow bit 
unsigned int mp_sub(const digit_t* a, const digit_t* b, digit_t* c, const unsigned int nwords);
digit_t mp_sub765x2_asm(const digit_t* a, const digit_t* b, digit_t* c);

// Double 2x765-bit multiprecision subtraction, c = c-a-b, where c > a and c > b
void mp_dblsub765x2_asm(const digit_t* a, const digit_t* b, digit_t* c);

// Multiprecision left shift
void mp_shiftleft(digit_t* x, unsigned int shift, const unsigned int nwords);

// Multiprecision right shift by one
void mp_shiftr1(digit_t* x, const unsigned int nwords);

// Multiprecision left right shift by one    
void mp_shiftl1(digit_t* x, const unsigned int nwords);

// Digit multiplication, digit * digit -> 2-digit result
void digit_x_digit(const digit_t a, const digit_t b, digit_t* c);

// Multiprecision comba multiply, c = a*b, where lng(a) = lng(b) = nwords.
void mp_mul(const digit_t* a, const digit_t* b, digit_t* c, const unsigned int nwords);

/************ Field arithmetic functions *************/

// Copy of a field element, c = a
void fpcopy765(const digit_t* a, digit_t* c);

// Zeroing a field element, a = 0
void fpzero765(digit_t* a);

// Non constant-time comparison of two field elements. If a = b return TRUE, otherwise, return FALSE
bool fpequal765_non_constant_time(const digit_t* a, const digit_t* b); 

// Modular addition, c = a+b mod p765
extern void fpadd765(const digit_t* a, const digit_t* b, digit_t* c);
extern void fpadd765_asm(const digit_t* a, const digit_t* b, digit_t* c);

// Modular subtraction, c = a-b mod p765
extern void fpsub765(const digit_t* a, const digit_t* b, digit_t* c);
extern void fpsub765_asm(const digit_t* a, const digit_t* b, digit_t* c);

// Modular negation, a = -a mod p765        
extern void fpneg765(digit_t* a);  

// Modular division by two, c = a/2 mod p765.
void fpdiv2_765(const digit_t* a, digit_t* c);

// Modular correction to reduce field element a in [0, 2*p765-1] to [0, p765-1].
void fpcorrection765(digit_t* a);

// 765-bit Montgomery reduction, c = a mod p
void rdc_mont(const digit_t* a, digit_t* c);
            
// Field multiplication using Montgomery arithmetic, c = a*b*R^-1 mod p765, where R=2^768
void fpmul765_mont(const digit_t* a, const digit_t* b, digit_t* c);
void mul765_asm(const digit_t* a, const digit_t* b, digit_t* c);
void rdc765_asm(const digit_t* ma, digit_t* mc);
   
// Field squaring using Montgomery arithmetic, c = a*b*R^-1 mod p765, where R=2^768
void fpsqr765_mont(const digit_t* ma, digit_t* mc);
void sqr765_asm(const digit_t* a, digit_t* b);

// Conversion to Montgomery representation
void to_mont(const digit_t* a, digit_t* mc);
    
// Conversion from Montgomery representation to standard representation
void from_mont(const digit_t* ma, digit_t* c);

// Field inversion, a = a^-1 in GF(p765)
void fpinv765_mont(digit_t* a);

// Field inversion, a = a^-1 in GF(p765) using the binary GCD 
void fpinv765_mont_bingcd(digit_t* a);

// Chain to compute (p765-3)/4 using Montgomery arithmetic
void fpinv765_chain_mont(digit_t* a);

/************ GF(p^2) arithmetic functions *************/
    
// Copy of a GF(p765^2) element, c = a
void fp2copy765(const f2elm_t a, f2elm_t c);

// Zeroing a GF(p765^2) element, a = 0
void fp2zero765(f2elm_t a);

// GF(p765^2) negation, a = -a in GF(p765^2)
void fp2neg765(f2elm_t a);

// GF(p765^2) addition, c = a+b in GF(p765^2)
extern void fp2add765(const f2elm_t a, const f2elm_t b, f2elm_t c);           

// GF(p765^2) subtraction, c = a-b in GF(p765^2)
extern void fp2sub765(const f2elm_t a, const f2elm_t b, f2elm_t c); 

// GF(p765^2) division by two, c = a/2  in GF(p765^2) 
void fp2div2_765(const f2elm_t a, f2elm_t c);

// Modular correction, a = a in GF(p765^2)
void fp2correction765(f2elm_t a);
            
// GF(p765^2) squaring using Montgomery arithmetic, c = a^2 in GF(p765^2)
void fp2sqr765_mont(const f2elm_t a, f2elm_t c);
 
// GF(p765^2) multiplication using Montgomery arithmetic, c = a*b in GF(p765^2)
void fp2mul765_mont(const f2elm_t a, const f2elm_t b, f2elm_t c);
    
// Conversion of a GF(p765^2) element to Montgomery representation
void to_fp2mont(const f2elm_t a, f2elm_t mc);

// Conversion of a GF(p765^2) element from Montgomery representation to standard representation
void from_fp2mont(const f2elm_t ma, f2elm_t c);

// GF(p765^2) inversion using Montgomery arithmetic, a = (a0-i*a1)/(a0^2+a1^2)
void fp2inv765_mont(f2elm_t a);

// GF(p765^2) inversion, a = (a0-i*a1)/(a0^2+a1^2), GF(p765) inversion done using the binary GCD 
void fp2inv765_mont_bingcd(f2elm_t a);

// n-way Montgomery inversion
void mont_n_way_inv(const f2elm_t* vec, const int n, f2elm_t* out);

/************ Elliptic curve and isogeny functions *************/

// Computes the j-invariant of a Montgomery curve with projective constant.
void j_inv(const f2elm_t A, const f2elm_t C, f2elm_t jinv);

// Simultaneous doubling and differential addition.
void xDBLADD(point_proj_t P, point_proj_t Q, const f2elm_t xPQ, const f2elm_t A24);

// Doubling of a Montgomery point in projective coordinates (X:Z).
void xDBL(const point_proj_t P, point_proj_t Q, const f2elm_t A24plus, const f2elm_t C24);

// Computes [2^e](X:Z) on Montgomery curve with projective constant via e repeated doublings.
void xDBLe(const point_proj_t P, point_proj_t Q, const f2elm_t A24plus, const f2elm_t C24, const int e);

// Differential addition.
void xADD(point_proj_t P, const point_proj_t Q, const f2elm_t xPQ);

// Computes the corresponding 4-isogeny of a projective Montgomery point (X4:Z4) of order 4.
void get_4_isog(const point_proj_t P, f2elm_t A24plus, f2elm_t C24, f2elm_t* coeff);

// Evaluates the isogeny at the point (X:Z) in the domain of the isogeny.
void eval_4_isog(point_proj_t P, f2elm_t* coeff);

// Tripling of a Montgomery point in projective coordinates (X:Z).
void xTPL(const point_proj_t P, point_proj_t Q, const f2elm_t A24minus, const f2elm_t A24plus);

// Computes [3^e](X:Z) on Montgomery curve with projective constant via e repeated triplings.
void xTPLe(const point_proj_t P, point_proj_t Q, const f2elm_t A24minus, const f2elm_t A24plus, const int e);

// Computes the corresponding 3-isogeny of a projective Montgomery point (X3:Z3) of order 3.
void get_3_isog(const point_proj_t P, f2elm_t A24minus, f2elm_t A24plus, f2elm_t* coeff);

// Computes the 3-isogeny R=phi(X:Z), given projective point (X3:Z3) of order 3 on a Montgomery curve and a point P with coefficients given in coeff.
void eval_3_isog(point_proj_t Q, const f2elm_t* coeff);

//Newfunctions for [5], get 5-isogeny and eval 5-isogeny
void xQPL(const point_proj_t P, point_proj_t Q, const f2elm_t A24p, const f2elm_t A24m,  const f2elm_t C24);
void xQPLe(const point_proj_t P, point_proj_t Q, const f2elm_t A24p, const f2elm_t A24m,  const f2elm_t C24, int e);
void eval_5_isog(const point_proj_t *KPs, point_proj_t Q);
void get_5_isog(const point_proj_t P,  point_proj_t *KPs,  f2elm_t A24plus, f2elm_t A24minus, f2elm_t C24);
void yDBL(const point_proj_t P, point_proj_t Q, const f2elm_t A24plus, const f2elm_t C24);

// //Newfunctions for [7], get 7-isogeny and eval 7-isogeny
// void xSPL(const point_proj_t P, point_proj_t Q, const f2elm_t A24p, const f2elm_t A24m,  const f2elm_t C24);
// void xSPLe(const point_proj_t P, point_proj_t Q, const f2elm_t A24p, const f2elm_t A24m,  const f2elm_t C24, int e);
// void eval_7_isog(const point_proj_t *KPs, point_proj_t Q);
// void get_7_isog(const point_proj_t P,  point_proj_t *KPs,  f2elm_t A24plus, f2elm_t A24minus, f2elm_t C24);

// 3-way simultaneous inversion
void inv_3_way(f2elm_t z1, f2elm_t z2, f2elm_t z3);

// Given the x-coordinates of P, Q, and R, returns the value A corresponding to the Montgomery curve E_A: y^2=x^3+A*x^2+x such that R=Q-P on E_A.
void get_A(const f2elm_t xP, const f2elm_t xQ, const f2elm_t xR, f2elm_t A);


#endif
