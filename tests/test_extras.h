/********************************************************************************************
* SIDH: an efficient supersingular isogeny cryptography library
*
* Abstract: utility header file for tests
*********************************************************************************************/  

#ifndef __TEST_EXTRAS_H__
#define __TEST_EXTRAS_H__
    
#include "config.h"

#define PASSED    0
#define FAILED    1


#if (TARGET == TARGET_ARM || TARGET == TARGET_ARM64)
    #define print_unit printf("nsec");
#else
    #define print_unit printf("cycles");
#endif

    
// Access system counter for benchmarking
int64_t cpucycles(void);

// Comparing "nword" elements, a=b? : (1) a!=b, (0) a=b
int compare_words(digit_t* a, digit_t* b, unsigned int nwords);

// Generating a pseudo-random field element in [0, p503-1] 
void fprandom503_test(digit_t* a);

// Generating a pseudo-random element in GF(p503^2)
void fp2random503_test(digit_t* a);

// Generating a pseudo-random field element in [0, p751-1] 
void fprandom751_test(digit_t* a);

// Generating a pseudo-random element in GF(p751^2)
void fp2random751_test(digit_t* a);

// Generating a pseudo-random field element in [0, p765-1] 
void fprandom765_test(digit_t* a);

// Generating a pseudo-random element in GF(p765^2)
void fp2random765_test(digit_t* a);

// Generating a pseudo-random field element in [0, p765-1] 
void fprandom765_7_5_test(digit_t* a);

// Generating a pseudo-random element in GF(p765^2)
void fp2random765_7_5_test(digit_t* a);

// Generating a pseudo-random field element in [0, p766-1] 
void fprandom766_test(digit_t* a);

// Generating a pseudo-random element in GF(p766^2)
void fp2random766_test(digit_t* a);

// Generating a pseudo-random field element in [0, p764-1] 
void fprandom764_test(digit_t* a);

// Generating a pseudo-random element in GF(p764^2)
void fp2random764_test(digit_t* a);

// Generating a pseudo-random element in GF(p509^2)
void fprandom509_test(digit_t* a);

// Generating a pseudo-random element in GF(p509^2)
void fp2random509_test(digit_t* a);

// Generating a pseudo-random element in GF(p546^2)
void fprandom546_test(digit_t* a);

// Generating a pseudo-random element in GF(p546^2)
void fp2random546_test(digit_t* a);

// Generating a pseudo-random element in GF(p557^2)
void fprandom557_test(digit_t* a);

// Generating a pseudo-random element in GF(p557^2)
void fp2random557_test(digit_t* a);

// Generating a pseudo-random element in GF(p443^2)
void fprandom443_test(digit_t* a);

// Generating a pseudo-random element in GF(p443^2)
void fp2random443_test(digit_t* a);

// Generating a pseudo-random element in GF(p434^2)
void fprandom434_test(digit_t* a);

// Generating a pseudo-random element in GF(p434^2)
void fp2random434_test(digit_t* a);

// Generating a pseudo-random element in GF(p1013^2)
void fprandom1013_test(digit_t* a);

// Generating a pseudo-random element in GF(p1013^2)
void fp2random1013_test(digit_t* a);

// Generating a pseudo-random element in GF(p1005^2)
void fprandom1005_test(digit_t* a);

// Generating a pseudo-random element in GF(p1005^2)
void fp2random1005_test(digit_t* a);

// Generating a pseudo-random element in GF(p1010^2)
void fprandom1010_test(digit_t* a);

// Generating a pseudo-random element in GF(p1010^2)
void fp2random1010_test(digit_t* a);

// Generating a pseudo-random element in GF(p1012^2)
void fprandom1012_test(digit_t* a);

// Generating a pseudo-random element in GF(p1012^2)
void fp2random1012_test(digit_t* a);

// Generating a pseudo-random element in GF(p507^2)
void fprandom507_test(digit_t* a);

// Generating a pseudo-random element in GF(p507^2)
void fp2random507_test(digit_t* a);

// Generating a pseudo-random element in GF(p508^2)
void fprandom508_test(digit_t* a);

// Generating a pseudo-random element in GF(p508^2)
void fp2random508_test(digit_t* a);


#endif
