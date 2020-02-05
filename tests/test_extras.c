/********************************************************************************************
* SIDH: an efficient supersingular isogeny cryptography library
*
* Abstract: utility functions for testing and benchmarking
*********************************************************************************************/

#include "test_extras.h"
#if (OS_TARGET == OS_WIN)
    #include <intrin.h>
#endif
#if (OS_TARGET == OS_LINUX) && (TARGET == TARGET_ARM || TARGET == TARGET_ARM64)
    #include <time.h>
#endif
#include <stdlib.h>


int64_t cpucycles(void)
{ // Access system counter for benchmarking
#if (OS_TARGET == OS_WIN) && (TARGET == TARGET_AMD64 || TARGET == TARGET_x86)
    return __rdtsc();
#elif (OS_TARGET == OS_WIN) && (TARGET == TARGET_ARM)
    return __rdpmccntr64();
#elif (OS_TARGET == OS_LINUX) && (TARGET == TARGET_AMD64 || TARGET == TARGET_x86)
    unsigned int hi, lo;

    asm volatile ("rdtsc\n\t" : "=a" (lo), "=d"(hi));
    return ((int64_t)lo) | (((int64_t)hi) << 32);
#elif (OS_TARGET == OS_LINUX) && (TARGET == TARGET_ARM || TARGET == TARGET_ARM64)
    struct timespec time;

    clock_gettime(CLOCK_REALTIME, &time);
    return (int64_t)(time.tv_sec*1e9 + time.tv_nsec);
#else
    return 0;            
#endif
}


int compare_words(digit_t* a, digit_t* b, unsigned int nwords)
{ // Comparing "nword" elements, a=b? : (1) a>b, (0) a=b, (-1) a<b
  // SECURITY NOTE: this function does not have constant-time execution. TO BE USED FOR TESTING ONLY.
    int i;

    for (i = nwords-1; i >= 0; i--)
    {
        if (a[i] > b[i]) return 1;
        else if (a[i] < b[i]) return -1;
    }

    return 0; 
}


static void sub_test(digit_t* a, digit_t* b, digit_t* c, unsigned int nwords)
{ // Subtraction without borrow, c = a-b where a>b
  // SECURITY NOTE: this function does not have constant-time execution. It is for TESTING ONLY.     
    unsigned int i;
    digit_t res, carry, borrow = 0;
  
    for (i = 0; i < nwords; i++)
    {
        res = a[i] - b[i];
        carry = (a[i] < b[i]);
        c[i] = res - borrow;
        borrow = carry || (res < borrow);
    } 
}

#ifdef P751

extern uint64_t p751[12];

void fprandom751_test(digit_t* a)
{ // Generating a pseudo-random field element in [0, p751-1] 
  // SECURITY NOTE: distribution is not fully uniform. TO BE USED FOR TESTING ONLY.
    unsigned int i, diff = 768-751, nwords = NBITS_TO_NWORDS(751);
    unsigned char* string = NULL;

    string = (unsigned char*)a;
    for (i = 0; i < sizeof(digit_t)*nwords; i++) {
        *(string + i) = (unsigned char)rand();              // Obtain 768-bit number
    }
    a[nwords-1] &= (((digit_t)(-1) << diff) >> diff);

    while (compare_words((digit_t*)p751, a, nwords) < 1) {  // Force it to [0, modulus-1]
        sub_test(a, (digit_t*)p751, a, nwords);
    }
}


void fp2random751_test(digit_t* a)
{ // Generating a pseudo-random element in GF(p751^2) 
  // SECURITY NOTE: distribution is not fully uniform. TO BE USED FOR TESTING ONLY.

    fprandom751_test(a);
    fprandom751_test(a+NBITS_TO_NWORDS(751));
}

#endif

#ifdef P765

extern uint64_t p765[12];

void fprandom765_test(digit_t* a)
{ // Generating a pseudo-random field element in [0, p751-1] 
  // SECURITY NOTE: distribution is not fully uniform. TO BE USED FOR TESTING ONLY.
    unsigned int i, diff = 768-765, nwords = NBITS_TO_NWORDS(765);
    unsigned char* string = NULL;

    string = (unsigned char*)a;
    for (i = 0; i < sizeof(digit_t)*nwords; i++) {
        *(string + i) = (unsigned char)rand();              // Obtain 768-bit number
    }
    a[nwords-1] &= (((digit_t)(-1) << diff) >> diff);

    while (compare_words((digit_t*)p765, a, nwords) < 1) {  // Force it to [0, modulus-1]
        sub_test(a, (digit_t*)p765, a, nwords);
    }
}

void fp2random765_test(digit_t* a)
{ // Generating a pseudo-random element in GF(p751^2) 
  // SECURITY NOTE: distribution is not fully uniform. TO BE USED FOR TESTING ONLY.

    fprandom765_test(a);
    fprandom765_test(a+NBITS_TO_NWORDS(765));
}

#endif

#ifdef P766

extern uint64_t p766[12];

void fprandom766_test(digit_t* a)
{ // Generating a pseudo-random field element in [0, p751-1] 
  // SECURITY NOTE: distribution is not fully uniform. TO BE USED FOR TESTING ONLY.
    unsigned int i, diff = 768-766, nwords = NBITS_TO_NWORDS(766);
    unsigned char* string = NULL;

    string = (unsigned char*)a;
    for (i = 0; i < sizeof(digit_t)*nwords; i++) {
        *(string + i) = (unsigned char)rand();              // Obtain 768-bit number
    }
    a[nwords-1] &= (((digit_t)(-1) << diff) >> diff);

    while (compare_words((digit_t*)p766, a, nwords) < 1) {  // Force it to [0, modulus-1]
        sub_test(a, (digit_t*)p766, a, nwords);
    }
}


void fp2random766_test(digit_t* a)
{ // Generating a pseudo-random element in GF(p751^2) 
  // SECURITY NOTE: distribution is not fully uniform. TO BE USED FOR TESTING ONLY.

    fprandom766_test(a);
    fprandom766_test(a+NBITS_TO_NWORDS(766));
}

#endif

#ifdef P764

extern uint64_t p764[12];

void fprandom764_test(digit_t* a)
{ // Generating a pseudo-random field element in [0, p751-1] 
  // SECURITY NOTE: distribution is not fully uniform. TO BE USED FOR TESTING ONLY.
    unsigned int i, diff = 768-764, nwords = NBITS_TO_NWORDS(764);
    unsigned char* string = NULL;

    string = (unsigned char*)a;
    for (i = 0; i < sizeof(digit_t)*nwords; i++) {
        *(string + i) = (unsigned char)rand();              // Obtain 768-bit number
    }
    a[nwords-1] &= (((digit_t)(-1) << diff) >> diff);

    while (compare_words((digit_t*)p764, a, nwords) < 1) {  // Force it to [0, modulus-1]
        sub_test(a, (digit_t*)p764, a, nwords);
    }
}

void fp2random764_test(digit_t* a)
{ // Generating a pseudo-random element in GF(p751^2) 
  // SECURITY NOTE: distribution is not fully uniform. TO BE USED FOR TESTING ONLY.

    fprandom764_test(a);
    fprandom764_test(a+NBITS_TO_NWORDS(764));
}

#endif

#ifdef P509

extern uint64_t p509[8];

void fprandom509_test(digit_t* a)
{ // Generating a pseudo-random field element in [0, p751-1] 
  // SECURITY NOTE: distribution is not fully uniform. TO BE USED FOR TESTING ONLY.
    unsigned int i, diff = 512-509, nwords = NBITS_TO_NWORDS(509);
    unsigned char* string = NULL;

    string = (unsigned char*)a;
    for (i = 0; i < sizeof(digit_t)*nwords; i++) {
        *(string + i) = (unsigned char)rand();              // Obtain 768-bit number
    }
    a[nwords-1] &= (((digit_t)(-1) << diff) >> diff);

    while (compare_words((digit_t*)p509, a, nwords) < 1) {  // Force it to [0, modulus-1]
        sub_test(a, (digit_t*)p509, a, nwords);
    }
}

void fp2random509_test(digit_t* a)
{ // Generating a pseudo-random element in GF(p751^2) 
  // SECURITY NOTE: distribution is not fully uniform. TO BE USED FOR TESTING ONLY.

    fprandom509_test(a);
    fprandom509_test(a+NBITS_TO_NWORDS(509));
}

#endif

#ifdef P546

extern uint64_t p546[9];

void fprandom546_test(digit_t* a)
{ // Generating a pseudo-random field element in [0, p751-1] 
  // SECURITY NOTE: distribution is not fully uniform. TO BE USED FOR TESTING ONLY.
    unsigned int i, diff = 576-546, nwords = NBITS_TO_NWORDS(546);
    unsigned char* string = NULL;

    string = (unsigned char*)a;
    for (i = 0; i < sizeof(digit_t)*nwords; i++) {
        *(string + i) = (unsigned char)rand();              // Obtain 768-bit number
    }
    a[nwords-1] &= (((digit_t)(-1) << diff) >> diff);

    while (compare_words((digit_t*)p546, a, nwords) < 1) {  // Force it to [0, modulus-1]
        sub_test(a, (digit_t*)p546, a, nwords);
    }
}

void fp2random546_test(digit_t* a)
{ // Generating a pseudo-random element in GF(p751^2) 
  // SECURITY NOTE: distribution is not fully uniform. TO BE USED FOR TESTING ONLY.

    fprandom546_test(a);
    fprandom546_test(a+NBITS_TO_NWORDS(546));
}

#endif

#ifdef P557

extern uint64_t p557[9];

void fprandom557_test(digit_t* a)
{ // Generating a pseudo-random field element in [0, p751-1] 
  // SECURITY NOTE: distribution is not fully uniform. TO BE USED FOR TESTING ONLY.
    unsigned int i, diff = 576-557, nwords = NBITS_TO_NWORDS(557);
    unsigned char* string = NULL;

    string = (unsigned char*)a;
    for (i = 0; i < sizeof(digit_t)*nwords; i++) {
        *(string + i) = (unsigned char)rand();              // Obtain 768-bit number
    }
    a[nwords-1] &= (((digit_t)(-1) << diff) >> diff);

    while (compare_words((digit_t*)p557, a, nwords) < 1) {  // Force it to [0, modulus-1]
        sub_test(a, (digit_t*)p557, a, nwords);
    }
}

void fp2random557_test(digit_t* a)
{ // Generating a pseudo-random element in GF(p751^2) 
  // SECURITY NOTE: distribution is not fully uniform. TO BE USED FOR TESTING ONLY.

    fprandom557_test(a);
    fprandom557_test(a+NBITS_TO_NWORDS(557));
}

#endif

#ifdef P443

extern uint64_t p443[7];

void fprandom443_test(digit_t* a)
{ // Generating a pseudo-random field element in [0, p751-1] 
  // SECURITY NOTE: distribution is not fully uniform. TO BE USED FOR TESTING ONLY.
    unsigned int i, diff = 448-443, nwords = NBITS_TO_NWORDS(443);
    unsigned char* string = NULL;

    string = (unsigned char*)a;
    for (i = 0; i < sizeof(digit_t)*nwords; i++) {
        *(string + i) = (unsigned char)rand();              // Obtain 768-bit number
    }
    a[nwords-1] &= (((digit_t)(-1) << diff) >> diff);

    while (compare_words((digit_t*)p443, a, nwords) < 1) {  // Force it to [0, modulus-1]
        sub_test(a, (digit_t*)p443, a, nwords);
    }
}

void fp2random443_test(digit_t* a)
{ // Generating a pseudo-random element in GF(p751^2) 
  // SECURITY NOTE: distribution is not fully uniform. TO BE USED FOR TESTING ONLY.

    fprandom443_test(a);
    fprandom443_test(a+NBITS_TO_NWORDS(443));
}

#endif

#ifdef P434

extern uint64_t p434[7];

void fprandom434_test(digit_t* a)
{ // Generating a pseudo-random field element in [0, p751-1] 
  // SECURITY NOTE: distribution is not fully uniform. TO BE USED FOR TESTING ONLY.
    unsigned int i, diff = 448-434, nwords = NBITS_TO_NWORDS(434);
    unsigned char* string = NULL;

    string = (unsigned char*)a;
    for (i = 0; i < sizeof(digit_t)*nwords; i++) {
        *(string + i) = (unsigned char)rand();              // Obtain 768-bit number
    }
    a[nwords-1] &= (((digit_t)(-1) << diff) >> diff);

    while (compare_words((digit_t*)p434, a, nwords) < 1) {  // Force it to [0, modulus-1]
        sub_test(a, (digit_t*)p434, a, nwords);
    }
}

void fp2random434_test(digit_t* a)
{ // Generating a pseudo-random element in GF(p751^2) 
  // SECURITY NOTE: distribution is not fully uniform. TO BE USED FOR TESTING ONLY.

    fprandom434_test(a);
    fprandom434_test(a+NBITS_TO_NWORDS(434));
}

#endif

#ifdef P1013

extern uint64_t p1013[16];

void fprandom1013_test(digit_t* a)
{ // Generating a pseudo-random field element in [0, p751-1] 
  // SECURITY NOTE: distribution is not fully uniform. TO BE USED FOR TESTING ONLY.
    unsigned int i, diff = 1024-1013, nwords = NBITS_TO_NWORDS(1013);
    unsigned char* string = NULL;

    string = (unsigned char*)a;
    for (i = 0; i < sizeof(digit_t)*nwords; i++) {
        *(string + i) = (unsigned char)rand();              // Obtain 768-bit number
    }
    a[nwords-1] &= (((digit_t)(-1) << diff) >> diff);

    while (compare_words((digit_t*)p1013, a, nwords) < 1) {  // Force it to [0, modulus-1]
        sub_test(a, (digit_t*)p1013, a, nwords);
    }
}

void fp2random1013_test(digit_t* a)
{ // Generating a pseudo-random element in GF(p751^2) 
  // SECURITY NOTE: distribution is not fully uniform. TO BE USED FOR TESTING ONLY.

    fprandom1013_test(a);
    fprandom1013_test(a+NBITS_TO_NWORDS(1013));
}

#endif

#ifdef P1005

extern uint64_t p1005[16];

void fprandom1005_test(digit_t* a)
{ // Generating a pseudo-random field element in [0, p751-1] 
  // SECURITY NOTE: distribution is not fully uniform. TO BE USED FOR TESTING ONLY.
    unsigned int i, diff = 1024-1005, nwords = NBITS_TO_NWORDS(1005);
    unsigned char* string = NULL;


    string = (unsigned char*)a;
    for (i = 0; i < sizeof(digit_t)*nwords; i++) {
        *(string + i) = (unsigned char)rand();              // Obtain 768-bit number
    }
    a[nwords-1] &= (((digit_t)(-1) << diff) >> diff);

    while (compare_words((digit_t*)p1005, a, nwords) < 1) {  // Force it to [0, modulus-1]
        sub_test(a, (digit_t*)p1005, a, nwords);
    }
}

void fp2random1005_test(digit_t* a)
{ // Generating a pseudo-random element in GF(p751^2) 
  // SECURITY NOTE: distribution is not fully uniform. TO BE USED FOR TESTING ONLY.

    fprandom1005_test(a);
    fprandom1005_test(a+NBITS_TO_NWORDS(1005));
}

#endif


#ifdef P1010

extern uint64_t p1010[16];

void fprandom1010_test(digit_t* a)
{ // Generating a pseudo-random field element in [0, p751-1] 
  // SECURITY NOTE: distribution is not fully uniform. TO BE USED FOR TESTING ONLY.
    unsigned int i, diff = 1024-1010, nwords = NBITS_TO_NWORDS(1010);
    unsigned char* string = NULL;


    string = (unsigned char*)a;
    for (i = 0; i < sizeof(digit_t)*nwords; i++) {
        *(string + i) = (unsigned char)rand();              // Obtain 768-bit number
    }
    a[nwords-1] &= (((digit_t)(-1) << diff) >> diff);

    while (compare_words((digit_t*)p1010, a, nwords) < 1) {  // Force it to [0, modulus-1]
        sub_test(a, (digit_t*)p1010, a, nwords);
    }
}

void fp2random1010_test(digit_t* a)
{ // Generating a pseudo-random element in GF(p751^2) 
  // SECURITY NOTE: distribution is not fully uniform. TO BE USED FOR TESTING ONLY.

    fprandom1010_test(a);
    fprandom1010_test(a+NBITS_TO_NWORDS(1010));
}

#endif

#ifdef P1012

extern uint64_t p1012[16];

void fprandom1012_test(digit_t* a)
{ // Generating a pseudo-random field element in [0, p751-1] 
  // SECURITY NOTE: distribution is not fully uniform. TO BE USED FOR TESTING ONLY.
    unsigned int i, diff = 1024-1012, nwords = NBITS_TO_NWORDS(1012);
    unsigned char* string = NULL;


    string = (unsigned char*)a;
    for (i = 0; i < sizeof(digit_t)*nwords; i++) {
        *(string + i) = (unsigned char)rand();              // Obtain 768-bit number
    }
    a[nwords-1] &= (((digit_t)(-1) << diff) >> diff);

    while (compare_words((digit_t*)p1012, a, nwords) < 1) {  // Force it to [0, modulus-1]
        sub_test(a, (digit_t*)p1012, a, nwords);
    }
}

void fp2random1012_test(digit_t* a)
{ // Generating a pseudo-random element in GF(p751^2) 
  // SECURITY NOTE: distribution is not fully uniform. TO BE USED FOR TESTING ONLY.

    fprandom1012_test(a);
    fprandom1012_test(a+NBITS_TO_NWORDS(1012));
}

#endif

#ifdef P503

extern uint64_t p503[8];

void fprandom503_test(digit_t* a)
{ // Generating a pseudo-random field element in [0, p503-1] 
  // SECURITY NOTE: distribution is not fully uniform. TO BE USED FOR TESTING ONLY.
    unsigned int i, diff = 512-503, nwords = NBITS_TO_NWORDS(503);
    unsigned char* string = NULL;

    string = (unsigned char*)a;
    for (i = 0; i < sizeof(digit_t)*nwords; i++) {
        *(string + i) = (unsigned char)rand();              // Obtain 512-bit number
    }
    a[nwords-1] &= (((digit_t)(-1) << diff) >> diff);

    while (compare_words((digit_t*)p503, a, nwords) < 1) {  // Force it to [0, modulus-1]
        sub_test(a, (digit_t*)p503, a, nwords);
    }
}


void fp2random503_test(digit_t* a)
{ // Generating a pseudo-random element in GF(p503^2) 
  // SECURITY NOTE: distribution is not fully uniform. TO BE USED FOR TESTING ONLY.

    fprandom503_test(a);
    fprandom503_test(a+NBITS_TO_NWORDS(503));
}

#endif


#ifdef P507

extern uint64_t p507[8];

void fprandom507_test(digit_t* a)
{ // Generating a pseudo-random field element in [0, p503-1] 
  // SECURITY NOTE: distribution is not fully uniform. TO BE USED FOR TESTING ONLY.
    unsigned int i, diff = 512-507, nwords = NBITS_TO_NWORDS(507);
    unsigned char* string = NULL;

    string = (unsigned char*)a;
    for (i = 0; i < sizeof(digit_t)*nwords; i++) {
        *(string + i) = (unsigned char)rand();              // Obtain 512-bit number
    }
    a[nwords-1] &= (((digit_t)(-1) << diff) >> diff);

    while (compare_words((digit_t*)p507, a, nwords) < 1) {  // Force it to [0, modulus-1]
        sub_test(a, (digit_t*)p507, a, nwords);
    }
}


void fp2random507_test(digit_t* a)
{ // Generating a pseudo-random element in GF(p503^2) 
  // SECURITY NOTE: distribution is not fully uniform. TO BE USED FOR TESTING ONLY.

    fprandom507_test(a);
    fprandom507_test(a+NBITS_TO_NWORDS(507));
}

#endif


#ifdef P508

extern uint64_t p508[8];

void fprandom508_test(digit_t* a)
{ // Generating a pseudo-random field element in [0, p503-1] 
  // SECURITY NOTE: distribution is not fully uniform. TO BE USED FOR TESTING ONLY.
    unsigned int i, diff = 512-508, nwords = NBITS_TO_NWORDS(508);
    unsigned char* string = NULL;

    string = (unsigned char*)a;
    for (i = 0; i < sizeof(digit_t)*nwords; i++) {
        *(string + i) = (unsigned char)rand();              // Obtain 512-bit number
    }
    a[nwords-1] &= (((digit_t)(-1) << diff) >> diff);

    while (compare_words((digit_t*)p508, a, nwords) < 1) {  // Force it to [0, modulus-1]
        sub_test(a, (digit_t*)p508, a, nwords);
    }
}


void fp2random508_test(digit_t* a)
{ // Generating a pseudo-random element in GF(p503^2) 
  // SECURITY NOTE: distribution is not fully uniform. TO BE USED FOR TESTING ONLY.

    fprandom508_test(a);
    fprandom508_test(a+NBITS_TO_NWORDS(508));
}

#endif


