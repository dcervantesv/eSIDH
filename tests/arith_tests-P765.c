/********************************************************************************************
* SIDH: an efficient supersingular isogeny cryptography library
*
* Abstract: testing code for field arithmetic, elliptic curve and isogeny functions
*********************************************************************************************/

#include "config.h"
#include "P765_internal.h"
#include "test_extras.h"
#include <stdio.h>


// Benchmark and test parameters  
#if defined(OPTIMIZED_GENERIC_IMPLEMENTATION) || (TARGET == TARGET_ARM) || (TARGET == TARGET_ARM64)
    #define BENCH_LOOPS           100       // Number of iterations per bench
    #define SMALL_BENCH_LOOPS     100       // Number of iterations per bench
    #define TEST_LOOPS             10       // Number of iterations per test
#else
    #define BENCH_LOOPS        500000 
    #define SMALL_BENCH_LOOPS  100000
    #define TEST_LOOPS            100  
#endif


bool fp_test()
{ // Tests for the field arithmetic
    bool OK = true;
    int n, passed;
    felm_t a, b, c, d, e, f, ma, mb, mc, md, me, mf;

    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("Testing field arithmetic over GF(p765): \n\n"); 

    // Field addition over the prime p765
    passed = 1;
    for (n=0; n<TEST_LOOPS; n++)
    {
        fprandom765_test(a); fprandom765_test(b); fprandom765_test(c); fprandom765_test(d); fprandom765_test(e); fprandom765_test(f); 

        fpadd765(a, b, d); fpadd765(d, c, e);                 // e = (a+b)+c
        fpadd765(b, c, d); fpadd765(d, a, f);                 // f = a+(b+c)
        fpcorrection765(e);
        fpcorrection765(f);
        if (compare_words(e, f, NWORDS_FIELD)!=0) { passed=0; break; }

        fpadd765(a, b, d);                                     // d = a+b 
        fpadd765(b, a, e);                                     // e = b+a
        fpcorrection765(d);
        fpcorrection765(e);
        if (compare_words(d, e, NWORDS_FIELD)!=0) { passed=0; break; }

        fpzero765(b);
        fpadd765(a, b, d);                                     // d = a+0 
        if (compare_words(a, d, NWORDS_FIELD)!=0) { passed=0; break; }
        
        fpzero765(b);
        fpcopy765(a, d);     
        fpneg765(d);                      
        fpadd765(a, d, e);                                     // e = a+(-a)
        fpcorrection765(e);
        if (compare_words(e, b, NWORDS_FIELD)!=0) { passed=0; break; }
    }
    if (passed==1) printf("  GF(p) addition tests ............................................ PASSED");
    else { printf("  GF(p) addition tests... FAILED"); printf("\n"); return false; }
    printf("\n");

    // Field subtraction over the prime p765
    passed = 1;
    for (n=0; n<TEST_LOOPS; n++)
    {
        fprandom765_test(a); fprandom765_test(b); fprandom765_test(c); fprandom765_test(d); fprandom765_test(e); fprandom765_test(f); 

        fpsub765(a, b, d); fpsub765(d, c, e);                 // e = (a-b)-c
        fpadd765(b, c, d); fpsub765(a, d, f);                 // f = a-(b+c)
        fpcorrection765(e);
        fpcorrection765(f);
        if (compare_words(e, f, NWORDS_FIELD)!=0) { passed=0; break; }

        fpsub765(a, b, d);                                     // d = a-b 
        fpsub765(b, a, e);                                         
        fpneg765(e);                                           // e = -(b-a)
        fpcorrection765(d);
        fpcorrection765(e);
        if (compare_words(d, e, NWORDS_FIELD)!=0) { passed=0; break; }

        fpzero765(b);
        fpsub765(a, b, d);                                     // d = a-0 
        if (compare_words(a, d, NWORDS_FIELD)!=0) { passed=0; break; }
        
        fpzero765(b);
        fpcopy765(a, d);                 
        fpsub765(a, d, e);                                     // e = a+(-a)
        fpcorrection765(e);
        if (compare_words(e, b, NWORDS_FIELD)!=0) { passed=0; break; }
    }
    if (passed==1) printf("  GF(p) subtraction tests ......................................... PASSED");
    else { printf("  GF(p) subtraction tests... FAILED"); printf("\n"); return false; }
    printf("\n");

    // Field multiplication over the prime p765
    passed = 1;
    for (n=0; n<TEST_LOOPS; n++)
    {    
        fprandom765_test(a); fprandom765_test(b); fprandom765_test(c);  
        fprandom765_test(ma); fprandom765_test(mb); fprandom765_test(mc); fprandom765_test(md); fprandom765_test(me); fprandom765_test(mf); 

        to_mont(a, ma);
        fpcopy765(ma, mc);
        from_mont(mc, c);
        if (compare_words(a, c, NWORDS_FIELD)!=0) { 
        passed=0; break; }
        
        to_mont(a, ma); to_mont(b, mb); to_mont(c, mc); 
        fpmul765_mont(ma, mb, md); fpmul765_mont(md, mc, me);                          // e = (a*b)*c
        fpmul765_mont(mb, mc, md); fpmul765_mont(md, ma, mf);                          // f = a*(b*c)
        from_mont(me, e);
        from_mont(mf, f);
        if (compare_words(e, f, NWORDS_FIELD)!=0) { passed=0; break; }
      
        to_mont(a, ma); to_mont(b, mb); to_mont(c, mc); 
        fpadd765(mb, mc, md); fpmul765_mont(ma, md, me);                               // e = a*(b+c)
        fpmul765_mont(ma, mb, md); fpmul765_mont(ma, mc, mf); fpadd765(md, mf, mf);    // f = a*b+a*c
        from_mont(me, e);
        from_mont(mf, f);
        if (compare_words(e, f, NWORDS_FIELD)!=0) { passed=0; break; }
       
        to_mont(a, ma); to_mont(b, mb);
        fpmul765_mont(ma, mb, md);                                                      // d = a*b 
        fpmul765_mont(mb, ma, me);                                                      // e = b*a 
        from_mont(md, d);
        from_mont(me, e);
        if (compare_words(d, e, NWORDS_FIELD)!=0) { passed=0; break; }
        
        to_mont(a, ma);
        fpzero765(b); b[0] = 1; to_mont(b, mb);
        fpmul765_mont(ma, mb, md);                                                      // d = a*1  
        from_mont(ma, a);
        from_mont(md, d);                
        if (compare_words(a, d, NWORDS_FIELD)!=0) { passed=0; break; }
        
        fpzero765(b); to_mont(b, mb);
        fpmul765_mont(ma, mb, md);                                                      // d = a*0  
        from_mont(mb, b);
        from_mont(md, d);                
        if (compare_words(b, d, NWORDS_FIELD)!=0) { passed=0; break; } 
    }
    if (passed==1) printf("  GF(p) multiplication tests ...................................... PASSED");
    else { printf("  GF(p) multiplication tests... FAILED"); printf("\n"); return false; }
    printf("\n");

    // Field squaring over the prime p765
    passed = 1;
    for (n=0; n<TEST_LOOPS; n++)
    {
        fprandom765_test(a);
        
        to_mont(a, ma);
        fpsqr765_mont(ma, mb);                                 // b = a^2
        fpmul765_mont(ma, ma, mc);                             // c = a*a 
        if (compare_words(mb, mc, NWORDS_FIELD)!=0) { passed=0; break; }

        fpzero765(a); to_mont(a, ma);
        fpsqr765_mont(ma, md);                                 // d = 0^2 
        if (compare_words(ma, md, NWORDS_FIELD)!=0) { passed=0; break; }
    }
    if (passed==1) printf("  GF(p) squaring tests............................................. PASSED");
    else { printf("  GF(p) squaring tests... FAILED"); printf("\n"); return false; }
    printf("\n");

    // Field inversion over the prime p765
    passed = 1;
    for (n=0; n<TEST_LOOPS; n++)
    {
        fprandom765_test(a); 
        to_mont(a, ma);
        fpzero765(d); d[0]=1; to_mont(d, md);
        fpcopy765(ma, mb);                            
        fpinv765_mont(ma);                                
        fpmul765_mont(ma, mb, mc);                             // c = a*a^-1 
//        fpcorrection765(mc);
        if (compare_words(mc, md, NWORDS_FIELD)!=0) { passed=0; break; }
    }
    if (passed==1) printf("  GF(p) inversion tests............................................ PASSED");
    else { printf("  GF(p) inversion tests... FAILED"); printf("\n"); return false; }
    printf("\n");
  
    return OK;
}


bool fp2_test()
{ // Tests for the quadratic extension field arithmetic
    bool OK = true;
    int n, passed;
    f2elm_t a, b, c, d, e, f, ma, mb, mc, md, me, mf;

    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("Testing quadratic extension arithmetic over GF(p765^2): \n\n"); 

    // Addition over GF(p765^2)
    passed = 1;
    for (n=0; n<TEST_LOOPS; n++)
    {
        fp2random765_test((digit_t*)a); fp2random765_test((digit_t*)b); fp2random765_test((digit_t*)c); fp2random765_test((digit_t*)d); fp2random765_test((digit_t*)e); fp2random765_test((digit_t*)f); 

        fp2add765(a, b, d); fp2add765(d, c, e);                 // e = (a+b)+c
        fp2add765(b, c, d); fp2add765(d, a, f);                 // f = a+(b+c)
        if (compare_words((digit_t*)e, (digit_t*)f, 2*NWORDS_FIELD)!=0) { passed=0; break; }

        fp2add765(a, b, d);                                     // d = a+b 
        fp2add765(b, a, e);                                     // e = b+a
        if (compare_words((digit_t*)d, (digit_t*)e, 2*NWORDS_FIELD)!=0) { passed=0; break; }

        fp2zero765(b);
        fp2add765(a, b, d);                                     // d = a+0 
        if (compare_words((digit_t*)a, (digit_t*)d, 2*NWORDS_FIELD)!=0) { passed=0; break; }
        
        fp2zero765(b);
        fp2copy765(a, d);     
        fp2neg765(d);                      
        fp2add765(a, d, e);                                     // e = a+(-a)
        if (compare_words((digit_t*)e, (digit_t*)b, 2*NWORDS_FIELD)!=0) { passed=0; break; }
    }
    if (passed==1) printf("  GF(p^2) addition tests .......................................... PASSED");
    else { printf("  GF(p^2) addition tests... FAILED"); printf("\n"); return false; }
    printf("\n");

    // Subtraction over GF(p765^2)
    passed = 1;
    for (n=0; n<TEST_LOOPS; n++)
    {
        fp2random765_test((digit_t*)a); fp2random765_test((digit_t*)b); fp2random765_test((digit_t*)c); fp2random765_test((digit_t*)d); fp2random765_test((digit_t*)e); fp2random765_test((digit_t*)f); 

        fp2sub765(a, b, d); fp2sub765(d, c, e);                 // e = (a-b)-c
        fp2add765(b, c, d); fp2sub765(a, d, f);                 // f = a-(b+c)
        if (compare_words((digit_t*)e, (digit_t*)f, 2*NWORDS_FIELD)!=0) { passed=0; break; }

        fp2sub765(a, b, d);                                     // d = a-b 
        fp2sub765(b, a, e);                                         
        fp2neg765(e);                                           // e = -(b-a)
        if (compare_words((digit_t*)d, (digit_t*)e, 2*NWORDS_FIELD)!=0) { passed=0; break; }

        fp2zero765(b);
        fp2sub765(a, b, d);                                     // d = a-0 
        if (compare_words((digit_t*)a, (digit_t*)d, 2*NWORDS_FIELD)!=0) { passed=0; break; }
        
        fp2zero765(b);
        fp2copy765(a, d);                 
        fp2sub765(a, d, e);                                     // e = a+(-a)
        if (compare_words((digit_t*)e, (digit_t*)b, 2*NWORDS_FIELD)!=0) { passed=0; break; }
    }
    if (passed==1) printf("  GF(p^2) subtraction tests ....................................... PASSED");
    else { printf("  GF(p^2) subtraction tests... FAILED"); printf("\n"); return false; }
    printf("\n");

    // Multiplication over GF(p765^2)
    passed = 1;
    for (n=0; n<TEST_LOOPS; n++)
    {    
        fp2random765_test((digit_t*)a); fp2random765_test((digit_t*)b); fp2random765_test((digit_t*)c);  
        fp2random765_test((digit_t*)ma); fp2random765_test((digit_t*)mb); fp2random765_test((digit_t*)mc); fp2random765_test((digit_t*)md); fp2random765_test((digit_t*)me); fp2random765_test((digit_t*)mf); 

        to_fp2mont(a, ma);
        fp2copy765(ma, mc);
        from_fp2mont(mc, c);
        if (compare_words((digit_t*)a, (digit_t*)c, 2*NWORDS_FIELD)!=0) { passed=0; break; }
        
        to_fp2mont(a, ma); to_fp2mont(b, mb); to_fp2mont(c, mc); 
        fp2mul765_mont(ma, mb, md); fp2mul765_mont(md, mc, me);                          // e = (a*b)*c
        fp2mul765_mont(mb, mc, md); fp2mul765_mont(md, ma, mf);                          // f = a*(b*c)
        from_fp2mont(me, e);
        from_fp2mont(mf, f);
        if (compare_words((digit_t*)e, (digit_t*)f, 2*NWORDS_FIELD)!=0) { passed=0; break; }
      
        to_fp2mont(a, ma); to_fp2mont(b, mb); to_fp2mont(c, mc); 
        fp2add765(mb, mc, md); fp2mul765_mont(ma, md, me);                               // e = a*(b+c)
        fp2mul765_mont(ma, mb, md); fp2mul765_mont(ma, mc, mf); fp2add765(md, mf, mf);   // f = a*b+a*c
        from_fp2mont(me, e);
        from_fp2mont(mf, f);
        if (compare_words((digit_t*)e, (digit_t*)f, 2*NWORDS_FIELD)!=0) { passed=0; break; }
       
        to_fp2mont(a, ma); to_fp2mont(b, mb);
        fp2mul765_mont(ma, mb, md);                                                      // d = a*b 
        fp2mul765_mont(mb, ma, me);                                                      // e = b*a 
        from_fp2mont(md, d);
        from_fp2mont(me, e);
        if (compare_words((digit_t*)d, (digit_t*)e, 2*NWORDS_FIELD)!=0) { passed=0; break; }
        
        to_fp2mont(a, ma);
        fp2zero765(b); b[0][0] = 1; to_fp2mont(b, mb);
        fp2mul765_mont(ma, mb, md);                                                      // d = a*1  
        from_fp2mont(md, d);               
        if (compare_words((digit_t*)a, (digit_t*)d, 2*NWORDS_FIELD)!=0) { passed=0; break; }
        
        fp2zero765(b); to_fp2mont(b, mb);
        fp2mul765_mont(ma, mb, md);                                                      // d = a*0 
        from_fp2mont(md, d);               
        if (compare_words((digit_t*)b, (digit_t*)d, 2*NWORDS_FIELD)!=0) { passed=0; break; } 
    }
    if (passed==1) printf("  GF(p^2) multiplication tests .................................... PASSED");
    else { printf("  GF(p^2) multiplication tests... FAILED"); printf("\n"); return false; }
    printf("\n");

    // Squaring over GF(p765^2)
    passed = 1;
    for (n=0; n<TEST_LOOPS; n++)
    {
        fp2random765_test((digit_t*)a);
        
        to_fp2mont(a, ma);
        fp2sqr765_mont(ma, mb);                                 // b = a^2
        fp2mul765_mont(ma, ma, mc);                             // c = a*a 
        from_fp2mont(mb, b);               
        from_fp2mont(mc, c);               
        if (compare_words((digit_t*)b, (digit_t*)c, 2*NWORDS_FIELD)!=0) { passed=0; break; }

        fp2zero765(a); to_fp2mont(a, ma);
        fp2sqr765_mont(ma, md);                                 // d = 0^2 
        from_fp2mont(md, d);               
        if (compare_words((digit_t*)a, (digit_t*)d, 2*NWORDS_FIELD)!=0) { passed=0; break; }
    }
    if (passed==1) printf("  GF(p^2) squaring tests........................................... PASSED");
    else { printf("  GF(p^2) squaring tests... FAILED"); printf("\n"); return false; }
    printf("\n");
    
    // Inversion over GF(p765^2)
    passed = 1;
    for (n=0; n<TEST_LOOPS; n++)
    {
        fp2random765_test((digit_t*)a);    
        
        to_fp2mont(a, ma);
        fp2zero765(d); d[0][0]=1; to_fp2mont(d, md);
        fp2copy765(ma, mb);                            
        fp2inv765_mont(ma);                                
        fp2mul765_mont(ma, mb, mc);                             // c = a*a^-1              
        from_fp2mont(mc, c);  
        if (compare_words((digit_t*)c, (digit_t*)d, 2*NWORDS_FIELD)!=0) { passed=0; break; }
    }
    if (passed==1) printf("  GF(p^2) inversion tests.......................................... PASSED");
    else { printf("  GF(p^2) inversion tests... FAILED"); printf("\n"); return false; }
    printf("\n");
    
    return OK;
}


bool fp_run()
{
    bool OK = true;
    int n;
    unsigned long long cycles, cycles1, cycles2;
    felm_t a, b, c;
    dfelm_t aa;
        
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("Benchmarking field arithmetic over GF(p765): \n\n"); 
        
    fprandom765_test(a); fprandom765_test(b); fprandom765_test(c);

    // GF(p) addition using p765
    cycles = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles(); 
        fpadd765(a, b, c);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  GF(p) addition runs in .......................................... %7lld ", cycles/BENCH_LOOPS); print_unit;
    printf("\n");

    // GF(p) subtraction using p765
    cycles = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles(); 
        fpsub765(a, b, c);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  GF(p) subtraction runs in ....................................... %7lld ", cycles/BENCH_LOOPS); print_unit;
    printf("\n");

    // GF(p) multiplication using p765
    cycles = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles(); 
        fpmul765_mont(a, b, c);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  GF(p) multiplication runs in .................................... %7lld ", cycles/BENCH_LOOPS); print_unit;
    printf("\n");

    // GF(p) reduction using p765
    cycles = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        mp_mul(a, b, aa, NWORDS_FIELD);

        cycles1 = cpucycles(); 
        rdc_mont(aa, c);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  GF(p) reduction runs in ......................................... %7lld ", cycles/BENCH_LOOPS); print_unit;
    printf("\n");

    // GF(p) inversion
    cycles = 0;
    for (n=0; n<SMALL_BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles(); 
        fpinv765_mont(a);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  GF(p) inversion (exponentiation) runs in ........................ %7lld ", cycles/SMALL_BENCH_LOOPS); print_unit;
    printf("\n");
    
    return OK;
}


bool fp2_run()
{
    bool OK = true;
    int n;
    unsigned long long cycles, cycles1, cycles2;
    f2elm_t a, b, c;
        
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("Benchmarking quadratic extension arithmetic over GF(p765^2): \n\n"); 
    
    fp2random765_test((digit_t*)a); fp2random765_test((digit_t*)b); fp2random765_test((digit_t*)c);

    // GF(p^2) addition
    cycles = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles(); 
        fp2add765(a, b, c);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  GF(p^2) addition runs in ........................................ %7lld ", cycles/BENCH_LOOPS); print_unit;
    printf("\n");

    // GF(p^2) subtraction
    cycles = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles(); 
        fp2sub765(a, b, c);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  GF(p^2) subtraction runs in ..................................... %7lld ", cycles/BENCH_LOOPS); print_unit;
    printf("\n");

    // GF(p^2) multiplication
    cycles = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles(); 
        fp2mul765_mont(a, b, c);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  GF(p^2) multiplication runs in .................................. %7lld ", cycles/BENCH_LOOPS); print_unit;
    printf("\n");

    // GF(p^2) squaring
    cycles = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles(); 
        fp2sqr765_mont(a, b);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  GF(p^2) squaring runs in ........................................ %7lld ", cycles/BENCH_LOOPS); print_unit;
    printf("\n");

    // GF(p^2) inversion
    cycles = 0;
    for (n=0; n<SMALL_BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles(); 
        fp2inv765_mont(a);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  GF(p^2) inversion (exponentiation) runs in ...................... %7lld ", cycles/SMALL_BENCH_LOOPS); print_unit;
    printf("\n"); 
    
    return OK;
}


bool ecisog_run()
{
    bool OK = true;
    int n;
    unsigned long long cycles, cycles1, cycles2;
    f2elm_t A24, C24, A24m, A24p, A4, A, C, coeff[5];
    point_proj_t P, Q, KPs[KERNEL_POINTS2];
        
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("Benchmarking elliptic curve and isogeny functions: \n\n"); 

    // Point doubling
    cycles = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        fp2random765_test((digit_t*)A24); fp2random765_test((digit_t*)C24);

        cycles1 = cpucycles(); 
        xDBL(P, Q, A24, C24);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  Point doubling runs in .......................................... %7lld ", cycles/BENCH_LOOPS); print_unit;
    printf("\n");

    // 4-isogeny of a projective point
    cycles = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        fp2random765_test((digit_t*)A); fp2random765_test((digit_t*)coeff[0]); fp2random765_test((digit_t*)coeff[1]); fp2random765_test((digit_t*)coeff[2]);

        cycles1 = cpucycles(); 
        get_4_isog(P, A, C, coeff);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  4-isogeny of projective point runs in ........................... %7lld ", cycles/BENCH_LOOPS); print_unit;
    printf("\n");

    // 4-isogeny evaluation at projective point
    cycles = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        fp2random765_test((digit_t*)A); fp2random765_test((digit_t*)coeff[0]); fp2random765_test((digit_t*)coeff[1]); fp2random765_test((digit_t*)coeff[2]);

        cycles1 = cpucycles(); 
        eval_4_isog(P, coeff);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  4-isogeny evaluation at projective point runs in ................ %7lld ", cycles/BENCH_LOOPS); print_unit;
    printf("\n");

    // Point tripling
    cycles = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        fp2random765_test((digit_t*)A4); fp2random765_test((digit_t*)C);

        cycles1 = cpucycles(); 
        xTPL(P, Q, A4, C);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  Point tripling runs in .......................................... %7lld ", cycles/BENCH_LOOPS); print_unit;
    printf("\n");

    // 3-isogeny of a projective point
    cycles = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        fp2random765_test((digit_t*)A); fp2random765_test((digit_t*)C);

        cycles1 = cpucycles(); 
        get_3_isog(P, A, C, coeff);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  3-isogeny of projective point runs in ........................... %7lld ", cycles/BENCH_LOOPS); print_unit;
    printf("\n");

    // 3-isogeny evaluation at projective point
    cycles = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles(); 
        eval_3_isog(Q, coeff);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  3-isogeny evaluation at projective point runs in ................ %7lld ", cycles/BENCH_LOOPS); print_unit;
    printf("\n");
    
    // Point quintupling
    cycles = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        fp2random765_test((digit_t*)A24m); fp2random765_test((digit_t*)A24p); fp2random765_test((digit_t*)C);

        cycles1 = cpucycles(); 
        xQPL(P, Q, A24p, A24m, C);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  Point quintupling runs in ....................................... %7lld ", cycles/BENCH_LOOPS); print_unit;
    printf("\n");

    // 5-isogeny of a projective point
    cycles = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        fp2random765_test((digit_t*)A24m); fp2random765_test((digit_t*)A24p); fp2random765_test((digit_t*)C);

        cycles1 = cpucycles(); 
        get_5_isog(P,  KPs,  A24p, A24m, C);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  5-isogeny of projective point runs in ........................... %7lld ", cycles/BENCH_LOOPS); print_unit;
    printf("\n");

    // 5-isogeny evaluation at projective point
    cycles = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles(); 
        eval_5_isog(KPs, Q);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  5-isogeny evaluation at projective point runs in ................ %7lld ", cycles/BENCH_LOOPS); print_unit;
    printf("\n");    
    
    return OK;
}

int main()
{
    bool OK = true;
    
    OK = OK && fp_test();          // Test field operations using p765
    OK = OK && fp_run();           // Benchmark field operations using p765

    OK = OK && fp2_test();         // Test arithmetic functions over GF(p765^2)
    OK = OK && fp2_run();          // Benchmark arithmetic functions over GF(p765^2)
   
    OK = OK && ecisog_run();       // Benchmark elliptic curve and isogeny functions

    return OK;
}
