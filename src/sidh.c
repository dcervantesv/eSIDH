/********************************************************************************************
* SIDH: an efficient supersingular isogeny cryptography library
*
* Abstract: ephemeral supersingular isogeny Diffie-Hellman key exchange (SIDH)
*********************************************************************************************/ 

#include "random/random.h"

#if defined PARALLEL

#include <omp.h>
#include<pthread.h>
#endif

static void clear_words(void* mem, digit_t nwords)
{ // Clear digits from memory. "nwords" indicates the number of digits to be zeroed.
  // This function uses the volatile type qualifier to inform the compiler not to optimize out the memory clearing.
    unsigned int i;
    volatile digit_t *v = mem; 

    for (i = 0; i < nwords; i++) {
        v[i] = 0;
    }
}


static void init_basis(digit_t *gen, f2elm_t XP, f2elm_t XQ, f2elm_t XR)
{ // Initialization of basis points

    fpcopy(gen,                  XP[0]);
    fpcopy(gen +   NWORDS_FIELD, XP[1]);
    fpcopy(gen + 2*NWORDS_FIELD, XQ[0]);
    fpzero(XQ[1]);
    fpcopy(gen + 3*NWORDS_FIELD, XR[0]);
    fpcopy(gen + 4*NWORDS_FIELD, XR[1]);
}

static void fp2_encode(const f2elm_t x, unsigned char *enc)
{ // Conversion of GF(p^2) element from Montgomery to standard representation, and encoding by removing leading 0 bytes
    unsigned int i;
    f2elm_t t;

    from_fp2mont(x, t);

    for (i = 0; i < FP2_ENCODED_BYTES / 2; i++) {
        enc[i] = ((unsigned char*)t)[i];
        enc[i + FP2_ENCODED_BYTES / 2] = ((unsigned char*)t)[i + MAXBITS_FIELD / 8];
    }
}


static void fp2_decode(const unsigned char *enc, f2elm_t x)
{ // Parse byte sequence back into GF(p^2) element, and conversion to Montgomery representation
    unsigned int i;

    for (i = 0; i < 2*(MAXBITS_FIELD / 8); i++) ((unsigned char *)x)[i] = 0;
    for (i = 0; i < FP2_ENCODED_BYTES / 2; i++) {
        ((unsigned char*)x)[i] = enc[i];
        ((unsigned char*)x)[i + MAXBITS_FIELD / 8] = enc[i + FP2_ENCODED_BYTES / 2];
    }
    to_fp2mont(x, x);
}

void random_mod_order_A(unsigned char* random_digits)
{  // Generation of Alice's secret key  
   // Outputs random value in [0, 2^eA - 1]
    unsigned long long nbytes = NBITS_TO_NBYTES(OALICE_BITS);

    clear_words((void*)random_digits, MAXWORDS_ORDER);
    randombytes(random_digits, nbytes);
    random_digits[nbytes-1] &= MASK_ALICE;    // Masking last byte 
    

}

#if ((BOB_PRIMES == 1))

void random_mod_order_B(unsigned char* random_digits)
{  // Generation of Bob's secret key  
   // Outputs random value in [0, 2^Floor(Log(2, oB)) - 1]
    unsigned long long nbytes = NBITS_TO_NBYTES(OBOB_BITS-1);

    clear_words((void*)random_digits, MAXWORDS_ORDER);
    randombytes(random_digits, nbytes);
    random_digits[nbytes-1] &= MASK_BOB;     // Masking last byte 
}

#else

void random_mod_order_B(unsigned char* random_digits)
{  // Generation of Bob's secret key  
   // Outputs random value in [0, 2^Floor(Log(2, oB)) - 1]
    unsigned long long nbytes = NBITS_TO_NBYTES(MAXWORDS_ORDER*64);

    clear_words((void*)random_digits, MAXWORDS_ORDER);
    randombytes(random_digits, nbytes);
    random_digits[(nbytes/2)-1] &= MASK_BOB;     // Masking last byte
    random_digits[nbytes-1] &= MASK_BOB1;             // Masking last byte
}

#endif





int EphemeralKeyGeneration_A(const unsigned char* PrivateKeyA, unsigned char* PublicKeyA)
{ // Alice's ephemeral public key generation
  // Input:  a private key PrivateKeyA in the range [0, 2^eA - 1]. 
  // Output: the public key PublicKeyA consisting of 3 elements in GF(p^2) which are encoded by removing leading 0 bytes.
    point_proj_t R, phiP = {0}, phiQ = {0}, phiR = {0}, pts[MAX_INT_POINTS_ALICE];
    f2elm_t XPA, XQA, XRA, coeff[3], A24plus = {0}, C24 = {0}, A = {0};
    unsigned int i, row, m, index = 0, pts_index[MAX_INT_POINTS_ALICE], npts = 0, ii = 0, Bound;
    #if defined PARALLEL
        int tid, cn;
        cn = omp_get_max_threads();
    #endif
        
    #if defined NEWSTRAT
        int b = MAX_Alice - strat_Alice[0], availablecores;
        point_proj_t R_tmp, Pts_tmp[4], Rt;
        f2elm_t coeffs[b+1][3];
        ladder_struct *ladder_stuff = malloc(4*sizeof(f2elm_t) + sizeof(digit_t) + sizeof(point_proj_t) + sizeof(unsigned int)); 
        void *ptr =NULL;
        int bt;
    #endif    
        
    // Initialize basis points
    init_basis((digit_t*)A_gen, XPA, XQA, XRA);
    init_basis((digit_t*)B_gen, phiP->X, phiQ->X, phiR->X);
    fpcopy((digit_t*)&Montgomery_one, (phiP->Z)[0]);
    fpcopy((digit_t*)&Montgomery_one, (phiQ->Z)[0]);
    fpcopy((digit_t*)&Montgomery_one, (phiR->Z)[0]);

    // Initialize constants
    fpcopy((digit_t*)&Montgomery_one, A24plus[0]);
    fp2add(A24plus, A24plus, C24);
#if defined NEWSTRAT
    availablecores = cn - 1;
    
    fp2copy(XPA, ladder_stuff->xP);
    fp2copy(XQA, ladder_stuff->xQ);
    fp2copy(XRA, ladder_stuff->xPQ);
    fp2copy(A, ladder_stuff->A);
    ladder_stuff->AliceOrBob = ALICE;
    *(ladder_stuff->R)=*R;
    ladder_stuff->m=(digit_t *)PrivateKeyA;
    pthread_t ladder_t;
    int lpt, j;
    //thread computing R
    lpt = pthread_create(&ladder_t, NULL, LADDER3PT_ptNS, (void *)ladder_stuff);
    //Processing St_h
            LADDER3PT_NS_prcmp(&PA_prcmp[2*NWORDS_FIELD * strat_Alice[0]], &RA_prcmp[2*NWORDS_FIELD * strat_Alice[0]], (digit_t*)PrivateKeyA, ALICE, R_tmp, A, 2*b);
            index = 0;  
            ii=1;
            fp2copy(R_tmp->X, pts[0]->X);
            fp2copy(R_tmp->Z, pts[0]->Z);
            pts_index[0] = 0;
            npts = 1;
            for (row = 0; row < b - 1; row++)
            {
                while (index < b - 1 - row) {
                    m = strat_Alice[ii];
                    index += m;
                    pts_index[npts] = index;
                    xDBLe(R_tmp, R_tmp, A24plus, C24, (int)(2*m));
                    fp2copy(R_tmp->X, pts[npts]->X);
                    fp2copy(R_tmp->Z, pts[npts]->Z);
                    ii +=1;
                    npts +=1;
                }
                get_4_isog(R_tmp, A24plus, C24, coeffs[row]);     
                #pragma omp parallel num_threads(availablecores) private(tid,i) shared(pts, npts, coeffs, availablecores) 
                {   
                    tid = omp_get_thread_num();
                    for (i = tid; i < npts - 1; i = i + availablecores) {
                        eval_4_isog(pts[i], coeffs[row]);
                    }
                }
                fp2copy(pts[npts-2]->X, R_tmp->X); 
                fp2copy(pts[npts-2]->Z, R_tmp->Z);
                index = pts_index[npts-2];
                npts -= 1;
            }
    bt = pthread_join(ladder_t, &ptr);
    //This could be improved
    fp2copy(ladder_stuff->R->X, Pts_tmp[0]->X);
    fp2copy(ladder_stuff->R->Z, Pts_tmp[0]->Z);
    fp2copy(phiP->X, Pts_tmp[1]->X);
    fp2copy(phiP->Z, Pts_tmp[1]->Z);
    fp2copy(phiQ->X, Pts_tmp[2]->X);
    fp2copy(phiQ->Z, Pts_tmp[2]->Z);
    fp2copy(phiR->X, Pts_tmp[3]->X);
    fp2copy(phiR->Z, Pts_tmp[3]->Z);

    
    //sequential step
    get_4_isog(R_tmp, A24plus, C24, coeffs[b-1]);
    for(i = 0; i<b; i++){
        #pragma omp parallel for
        for(j=0; j < 4; j++){
            eval_4_isog(Pts_tmp[j], coeffs[i]);
    }
    }
    fp2copy(Pts_tmp[0]->X, R->X );
    fp2copy(Pts_tmp[0]->Z, R->Z );
    fp2copy(Pts_tmp[1]->X, phiP->X);
    fp2copy(Pts_tmp[1]->Z, phiP->Z);
    fp2copy(Pts_tmp[2]->X, phiQ->X);
    fp2copy(Pts_tmp[2]->Z, phiQ->Z);
    fp2copy(Pts_tmp[3]->X, phiR->X);
    fp2copy(Pts_tmp[3]->Z, phiR->Z);
    ii=b;
    Bound = strat_Alice[0];
#else
    // Retrieve kernel point
    LADDER3PT(XPA, XQA, XRA, (digit_t*)PrivateKeyA, ALICE, R, A);
    Bound = MAX_Alice;
#endif
    // Traverse tree
    index = 0;     
    fp2copy(R->X, pts[0]->X);
    fp2copy(R->Z, pts[0]->Z);
    pts_index[0] = 0;
    npts = 1;
    for (row = 0; row < Bound - 1; row++)
    {
        while (index < Bound - 1 - row) {
            m = strat_Alice[ii];
            ii +=1;
            index += m;
            pts_index[npts] = index;
            xDBLe(R, R, A24plus, C24, (int)(2*m));
            fp2copy(R->X, pts[npts]->X);
            fp2copy(R->Z, pts[npts]->Z);
            npts +=1;
        }
        get_4_isog(R, A24plus, C24, coeff);        
#if defined PARALLEL
        #pragma omp parallel num_threads(cn) private(tid,i) shared(pts, npts, coeff, cn) 
        {   
            tid = omp_get_thread_num();
            for (i = tid; i < npts-1; i = i + cn) {
                eval_4_isog(pts[i], coeff);
            }
        }
#else
        for (i = 0; i < npts-1; i++)
        {
            eval_4_isog(pts[i], coeff);
        }
#endif
        eval_4_isog(phiP, coeff);
        eval_4_isog(phiQ, coeff);
        eval_4_isog(phiR, coeff);

        fp2copy(pts[npts-2]->X, R->X); 
        fp2copy(pts[npts-2]->Z, R->Z);
        index = pts_index[npts-2];
        npts -= 1;
    }

    get_4_isog(R, A24plus, C24, coeff); 
    eval_4_isog(phiP, coeff);
    eval_4_isog(phiQ, coeff);
    eval_4_isog(phiR, coeff);

    inv_3_way(phiP->Z, phiQ->Z, phiR->Z);
    fp2mul_mont(phiP->X, phiP->Z, phiP->X);
    fp2mul_mont(phiQ->X, phiQ->Z, phiQ->X);
    fp2mul_mont(phiR->X, phiR->Z, phiR->X);
                
    // Format public key                   
    fp2_encode(phiP->X, PublicKeyA);
    fp2_encode(phiQ->X, PublicKeyA + FP2_ENCODED_BYTES);
    fp2_encode(phiR->X, PublicKeyA + 2*FP2_ENCODED_BYTES);

    return 0;
}

int EphemeralSecretAgreement_A(const unsigned char* PrivateKeyA, const unsigned char* PublicKeyB, unsigned char* SharedSecretA)
{ 
 // Alice's ephemeral shared secret computation
  // It produces a shared secret key SharedSecretA using her secret key PrivateKeyA and Bob's public key PublicKeyB
  // Inputs: Alice's PrivateKeyA is an integer in the range [0, oA-1]. 
  //         Bob's PublicKeyB consists of 3 elements in GF(p^2) encoded by removing leading 0 bytes.
  // Output: a shared secret SharedSecretA that consists of one element in GF(p^2) encoded by removing leading 0 bytes.  
    point_proj_t R, Rt, pts[MAX_INT_POINTS_ALICE];
    f2elm_t coeff[3], PKB[3], jinv;
    f2elm_t A24plus = {0}, C24 = {0}, A = {0};
    unsigned int i, row, m, index = 0, pts_index[MAX_INT_POINTS_ALICE], npts = 0, ii = 0, Bound;
    #if defined PARALLEL
        int tid, cn;
        cn = omp_get_max_threads();
    #endif
    #if defined NEWSTRAT
        int b = MAX_Alice - strat_Alice[0], availablecores;
        point_proj_t R_tmp;
        f2elm_t coeffs[b][3];
        ladder_struct *ladder_stuff = malloc(4*sizeof(f2elm_t) + sizeof(digit_t) + sizeof(point_proj_t) + sizeof(unsigned int)); 
        void *ptr =NULL;
        int bt;
    #endif
    
    // Initialize images of Bob's basis
    fp2_decode(PublicKeyB, PKB[0]);
    fp2_decode(PublicKeyB + FP2_ENCODED_BYTES, PKB[1]);
    fp2_decode(PublicKeyB + 2*FP2_ENCODED_BYTES, PKB[2]);

    // Initialize constants
    get_A(PKB[0], PKB[1], PKB[2], A); // TODO: Can return projective A?
    fpadd((digit_t*)&Montgomery_one, (digit_t*)&Montgomery_one, C24[0]);
    fp2add(A, C24, A24plus);
    fpadd(C24[0], C24[0], C24[0]);

    // Retrieve kernel point
#if defined NEWSTRAT
    availablecores = cn - 1;
    fp2copy(PKB[0], ladder_stuff->xP);
    fp2copy(PKB[1], ladder_stuff->xQ);
    fp2copy(PKB[2], ladder_stuff->xPQ);
    fp2copy(A, ladder_stuff->A);
    ladder_stuff->AliceOrBob = ALICE;
    *(ladder_stuff->R)=*R;
    ladder_stuff->m=(digit_t *)PrivateKeyA;
    pthread_t ladder_t;
    int lpt;
    lpt = pthread_create(&ladder_t, NULL, LADDER3PT_pt, (void *)ladder_stuff);
            LADDER3PT_NS(PKB[0], PKB[1], PKB[2], (digit_t*)PrivateKeyA, ALICE, R_tmp, A, 2*b);
            xDBLe(R_tmp, R_tmp, A24plus, C24, (int)(2*strat_Alice[0]));
            // Traverse tree
            index = 0;  
            ii=1;
            fp2copy(R_tmp->X, pts[0]->X);
            fp2copy(R_tmp->Z, pts[0]->Z);
            pts_index[0] = 0;
            npts = 1;
            for (row = 0; row < b - 1; row++)
            {
                while (index < b - 1 - row) {
                    m = strat_Alice[ii];
                    index += m;
                    pts_index[npts] = index;
                    xDBLe(R_tmp, R_tmp, A24plus, C24, (int)(2*m));
                    fp2copy(R_tmp->X, pts[npts]->X);
                    fp2copy(R_tmp->Z, pts[npts]->Z);
                    ii +=1;
                    npts +=1;

                }
                get_4_isog(R_tmp, A24plus, C24, coeffs[row]);     
                #pragma omp parallel num_threads(availablecores) private(tid,i) shared(pts, npts, coeffs, availablecores) 
                {   
                    tid = omp_get_thread_num();
                    for (i = tid; i < npts - 1; i = i + availablecores) {
                        eval_4_isog(pts[i], coeffs[row]);
                    }
                }
                fp2copy(pts[npts-2]->X, R_tmp->X); 
                fp2copy(pts[npts-2]->Z, R_tmp->Z);
                index = pts_index[npts-2];
                npts -= 1;
            }
    bt = pthread_join(ladder_t, &ptr);

    fp2copy(ladder_stuff->R->X, R->X);
    fp2copy(ladder_stuff->R->Z, R->Z);
    for(i= 0; i<b-1; i++)
        eval_4_isog(R,coeffs[i]);
    get_4_isog(R_tmp, A24plus, C24, coeff);
    eval_4_isog(R,coeff);
    ii=b;
    Bound = strat_Alice[0];
#else
    LADDER3PT(PKB[0], PKB[1], PKB[2], (digit_t*)PrivateKeyA, ALICE, R, A);
    Bound = MAX_Alice;
#endif
    // Traverse tree
    index = 0;   
    fp2copy(R->X, pts[0]->X);
    fp2copy(R->Z, pts[0]->Z);
    pts_index[0] = 0;
    npts = 1;
    for (row = 0; row < Bound - 1; row++)
    {
        while (index < Bound - 1 - row) {
            m = strat_Alice[ii];
            ii +=1;
            index += m;
            pts_index[npts] = index;
            xDBLe(R, R, A24plus, C24, (int)(2*m));
            fp2copy(R->X, pts[npts]->X);
            fp2copy(R->Z, pts[npts]->Z);
            npts +=1;
        }
        get_4_isog(R, A24plus, C24, coeff);        

#if defined PARALLEL
        #pragma omp parallel num_threads(cn) private(tid,i) shared(pts, npts, coeff, cn) 
        {   
            tid = omp_get_thread_num();
            for (i = tid; i < npts-1; i = i + cn) {
                eval_4_isog(pts[i], coeff);
            }
        }
#else
        for (i = 0; i < npts-1; i++)
        {
            eval_4_isog(pts[i], coeff);
        }
#endif

        fp2copy(pts[npts-2]->X, R->X); 
        fp2copy(pts[npts-2]->Z, R->Z);
        index = pts_index[npts-2];
        npts -= 1;
    }
    get_4_isog(R, A24plus, C24, coeff); 
    fp2div2(C24, C24);                                                
    fp2sub(A24plus, C24, A24plus);                              
    fp2div2(C24, C24);                               
    j_inv(A24plus, C24, jinv);
    fp2_encode(jinv, SharedSecretA);    // Format shared secret

    return 0;
}

 #if ((BOB_PRIMES == 1))
int EphemeralKeyGeneration_B(const unsigned char* PrivateKeyB, unsigned char* PublicKeyB)
{ // Bob's ephemeral public key generation
  // Input:  a private key PrivateKeyB in the range [0, 2^Floor(Log(2,oB)) - 1]. 
  // Output: the public key PublicKeyB consisting of 3 elements in GF(p^2) which are encoded by removing leading 0 bytes.
    point_proj_t R, phiP = {0}, phiQ = {0}, phiR = {0}, pts[MAX_INT_POINTS_BOB];
    f2elm_t XPB, XQB, XRB, coeff[3], A24plus = {0}, A24minus = {0}, A = {0};
    unsigned int i, row, m, index = 0, pts_index[MAX_INT_POINTS_BOB], npts = 0, ii = 0;
    
    #if defined PARALLEL
        int tid, cn;
        cn = omp_get_max_threads();
    #endif
    
    // Initialize basis points
    init_basis((digit_t*)B_gen, XPB, XQB, XRB);
    init_basis((digit_t*)A_gen, phiP->X, phiQ->X, phiR->X);
    
    
    fpcopy((digit_t*)&Montgomery_one, (phiP->Z)[0]);
    fpcopy((digit_t*)&Montgomery_one, (phiQ->Z)[0]);
    fpcopy((digit_t*)&Montgomery_one, (phiR->Z)[0]);

    // Initialize constants
    fpcopy((digit_t*)&Montgomery_one, A24plus[0]);
    fp2add(A24plus, A24plus, A24plus);
    fp2copy(A24plus, A24minus);
    fp2neg(A24minus);

    // Retrieve kernel point
    LADDER3PT(XPB, XQB, XRB, (digit_t*)PrivateKeyB, BOB, R, A);
    

    index = 0; 
    fp2copy(R->X, pts[0]->X);
    fp2copy(R->Z, pts[0]->Z);
    pts_index[0] = 0;
    npts = 1;
    for (row = 0; row < MAX_Bob - 1; row++) {
         while (index < MAX_Bob - 1 - row) {
            m = strat_Bob[ii];
            ii +=1;
            index += m;
            pts_index[npts] = index;
            xMUL1e(R, R, A24minus, A24plus, (int)m);
            fp2copy(R->X, pts[npts]->X);
            fp2copy(R->Z, pts[npts]->Z);
            npts +=1;
        }
        
        get_d_isog(R, A24minus, A24plus, coeff);
#if defined PARALLEL

        #pragma omp parallel for num_threads(cn)
        for (i = 0; i < npts-1; i++) {
            eval_d_isog(pts[i], coeff);
        }

#else
        for (i = 0; i < npts-1; i++) {
            eval_d_isog(pts[i], coeff);
        } 
#endif  
        eval_d_isog(phiP, coeff);
        eval_d_isog(phiQ, coeff);
        eval_d_isog(phiR, coeff);
        fp2copy(pts[npts-2]->X, R->X); 
        fp2copy(pts[npts-2]->Z, R->Z);
        index = pts_index[npts-2];
        npts -= 1;
    }
    get_d_isog(R, A24minus, A24plus, coeff);
    eval_d_isog(phiP, coeff);
    eval_d_isog(phiQ, coeff);
    eval_d_isog(phiR, coeff);

    inv_3_way(phiP->Z, phiQ->Z, phiR->Z);
    fp2mul_mont(phiP->X, phiP->Z, phiP->X);
    fp2mul_mont(phiQ->X, phiQ->Z, phiQ->X);
    fp2mul_mont(phiR->X, phiR->Z, phiR->X);

    // Format public key
    fp2_encode(phiP->X, PublicKeyB);
    fp2_encode(phiQ->X, PublicKeyB + FP2_ENCODED_BYTES);
    fp2_encode(phiR->X, PublicKeyB + 2*FP2_ENCODED_BYTES);

    return 0;
}

int EphemeralSecretAgreement_B(const unsigned char* PrivateKeyB, const unsigned char* PublicKeyA, unsigned char* SharedSecretB)
{ // Bob's ephemeral shared secret computation
  // It produces a shared secret key SharedSecretB using his secret key PrivateKeyB and Alice's public key PublicKeyA
  // Inputs: Bob's PrivateKeyB is an integer in the range [0, 2^Floor(Log(2,oB)) - 1]. 
  //         Alice's PublicKeyA consists of 3 elements in GF(p^2) encoded by removing leading 0 bytes.
  // Output: a shared secret SharedSecretB that consists of one element in GF(p^2) encoded by removing leading 0 bytes.  
    point_proj_t R, pts[MAX_INT_POINTS_BOB];
    f2elm_t coeff[3], PKB[3], jinv;
    f2elm_t A24plus = {0}, A24minus = {0}, A = {0};
    unsigned int i, row, m, index = 0, pts_index[MAX_INT_POINTS_BOB], npts = 0, ii = 0;
    #if defined PARALLEL
        int tid, cn;
        cn = omp_get_max_threads();
    #endif
    
    // Initialize images of Alice's basis
    fp2_decode(PublicKeyA, PKB[0]);
    fp2_decode(PublicKeyA + FP2_ENCODED_BYTES, PKB[1]);
    fp2_decode(PublicKeyA + 2*FP2_ENCODED_BYTES, PKB[2]);

    // Initialize constants
    get_A(PKB[0], PKB[1], PKB[2], A); // TODO: Can return projective A?
    fpadd((digit_t*)&Montgomery_one, (digit_t*)&Montgomery_one, A24minus[0]);
    fp2add(A, A24minus, A24plus);
    fp2sub(A, A24minus, A24minus);

    // Retrieve kernel point
    LADDER3PT(PKB[0], PKB[1], PKB[2], (digit_t*)PrivateKeyB, BOB, R, A);
    
//     Traverse tree
    index = 0; 
    fp2copy(R->X, pts[0]->X);
    fp2copy(R->Z, pts[0]->Z);
    pts_index[0] = 0;
    npts = 1;
    for (row = 0; row < MAX_Bob - 1; row++) {
         while (index < MAX_Bob - 1 - row) {
            m = strat_Bob[ii];
            ii +=1;
            index += m;
            pts_index[npts] = index;
            xMUL1e(R, R, A24minus, A24plus, (int)m);
            fp2copy(R->X, pts[npts]->X);
            fp2copy(R->Z, pts[npts]->Z);
            npts +=1;
        }
        
        get_d_isog(R, A24minus, A24plus, coeff);
#if defined PARALLEL

        #pragma omp parallel for num_threads(cn)
        for (i = 0; i < npts-1; i++) {
            eval_d_isog(pts[i], coeff);
        }

#else
        for (i = 0; i < npts-1; i++) {
            eval_d_isog(pts[i], coeff);
        } 
#endif
        fp2copy(pts[npts-2]->X, R->X); 
        fp2copy(pts[npts-2]->Z, R->Z);
        index = pts_index[npts-2];
        npts -= 1;
    }
    get_d_isog(R, A24minus, A24plus, coeff);    
    fp2add(A24plus, A24minus, A);                 
    fp2add(A, A, A);
    fp2sub(A24plus, A24minus, A24plus);                   
    j_inv(A, A24plus, jinv);
    fp2_encode(jinv, SharedSecretB);    // Format shared secret

    return 0;
}

#else

int EphemeralKeyGeneration_B(const unsigned char* PrivateKeyB, unsigned char* PublicKeyB)
{ // Bob's ephemeral public key generation
  // Input:  a private key PrivateKeyB in the range [0, 2^Floor(Log(2,oB)) - 1]. 
  // Output: the public key PublicKeyB consisting of 3 elements in GF(p^2) which are encoded by removing leading 0 bytes.
    point_proj_t R1, public_pts[4] = {0}, pts1[MAX_INT_POINTS_BOB1], pts2[MAX_INT_POINTS_BOB2], KernelPoints2[KERNEL_POINTS2], test;
    f2elm_t  A24plus = {0}, A24minus = {0}, A = {0}, C24 = {0}, XPB1, XQB1, XRB1, XPB2, XQB2, XRB2 ;
    unsigned int i, row, m, index = 0, pts_index1[MAX_INT_POINTS_BOB1], pts_index2[MAX_INT_POINTS_BOB2], npts = 0, ii = 0;                 
    #if defined PARALLEL
        int tid, cn;
        cn = omp_get_max_threads();
    #endif
    
    // Initialize basis points
    init_basis((digit_t*)A_gen, public_pts[0]->X, public_pts[1]->X, public_pts[2]->X);
    init_basis((digit_t*)B1_gen, XPB1, XQB1, XRB1);
    init_basis((digit_t*)B2_gen, XPB2, XQB2, XRB2);
    
    
    fpcopy((digit_t*)&Montgomery_one, (public_pts[0]->Z)[0]);
    fpcopy((digit_t*)&Montgomery_one, (public_pts[1]->Z)[0]);
    fpcopy((digit_t*)&Montgomery_one, (public_pts[2]->Z)[0]);

    // Initialize constants
    fpcopy((digit_t*)&Montgomery_one, A24plus[0]);
    fp2add(A24plus, A24plus, A24plus);
    fp2copy(A24plus, A24minus);
    fp2neg(A24minus);

#if defined PARALLEL

    // Retrieve kernel point
    #pragma omp parallel num_threads(2)
    {
        int i = omp_get_thread_num();
         if (i == 0) LADDER3PT(XPB1, XQB1, XRB1, (digit_t*)PrivateKeyB, BOB, R1, A);
         if (i == 1) LADDER3PT(XPB2, XQB2, XRB2, &((digit_t*)PrivateKeyB)[Index_Bob1], BOB1, public_pts[3], A);
    }

#else

    // Retrieve kernel point
    LADDER3PT(XPB1, XQB1, XRB1, (digit_t*)PrivateKeyB, BOB, R1, A);    
    LADDER3PT(XPB2, XQB2, XRB2, &((digit_t*)PrivateKeyB)[Index_Bob1], BOB1, public_pts[3], A);
   
#endif
    
    // Traverse tree for 3 isogenies===========
    index = 0; 
    fp2copy(R1->X, pts1[0]->X);
    fp2copy(R1->Z, pts1[0]->Z);
    pts_index1[0] = 0;
    npts = 1;
#if (KERNEL_POINTS1 == 1)
    f2elm_t coeff[3];
    for (row = 0; row < MAX_Bob1 -1 ; row++) {
        while (index < MAX_Bob1 - 1 - row) {
            m = strat_Bob1[ii];
            ii +=1;
            index += m;
            pts_index1[npts] = index;
            xTPLe(R1, R1, A24minus, A24plus, (int)m);
            fp2copy(R1->X, pts1[npts]->X);
            fp2copy(R1->Z, pts1[npts]->Z);
            npts +=1;
        }
        get_3_isog(R1, A24minus, A24plus, coeff);

    #if defined PARALLEL

        #pragma omp parallel for num_threads(cn)
        for (i = 0; i < (npts+3); i++) {
            if(i < npts - 1)
                eval_3_isog(pts1[i], coeff);
            else
                eval_3_isog(public_pts[i%4], coeff);
        }

    #else

        for (i = 0; i < npts-1; i++) {
            eval_3_isog(pts1[i], coeff);
        }     
        eval_3_isog(public_pts[0], coeff);
        eval_3_isog(public_pts[1], coeff);
        eval_3_isog(public_pts[2], coeff);
        eval_3_isog(public_pts[3], coeff);
        
    #endif        

        fp2copy(pts1[npts-2]->X, R1->X); 
        fp2copy(pts1[npts-2]->Z, R1->Z);
        index = pts_index1[npts-2];
        npts -= 1;
    }
   
    get_3_isog(R1, A24minus, A24plus, coeff);
    fp2sub(A24plus, A24minus, C24);  

    #if defined PARALLEL

    #pragma omp parallel for num_threads(3)
    for (i = 0; i < 4; i++) {
        eval_3_isog(public_pts[i], coeff);
    } 

    #else

    eval_3_isog(public_pts[0], coeff);
    eval_3_isog(public_pts[1], coeff);
    eval_3_isog(public_pts[2], coeff);
    eval_3_isog(public_pts[3], coeff);

    #endif    
#else
    //Traverse the tree for d1 isogeny
    point_proj_t KernelPoints1[KERNEL_POINTS1];
    fp2sub(A24plus, A24minus, C24); 
    for (row = 1; row < MAX_Bob1; row++) {
        while (index < MAX_Bob1 - row) {
            fp2copy(R1->X, pts1[npts]->X);
            fp2copy(R1->Z, pts1[npts]->Z);
            pts_index1[npts] = index;
            npts += 1;
            m = splits_Bob1[MAX_Bob1 - index - row];
            xMULe1(R1, R1, A24minus, A24plus, C24, (int)m);
            index += m;
        }
        get_d1_isog(R1, KernelPoints1, A24minus, A24plus, C24);

    #if defined PARALLEL

        #pragma omp parallel for num_threads(3)
        for (i = 0; i < (npts+4); i++) {
            if(i < npts)
                eval_d1_isog(KernelPoints1, pts1[i]);
            else
                eval_d1_isog(KernelPoints1, public_pts[i%4];
        }

    #else

        for (i = 0; i < npts; i++) {
            eval_d1_isog(KernelPoints1, pts1[i];
        }     
        eval_d1_isog(KernelPoints1,public_pts[0]);
        eval_d1_isog(KernelPoints1,public_pts[1]);
        eval_d1_isog(KernelPoints1,public_pts[2]);
        eval_d1_isog(KernelPoints1,public_pts[3]);
        
    #endif        

        fp2copy(pts1[npts-1]->X, R1->X); 
        fp2copy(pts1[npts-1]->Z, R1->Z);
        index = pts_index1[npts-1];
        npts -= 1;
    }
   
    get_d1_isog(R1, KernelPoints1, A24minus, A24plus, C24); 

    #if defined PARALLEL

    #pragma omp parallel for num_threads(cn)
    for (i = 0; i < 4; i++) {
        eval_d1_isog(KernelPoints1, public_pts[i]);
    } 

    #else

    eval_d1_isog(KernelPoints1, public_pts[0]);
    eval_d1_isog(KernelPoints1, public_pts[1]);
    eval_d1_isog(KernelPoints1, public_pts[2]);
    eval_d1_isog(KernelPoints1, public_pts[3]);
    #endif

#endif
    // Traverse tree for d2 isogeny//=================
    index = 0; 
    fp2copy(public_pts[3]->X, pts2[0]->X);
    fp2copy(public_pts[3]->Z, pts2[0]->Z);
    pts_index2[0] = 0;
    npts = 1;
    ii=0;
    for (row = 0; row < MAX_Bob2-1; row++) {
        while (index < MAX_Bob2 - 1 - row) {
            m = strat_Bob2[ii];
            ii +=1;
            index += m;
            pts_index2[npts] = index;
            xMULe2(public_pts[3], public_pts[3], A24plus, A24minus, C24, (int)m);
            fp2copy(public_pts[3]->X, pts2[npts]->X);
            fp2copy(public_pts[3]->Z, pts2[npts]->Z);
            npts += 1;
        }
        get_d2_isog(public_pts[3], KernelPoints2, A24plus, A24minus, C24);

    #if defined PARALLEL
 
        #pragma omp parallel for num_threads(cn)
        for (i = 0; i < (npts + 2); i++) {
            if(i < npts - 1)
                eval_d2_isog(KernelPoints2, pts2[i]);
            else
                eval_d2_isog(KernelPoints2, public_pts[i%3]);
        }
       
    #else        
        
        for (i = 0; i < npts - 1; i++) {
            eval_d2_isog(KernelPoints2, pts2[i]);
        }     
        eval_d2_isog(KernelPoints2, public_pts[0]);
        eval_d2_isog(KernelPoints2, public_pts[1]);
        eval_d2_isog(KernelPoints2, public_pts[2]);

    #endif

        fp2copy(pts2[npts-2]->X, public_pts[3]->X); 
        fp2copy(pts2[npts-2]->Z, public_pts[3]->Z);
        index = pts_index2[npts-2];
        npts -= 1;
    }

    get_d2_isog(public_pts[3], KernelPoints2, A24plus, A24minus, C24);

    #if defined PARALLEL

    #pragma omp parallel for num_threads(cn)
    for(i = 0; i < 3; i++){
        eval_d2_isog(KernelPoints2, public_pts[i]);
    }

    #else

    eval_d2_isog(KernelPoints2, public_pts[0]);
    eval_d2_isog(KernelPoints2, public_pts[1]);
    eval_d2_isog(KernelPoints2, public_pts[2]);

    #endif    

    inv_3_way(public_pts[0]->Z, public_pts[1]->Z, public_pts[2]->Z);
    fp2mul_mont(public_pts[0]->X, public_pts[0]->Z, public_pts[0]->X);
    fp2mul_mont(public_pts[1]->X, public_pts[1]->Z, public_pts[1]->X);
    fp2mul_mont(public_pts[2]->X, public_pts[2]->Z, public_pts[2]->X);

    // Format public key
    fp2_encode(public_pts[0]->X, PublicKeyB);
    fp2_encode(public_pts[1]->X, PublicKeyB + FP2_ENCODED_BYTES);    
    fp2_encode(public_pts[2]->X, PublicKeyB + 2*FP2_ENCODED_BYTES);

    return 0;
}

int EphemeralSecretAgreement_B(const unsigned char* PrivateKeyB, const unsigned char* PublicKeyA, unsigned char* SharedSecretB)
{ // Bob's ephemeral shared secret computation
  // It produces a shared secret key SharedSecretB using his secret key PrivateKeyB and Alice's public key PublicKeyA
  // Inputs: Bob's PrivateKeyB is an integer in the range [0, 2^Floor(Log(2,oB)) - 1]. 
  //         Alice's PublicKeyA consists of 3 elements in GF(p^2) encoded by removing leading 0 bytes.
  // Output: a shared secret SharedSecretB that consists of one element in GF(p^2) encoded by removing leading 0 bytes.  
    point_proj_t R1, R2,  pts1[MAX_INT_POINTS_BOB1], pts2[MAX_INT_POINTS_BOB2], KernelPoints2[KERNEL_POINTS2];
    f2elm_t PKB[3], jinv;
    f2elm_t A, A24plus = {0}, A24minus = {0}, C24={0};
    unsigned int i, row, m, index = 0, pts_index1[MAX_INT_POINTS_BOB1], pts_index2[MAX_INT_POINTS_BOB2], npts = 0, ii = 0;
    #if defined PARALLEL
        int tid, cn;
        cn = omp_get_max_threads();
    #endif
    // Initialize images of Alice's basis
    fp2_decode(PublicKeyA, PKB[0]);
    fp2_decode(PublicKeyA + FP2_ENCODED_BYTES, PKB[1]);
    fp2_decode(PublicKeyA + 2*FP2_ENCODED_BYTES, PKB[2]);

    // Initialize constants
    get_A(PKB[0], PKB[1], PKB[2], A); // TODO: Can return projective A?

    fpadd((digit_t*)&Montgomery_one, (digit_t*)&Montgomery_one, A24minus[0]);
    fp2add(A, A24minus, A24plus);
    fp2sub(A, A24minus, A24minus);
    fp2sub(A24plus, A24minus, C24);

#if defined PARALLEL

    // Retrieve kernel point
    #pragma omp parallel num_threads(2)
    {
        int i = omp_get_thread_num();
        if (i == 0){
            LADDER3PT(PKB[0], PKB[1], PKB[2], (digit_t*)PrivateKeyB, BOB, R1, A);
            xMULe2(R1, R1, A24plus, A24minus, C24, MAX_Bob2);
        }
        if (i == 1){
            LADDER3PT(PKB[0], PKB[1], PKB[2], &((digit_t*)PrivateKeyB)[Index_Bob1], BOB1, R2, A);
            #if (KERNEL_POINTS1 == 1)
                xTPLe(R2, R2, A24minus, A24plus,  MAX_Bob1);
            #else
                xMULe1(R2, R2, A24minus, A24plus,  MAX_Bob1);
            #endif
        }
    }

#else

    // Retrieve kernel point
    LADDER3PT(PKB[0], PKB[1], PKB[2], (digit_t*)PrivateKeyB, BOB, R1, A);
    xMULe2(R1, R1, A24plus, A24minus, C24, MAX_Bob2);
    LADDER3PT(PKB[0], PKB[1], PKB[2], &((digit_t*)PrivateKeyB)[Index_Bob1], BOB1, R2, A);
    #if (KERNEL_POINTS1 == 1)
        xTPLe(R2, R2, A24minus, A24plus,  MAX_Bob1);
    #else
        xMULe1(R2, R2, A24minus, A24plus,  MAX_Bob1);
    #endif

#endif
 
    // Traverse tree for 3-isogeny
    index = 0; 
    fp2copy(R1->X, pts1[0]->X);
    fp2copy(R1->Z, pts1[0]->Z);
    pts_index1[0] = 0;
    npts = 1;
#if (KERNEL_POINTS1 == 1)
    f2elm_t coeff[3];
    for (row = 0; row < MAX_Bob1 -1 ; row++) {
        while (index < MAX_Bob1 - 1 - row) {
            m = strat_Bob1[ii];
            ii +=1;
            index += m;
            pts_index1[npts] = index;
            xTPLe(R1, R1, A24minus, A24plus, (int)m);
            fp2copy(R1->X, pts1[npts]->X);
            fp2copy(R1->Z, pts1[npts]->Z);
            npts +=1;
        }
        get_3_isog(R1, A24minus, A24plus, coeff);

    #if defined PARALLEL

        #pragma omp parallel for num_threads(cn)
        for (i = 0; i < npts; i++) {
            if(i < npts - 1)
                eval_3_isog(pts1[i], coeff);
            else
                eval_3_isog(R2, coeff);
        }

    #else

        for (i = 0; i < npts-1; i++) {
            eval_3_isog(pts1[i], coeff);
        }
        eval_3_isog(R2, coeff);

    #endif

        fp2copy(pts1[npts-2]->X, R1->X); 
        fp2copy(pts1[npts-2]->Z, R1->Z);
        index = pts_index1[npts-2];
        npts -= 1;
    }
    
    get_3_isog(R1, A24minus, A24plus, coeff);
    fp2sub(A24plus, A24minus, C24);
    eval_3_isog(R2, coeff);
//end 3-isogeny
#else
    // Traverse tree for d1 isogeny=======================
    point_proj_t KernelPoints2[KERNEL_POINTS1];
    index = 0; npts = 0;
    for (row = 1; row < MAX_Bob1; row++) {
        while (index < MAX_Bob1 - row) {
            fp2copy(R1->X, pts1[npts]->X);
            fp2copy(R1->Z, pts1[npts]->Z);
            pts_index1[npts] = index;
            npts += 1;
            m = splits_Bob1[MAX_Bob1 - index - row];
            xMULe1(R1, R1, A24minus, A24plus, (int)m);
            index += m;
        }
        get_d1_isog(R1, KernelPoints1,  A24minus, A24plus, C24);

    #if defined PARALLEL

        #pragma omp parallel for num_threads(3)
        for (i = 0; i < (npts+1); i++) {
            if(i < npts)
                eval_d1_isog(KernelPoints1, pts3[i]);
            else
                eval_d1_isog(KernelPoints1, R2);
        }

    #else

        for (i = 0; i < npts; i++) {
            eval_d1_isog(KernelPoints1, pts3[i]);
        }
        eval_d1_isog(KernelPoints1, R2);

    #endif

        fp2copy(pts1[npts-1]->X, R1->X); 
        fp2copy(pts1[npts-1]->Z, R1->Z);
        index = pts_index1[npts-1];
        npts -= 1;
    }
    
    get_d1_isog(R3, KernelPoints1, A24minus, A24plus, C24);
    eval_d1_isog(KernelPoints1, R2);
#endif
    
    // Traverse tree for d2 isogeny=======================
    index = 0;
    fp2copy(R2->X, pts2[0]->X);
    fp2copy(R2->Z, pts2[0]->Z);
    pts_index2[0] = 0;
    npts = 1;
    ii=0;
    for (row = 0; row < MAX_Bob2-1; row++) {
        while (index < MAX_Bob2 - 1 - row) {
            m = strat_Bob2[ii];
            ii +=1;
            index += m;
            pts_index2[npts] = index;
            xMULe2(R2, R2, A24plus, A24minus, C24, (int)m);
            fp2copy(R2->X, pts2[npts]->X);
            fp2copy(R2->Z, pts2[npts]->Z);
            npts += 1;
        }
        get_d2_isog(R2, KernelPoints2, A24plus, A24minus, C24);

    #if defined PARALLEL

        #pragma omp parallel for num_threads(cn)
        for (i = 0; i < npts-1; i++) {
            eval_d2_isog(KernelPoints2, pts2[i]);
        }

    #else
        
        for (i = 0; i < npts-1; i++) {
            eval_d2_isog(KernelPoints2, pts2[i]);
        }     

    #endif
        
        fp2copy(pts2[npts-2]->X, R2->X); 
        fp2copy(pts2[npts-2]->Z, R2->Z);
        index = pts_index2[npts-2];
        npts -= 1;
    }
    
    get_d2_isog(R2, KernelPoints2, A24plus, A24minus, C24);
    fp2add(A24plus, A24minus, A24plus);                 
    fp2add(A24plus, A24plus, A24plus);
    j_inv(A24plus, C24, jinv);
    
    fp2_encode(jinv, SharedSecretB);    // Format shared secret

    return 0;
}
#endif
