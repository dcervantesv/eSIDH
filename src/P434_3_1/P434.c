/********************************************************************************************
* SIDH: an efficient supersingular isogeny cryptography library
*
* Abstract: supersingular isogeny parameters and generation of functions for P434
*********************************************************************************************/  

#include "P434_api.h" 
#include "P434_internal.h"


// Encoding of field elements, elements over Z_order, elements over GF(p^2) and elliptic curve points:
// --------------------------------------------------------------------------------------------------
// Elements over GF(p) and Z_order are encoded with the least significant octet (and digit) located at the leftmost position (i.e., little endian format). 
// Elements (a+b*i) over GF(p^2), where a and b are defined over GF(p), are encoded as {a, b}, with a in the least significant position.
// Elliptic curve points P = (x,y) are encoded as {x, y}, with x in the least significant position. 
// Internally, the number of digits used to represent all these elements is obtained by approximating the number of bits to the immediately greater multiple of 32.
// For example, a 434-bit field element is represented with Ceil(434 / 64) = 12 64-bit digits or Ceil(434 / 32) = 24 32-bit digits.

//
// Curve isogeny system "SIDHp434". Base curve: Montgomery curve By^2 = Cx^3 + Ax^2 + Cx defined over GF(p434^2), where A=0, B=1, C=1 and p434 = 2^372*3^239-1
//
         
const uint64_t p434[NWORDS64_FIELD]              = { 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 
0xFDC1767AE2FFFFFF, 0x7BC65C783158AEA3, 0x6CFC5FD681C52056, 0x2341F27177344 };
const uint64_t p434p1[NWORDS64_FIELD]            = { 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0xFDC1767AE3000000, 0x7BC65C783158AEA3, 0x6CFC5FD681C52056, 0x2341F27177344 };
const uint64_t p434x2[NWORDS64_FIELD]            = { 0xFFFFFFFFFFFFFFFE, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 
0xFB82ECF5C5FFFFFF, 0xF78CB8F062B15D47, 0xD9F8BFAD038A40AC, 0x4683E4E2EE688 }; 
// Order of Alice's subgroup
const uint64_t Alice_order[NWORDS64_ORDER]       = { 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x1000000 }; //2^216
// Order of Bob's subgroup
const uint64_t Bob_order[NWORDS64_ORDER]         = { 0x58AEA3FDC1767AE3, 0xC520567BC65C7831, 0x1773446CFC5FD681, 0x2341F27 }; //3^137


// Alice's generator values {XPA0 + XPA1*i, XQA0, XRA0 + XRA1*i} in GF(p434^2), expressed in Montgomery representation
const uint64_t A_gen[5 * NWORDS64_FIELD]  = {
0xC148345E1A874621,0xAD0DF58AD5174982,0xBC031436B9826478,0x850B0F316FABD4BC,0xC604F2F140209357,0xA1418A010704614C,0x19D97C9ECCA74,
0x359F66A6431A8CEB,0xAC4F1437446CDAEB,0x65E79636ECC2E802,0xEC9F349E034E0EEF,0x8E22FF1B719ABC6E,0xDCAE79C21BA04037,0x2179AC860DF03,
0xBA6B4C181E2E1385,0x29C89ED743125AB3,0x74C9798B9E8CC934,0x5F0A4BBDF56C3911,0x5813CDFC65DC9483,0x6ECC1A9B0858C876,0x114445B2CF18A,
0x9FCB26941E81BEFD,0xCE0C786411E4E7AA,0x5B167FAAD7DE9A23,0xFE4264B07B429D07,0x80137E77C875E284,0x34D2F18743CC2975,0x309C355602D4,
0x96C8D7019584A4E7,0x669ABCE8ADD32E5,0x7D575D74BF456ABA,0x3864E0639C211414,0x25D7DC021A1CEFE1,0x1B5016D98CE59116,0x194EC745F6589};

// Bob's generator values {XPB0 + XPB1*i, XQB0, XRB0 + XRB1*i} in GF(p434^2), expressed in Montgomery representation
const uint64_t B_gen[5 * NWORDS64_FIELD]  = {
0x1EC9D18CED36993C,0xB986E81FE5A4AE16,0xB9ADA1757120E43C,0xC33FC934711893C2,0x7C5D154E1C436F02,0xA43320196C86CC2,0x13BC8C60CE2EB,
0x5CAFE1733050A9D6,0x762BFD61928C1086,0x1457C70149DA64E2,0xB0B5A28415AA8BA4,0xCAF0A978B415C385,0xB5B3618CFAB0D756,0x37F1CA961BAA,
0x11980E8529758837,0x6C818184E5D5F21,0xF5D8D6546B186E4B,0xAA447B20F909DBA4,0x3EF8A54DE1679A86,0xD0E63F5B36770B07,0x15791AE6503A6,
0xDA66E3B3A5EF5947,0x167B48226AC3D10B,0x6AB8BBA82D6B7A05,0xE5A10D2C72F02235,0xA2BDF72CC483D280,0x202D6BEE81FA39DC,0x118CFA82AE3C3,
0xAE58A2A52525BB38,0xCA1AA85583D40792,0xA1F262A040764AC,0xF878515B459DF066,0xF18B001580264060,0x62EA6E715FE3B20C,0x12851D0F56390};


// Montgomery constant Montgomery_R2 = (2^768)^2 mod p434
const uint64_t Montgomery_R2[NWORDS64_FIELD]     = { 0x28E55B65DCD69B30, 0xACEC7367768798C2, 0xAB27973F8311688D, 0x175CC6AF8D6C7C0B, 0xABCD92BF2DDE347E, 0x69E16A61C7686D9A, 0x000025A89BCDD12A };                                                    
// Value one in Montgomery representation 
const uint64_t Montgomery_one[NWORDS64_FIELD]    = { 0x000000000000742C, 0x0000000000000000, 0x0000000000000000, 0xB90FF404FC000000, 0xD801A4FB559FACD4, 0xE93254545F77410C, 0x0000ECEEA7BD2EDA };     
// Value (2^384)^2 mod 3^239                                                   
const uint64_t Montgomery_Rprime[NWORDS64_ORDER] = { 0x1A55482318541298, 0x070A6370DFA12A03, 0xCB1658E0E3823A40, 0xB3B7384EB5DEF3F9 };
// Value -(3^239)^-1 mod 2^384 
const uint64_t Montgomery_rprime[NWORDS64_ORDER] = { 0x48062A91D3AB563D, 0x6CE572434303C2F5, 0x5D1319F3F160EC9D, 0xE35554E8C2D5623A };                                           


// Fixed parameters for isogeny tree computation
const unsigned int strat_Alice[MAX_Alice-1] = 
 {48, 27, 15, 8, 4, 2, 2, 1, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 21, 12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1,1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1 };

const unsigned int strat_Bob[MAX_Bob-1] = 
{49, 33, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 
1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 12, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 4, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 16, 12, 8, 5, 3, 
2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 4, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 5, 3, 3, 2, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1};

// Setting up macro defines and including GF(p), GF(p^2), curve, isogeny and kex functions
#define fpcopy                        fpcopy434
#define fpzero                        fpzero434
#define fpadd                         fpadd434
#define fpsub                         fpsub434
#define fpneg                         fpneg434
#define fpdiv2                        fpdiv2_434
#define sqr_asm                       sqr434_asm
#define fpcorrection                  fpcorrection434
#define fpmul_mont                    fpmul434_mont
#define fpsqr_mont                    fpsqr434_mont
#define fpinv_mont                    fpinv434_mont
#define fpinv_chain_mont              fpinv434_chain_mont
#define fpinv_mont_bingcd             fpinv434_mont_bingcd
#define fp2copy                       fp2copy434
#define fp2zero                       fp2zero434
#define fp2add                        fp2add434
#define fp2sub                        fp2sub434
#define fp2neg                        fp2neg434
#define fp2div2                       fp2div2_434
#define fp2correction                 fp2correction434
#define fp2mul_mont                   fp2mul434_mont
#define fp2sqr_mont                   fp2sqr434_mont
#define fp2inv_mont                   fp2inv434_mont
#define fp2inv_mont_bingcd            fp2inv434_mont_bingcd
#define fpequal_non_constant_time     fpequal434_non_constant_time
#define mp_add_asm                    mp_add434_asm
#define mp_subx2_asm                  mp_sub434x2_asm
#define mp_dblsubx2_asm               mp_dblsub434x2_asm
#define crypto_kem_keypair            crypto_kem_keypair_SIKEp434
#define crypto_kem_enc                crypto_kem_enc_SIKEp434
#define crypto_kem_dec                crypto_kem_dec_SIKEp434
#define random_mod_order_A            random_mod_order_A_SIDHp434
#define random_mod_order_B            random_mod_order_B_SIDHp434
#define EphemeralKeyGeneration_A      EphemeralKeyGeneration_A_SIDHp434
#define EphemeralKeyGeneration_B      EphemeralKeyGeneration_B_SIDHp434
#define EphemeralSecretAgreement_A    EphemeralSecretAgreement_A_SIDHp434
#define EphemeralSecretAgreement_B    EphemeralSecretAgreement_B_SIDHp434

//Setting up xMUL and Isogeny functions
#define xMUL1e                          xTPLe
#define get_d_isog                      get_3_isog
#define eval_d_isog                     eval_3_isog

#include "../fpx.c"
#include "../ec_isogeny.c"
#include "../sidh.c"
#include "../sike.c"
