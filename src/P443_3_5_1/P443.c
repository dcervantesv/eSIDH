/********************************************************************************************
* SIDH: an efficient supersingular isogeny cryptography library
*
* Abstract: supersingular isogeny parameters and generation of functions for P443
*********************************************************************************************/  

#include "P443_api.h" 
#include "P443_internal.h"


// Encoding of field elements, elements over Z_order, elements over GF(p^2) and elliptic curve points:
// --------------------------------------------------------------------------------------------------
// Elements over GF(p) and Z_order are encoded with the least significant octet (and digit) located at the leftmost position (i.e., little endian format). 
// Elements (a+b*i) over GF(p^2), where a and b are defined over GF(p), are encoded as {a, b}, with a in the least significant position.
// Elliptic curve points P = (x,y) are encoded as {x, y}, with x in the least significant position. 
// Internally, the number of digits used to represent all these elements is obtained by approximating the number of bits to the immediately greater multiple of 32.
// For example, a 443-bit field element is represented with Ceil(443 / 64) = 12 64-bit digits or Ceil(443 / 32) = 24 32-bit digits.

//
// Curve isogeny system "SIDHp443". Base curve: Montgomery curve By^2 = Cx^3 + Ax^2 + Cx defined over GF(p443^2), where A=0, B=1, C=1 and p443 = 2^222*3^73*5^45-1
//
         
const uint64_t p443[NWORDS64_FIELD]              = { 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0x1D236727BFFFFFFF, 0xD9EBCCB301480652, 0xEA6695D5C29E0A13, 0x048F5AB5AA19ECE9 };
const uint64_t p443p1[NWORDS64_FIELD]            = { 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x1D236727C0000000, 0xD9EBCCB301480652, 0xEA6695D5C29E0A13, 0x048F5AB5AA19ECE9 };
const uint64_t p443x2[NWORDS64_FIELD]            = { 0xFFFFFFFFFFFFFFFE, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0x3A46CE4F7FFFFFFF, 0xB3D7996602900CA4, 0xD4CD2BAB853C1427, 0x91EB56B5433D9D3 }; 
// Order of Alice's subgroup
const uint64_t Alice_order[NWORDS64_ORDER]       = { 0x0, 0x0, 0x0, 0x1000000 }; 
// Order of Bob's subgroup
const uint64_t Bob_order[NWORDS64_ORDER]         = { 0xC968549F878A8EEB, 0x59B1A13F7CC76E3E, 0xE9867D6EBE876DA9, 0x2B5045CB25748 };

// Alice's generator values {XPA0 + XPA1*i, XQA0, XRA0 + XRA1*i} in GF(p443^2), expressed in Montgomery representation
const uint64_t A_gen[5 * NWORDS64_FIELD]  = {
0x5A6F62AE4DFFF457,0x53952DCB6814C58E,0xFE8AEE90532F5CF,0xF4628EAF40F1D020,0xD0369579AE909E96,0xE4F4C6AD4BB28844,0x91E94CB1B22B1,
0x160190D71E70337F,0x902357A514FB0795,0x5B0AEA6CBF099080,0x76C5BD5E900DCEC2,0x5DC59E0C55532687,0x629C04E85AD91E2F,0x461DC70307B0C05,
0xE5F95ADEC3B7A01D,0x19DA0C226A18965B,0xF12F5A816B2234BE,0x2E3382A17A684244,0x7555AC858FC3C443,0xC79B52DEC0254C,0x1BFB5C7C0C0AF24,
0x86149463E5CB3F06,0xBF72EF3925EDFBD1,0x2B7EEE4F85552037,0x11624318C4212944,0xC919BBF83F4853B6,0xAA70241409AFB77A,0x13A44EC67FE144B,
0x47C3970EC2715448,0xD2858D3867F0DB14,0xB8D0843E38F3B072,0xB78BF1AC386A8136,0x5F33B50C2CF61CB3,0xC249510B27EAA568,0x279A9EBC4F941AC};

// Bob's generator values {XPB0 + XPB1*i, XQB0, XRB0 + XRB1*i} in GF(p443^2), expressed in Montgomery representation
const uint64_t B_gen[5 * NWORDS64_FIELD]  = {
0xEAF3596A49CBDD5D,0x89A4B44445C9CB04,0x857FD0D871E3AA87,0x2EF14721D0688E4A,0x6E5D9A1FD37C9CCE,0xAEA42315CF833BB5,0x1626C760E0DC2C,
0x3C7226DEBD26798B,0x74B3729B80EE024B,0xCBB5E6B631D7B99,0x50C79B256607672,0x8098F4D609C8ED11,0x97D65EC4DF40D9F3,0x1E6DC4F3655C185,
0x203224E95C0AE957,0xD41BC35E58279153,0x78194BE21B3DAB3C,0xD18C2B388D4C0EA2,0xA8328CDD79F29F36,0xB2E53D30EBBB4800,0x14067030510DA0F,
0xA8DDA43DBAA0B89F,0x271CAB0DDCBC2930,0x266CFB6CB3537242,0xF71C34E4A7B19BAC,0xD21760087883FD15,0x3DD2C0F05E42084,0x2EA2E618B045344,
0xBD34DBD8C012925A,0xAC453CB3CCF93247,0x25FB7B92E219A516,0x1B8C3F5852E96654,0xC536255FF992EC44,0x5E9D5F6D3FE50234,0x1350BBCADF8679B};

const uint64_t B1_gen[5 * NWORDS64_FIELD]  = {
0xE7B2A558C4DD3EDF,0xB4F2874A32199728,0xB4AE6A1267E5050A,0x4CDE5216DE1741FD,0xC79C859F327A4584,0xFA37600FAF7DB3A1,0x1A3A9153A40356A,
0xD7B9165DC73D2DC6,0xC4A59A37A047611D,0x6B37C054B65F529A,0xE3B431D7AACFD558,0x581A0B90B8B9CA0C,0x1881AB33D6326C47,0x2A54024CCC60FCB,
0x627C4BB7EF6E602E,0xBDEC0BDFEDB6CF80,0x110166FC7D22689,0xEC15685572850E32,0x30759A1542CF339B,0x17C531767820FB48,0x17E03A9742B0B16,
0x5AB04A447058171B,0x6DC4C7C8AA5C5599,0x91705E8FC66D3D69,0xB79B838CC39C3251,0xDF196E27FD9FA2DE,0x88D409E2F7AC5B96,0x39128697D91E288,
0x49ADA8FAA1F8E3DB,0x85CF01EFDEAD3AB1,0x29718ADA735D72FB,0x448CD1D945D1B704,0x61A75F5AD00DB0D,0x414BE3D07D2F19B8,0x5CE63584CD2AB7};

const uint64_t B2_gen[5 * NWORDS64_FIELD]  = {
0xB0B45EB141FB4D31,0x403DFD5E8DDBED96,0x62869E6D56E489C9,0xB55025E132839ED3,0x9151297C447439F4,0xCB4BD1A17ECC2F9E,0x2C81A5162A03BA9,
0xA4A05B6A181E08E4,0x117E6CF5D251CE42,0x10CEB1A9B02CD419,0x386B28B3EC186648,0xBA57185501AE699,0x3A50FBD72EFEC581,0x2E4EF8ED301DB7A,
0x3EF463DDC7E5C6EE,0x6EBAEB5C9AEB2D,0x5DCCD343ECC99AC9,0xC840251D18D39210,0x26A0DF85CD6D114B,0x8EF05D38077E152E,0x1AFBEEDD83D2777,
0xDD8364E2CBBC21C0,0xC77F9971E3258DFE,0x4136614B8CA9A391,0x4DBB9F2BEF96315E,0xAA76D5CEE4DADAB0,0x58E2C1FB4BE730B5,0x34E2F1ED2D2B68D,
0xDCA7C512C4BBA413,0x9FB739689F15C3FA,0x34BF1EE988790754,0x841A67EA9CD276A5,0x2261E191FA4E09CC,0x371EB192DA649177,0xCC05ED3D86BFF2};

// Montgomery constant Montgomery_R2 = (2^768)^2 mod p443
const uint64_t Montgomery_R2[NWORDS64_FIELD]     = { 0xD3D00918DD40DF7E, 0x8E24046C952A32AE, 0x86BE1D24DB247F54, 0x10D06D0DF64591D2, 0xA92C77329BD6D0D9, 0x780222C835EF3FC4, 0x46726A189CEEB7B };   

// Value one in Montgomery representation 
const uint64_t Montgomery_one[NWORDS64_FIELD]    = { 0x0000000000000038, 0x0000000000000000, 0x0000000000000000, 0xA0416F4E00000000, 0x546B38D7B83E9E09, 0xB98F393D6D6DCBA8, 0xA42842CA542CD4 };     
  
const unsigned int strat_Alice[MAX_Alice-1] = 
 {48, 27, 15, 8, 4, 2, 2, 1, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 21, 12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1,1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1 };

const unsigned int strat_Bob1[MAX_Bob1] = 
{22, 18, 12, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 4, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 6, 4, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 1, 2,
1, 1, 1, 3, 2, 1, 1, 1, 1, 1  };

const unsigned int strat_Bob2[MAX_Bob2] = 
{13, 10, 7, 5, 3, 2, 2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 2, 2, 1, 1, 1, 1, 3, 2, 2, 1, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1};

// Setting up macro defines and including GF(p), GF(p^2), curve, isogeny and kex functions
#define fpcopy                        fpcopy443
#define fpzero                        fpzero443
#define fpadd                         fpadd443
#define fpsub                         fpsub443
#define fpneg                         fpneg443
#define fpdiv2                        fpdiv2_443
#define fpcorrection                  fpcorrection443
#define fpmul_mont                    fpmul443_mont
#define fpsqr_mont                    fpsqr443_mont
#define fpinv_mont                    fpinv443_mont
#define fpinv_chain_mont              fpinv443_chain_mont
#define fpinv_mont_bingcd             fpinv443_mont_bingcd
#define fp2copy                       fp2copy443
#define fp2zero                       fp2zero443
#define fp2add                        fp2add443
#define fp2sub                        fp2sub443
#define fp2neg                        fp2neg443
#define fp2div2                       fp2div2_443
#define fp2correction                 fp2correction443
#define fp2mul_mont                   fp2mul443_mont
#define fp2sqr_mont                   fp2sqr443_mont
#define fp2inv_mont                   fp2inv443_mont
#define fp2inv_mont_bingcd            fp2inv443_mont_bingcd
#define fpequal_non_constant_time     fpequal443_non_constant_time
#define sqr_asm                       sqr443_asm
#define mp_add_asm                    mp_add443_asm
#define mp_subx2_asm                  mp_sub443x2_asm
#define mp_dblsubx2_asm               mp_dblsub443x2_asm
#define crypto_kem_keypair            crypto_kem_keypair_SIKEp443
#define crypto_kem_enc                crypto_kem_enc_SIKEp443
#define crypto_kem_dec                crypto_kem_dec_SIKEp443
#define random_mod_order_A            random_mod_order_A_SIDHp443
#define random_mod_order_B            random_mod_order_B_SIDHp443
#define EphemeralKeyGeneration_A      EphemeralKeyGeneration_A_SIDHp443
#define EphemeralKeyGeneration_B      EphemeralKeyGeneration_B_SIDHp443
#define EphemeralSecretAgreement_A    EphemeralSecretAgreement_A_SIDHp443
#define EphemeralSecretAgreement_B    EphemeralSecretAgreement_B_SIDHp443

//Setting up xMUL and Isogeny functions
#define xMULe2                        xQPLe
#define get_d2_isog                   get_5_isog
#define eval_d2_isog                  eval_5_isog
    
#include "../fpx.c"
#include "../ec_isogeny.c"
#include "../sidh.c"
#include "../sike.c"
