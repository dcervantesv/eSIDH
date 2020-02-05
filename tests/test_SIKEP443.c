/********************************************************************************************
* SIDH: an efficient supersingular isogeny cryptography library
*
* Abstract: benchmarking/testing isogeny-based key encapsulation mechanism SIKEp443
*********************************************************************************************/ 

#include <stdio.h>
#include <string.h>
#include "test_extras.h"
#include "P443_api.h"


#define SCHEME_NAME    "SIKEp443"

#define crypto_kem_keypair            crypto_kem_keypair_SIKEp443
#define crypto_kem_enc                crypto_kem_enc_SIKEp443
#define crypto_kem_dec                crypto_kem_dec_SIKEp443

#include "test_sike.c"
