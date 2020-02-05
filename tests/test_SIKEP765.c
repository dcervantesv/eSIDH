/********************************************************************************************
* SIDH: an efficient supersingular isogeny cryptography library
*
* Abstract: benchmarking/testing isogeny-based key encapsulation mechanism SIKEp765
*********************************************************************************************/ 

#include <stdio.h>
#include <string.h>
#include "test_extras.h"
#include "P765_api.h"


#define SCHEME_NAME    "SIKEp765"

#define crypto_kem_keypair            crypto_kem_keypair_SIKEp765
#define crypto_kem_enc                crypto_kem_enc_SIKEp765
#define crypto_kem_dec                crypto_kem_dec_SIKEp765

#include "test_sike.c"
