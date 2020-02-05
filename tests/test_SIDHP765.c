/********************************************************************************************
* SIDH: an efficient supersingular isogeny cryptography library
*
* Abstract: benchmarking/testing isogeny-based key exchange SIDHp765
*********************************************************************************************/ 

#include <stdio.h>
#include <string.h>
#include "test_extras.h"
#include "P765_api.h"


#define SCHEME_NAME    "SIDHp765"

#define random_mod_order_A            random_mod_order_A_SIDHp765
#define random_mod_order_B            random_mod_order_B_SIDHp765
#define EphemeralKeyGeneration_A      EphemeralKeyGeneration_A_SIDHp765
#define EphemeralKeyGeneration_B      EphemeralKeyGeneration_B_SIDHp765
#define EphemeralSecretAgreement_A    EphemeralSecretAgreement_A_SIDHp765
#define EphemeralSecretAgreement_B    EphemeralSecretAgreement_B_SIDHp765

#include "test_sidh.c"
