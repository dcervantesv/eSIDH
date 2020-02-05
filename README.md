# eSIDH v1.0


This code is complementary to the the IACR eprint paper [eSIDH: the revenge of the SIDH](https://ia.cr/2020/021) 

## Introduction
 The eSIDH is a small variant of the protocol SIDH 
 which make use of two small primes for Bob isogeny
 computations instead of the original setting of only one prime.
 
 This new kind of primes allows us (in some cases) to find 
 better Montgomery-friendly primes than the original SIDH 

## Features
This version includes the primes

- p_{443} = 4^{222} 3^{73} 5^{45}
- p_{765} = 4^{391} 3^{119} 5^{81}

## Compilation and Usage
To create test files for Field and Quadratic Field arithmetic, SIDH and SIKE
protocols run the following
```Shell
make PRIME=P443_3_5_1 ARCH=x64 CC=clang OPT_LEVEL=FAST USE_MULX=TRUE USE_ADX=TRUE SET=EXTENDED
```
