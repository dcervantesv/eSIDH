# eSIDH v2.0 (C version)

This code is complementary to the the IACR eprint paper [eSIDH: the revenge of the SIDH](https://ia.cr/2020/021)   
and the paper [Parallel strategies for SIDH: Towards computing SIDH twice as fast](https://eprint.iacr.org/2020/383)  

The extended SIDH (eSIDH) is a proposal by : 

- Daniel Cervantes-Vázquez (dcervantes@computacion.cs.cinvestav.mx)  
- Eduardo Ochoa-Jiménez  
- Francisco Rodríguez-Henríquez  
***
## Introduction
 The eSIDH is a small variant of the protocol [SIDH](https://sike.org/) 
 which make use of two small primes for Bob isogeny
 computations instead of the original setting of only one prime.
 
 This new kind of primes allows us (in some cases) to find 
 better Montgomery-friendly primes than the original SIDH 
 
 This new eSIDH like primes takes advantage of parallel computing by the use
 of two private keys for Bob and by the fact that both can be computed independently
 of the other. As both keys have the half of the bitsize of the Bob's SIDH 
 private keys then we expect an acceleration on the protocol. 
 Moreover as eSIDH-like primes are expected to be "more" Montgomery Friendly
 than the SIDH primes, then we expect an acceleration on the Field and 
 Quadratic Field arithmetic.
 Both features makes eSIDH a competitive version of SIDH.
***
## Features of this version


This version includes the eSIDH primes

- p_{443} = 2^{222} 3^{73} 5^{45}
- p_{765} = 2^{391} 3^{119} 5^{81}

and the SIDH primes

- p_{434} = 2^{216} 3^{137}
- p_{751} = 2^{372} 3^{239}

This library allows the creation of test files for SIDH and SIKE protocols
using SIDH and eSIDH primes. Also allows to test the Field and Quadratic Field
arithmetic.

For now on, only the Parallel version of eSIDH is available.

### New in this version

- We include the parallel computing of Alice's secret key and a portion of the strategy to compute her 4^ea isogeny.
- We include the use of suitable strategies for this task and to compute strategies in parallel (software oriented).
- We include the Montgomery  ladder for fixed point introduced in [A Faster Software Implementation of the Supersingular Isogeny Diffie-Hellman Key Exchange Protocol](https://eprint.iacr.org/2017/1015) 
- We include magma code of algorithms to compute parallel strategies and the new strategy including the computation in parallel of Alice's secret point. Available in Magma directory. 

We only test this library under a x64-Linux distribution and a compatible
compiler supporting the [Open MP api](https://www.openmp.org/).

>Our library is based on [Microsoft SIDH library](https://github.com/Microsoft/PQCrypto-SIDH)

***
## Compilation and Usage
To create test files for Field and Quadratic Field arithmetic, (single core) SIDH and SIKE
protocols run the following.
>`make PRIME=prime ARCH=x64 CC=compiler OPT_LEVEL=FAST USE_MULX=TRUE/FALSE USE_ADX=TRUE/FALSE SET=EXTENDED`

where `prime` could be one of the following:

- `P443_3_5_1`
- `P765_3_5_1`
- `P751_3_1`
- `P434_3_1`

This create the executable files 

- `arith_tests-prime`
- `sidhprime/test_SIDH`
- `sikeprime/test_SIKE`

To create executble files for the parallel version of eSIDH it is
necessary to add the compilation flag `MODE=PARALLEL` and set the
linux environment variable `OMP_NUM_THREADS` with the desire number of threads

In version 2.0 we include the compiler flag `MODE=NEWSTRAT` which allows the use of tricks presented
in `Parallel strategies for SIDH: Towards computing SIDH twice as fast (to be published)`.

In this version, the user must manually selects the strategy to be used (as strategies for PARALLEL and NEWSTRAT are different).
This can be achieved as follows.
Go to the path `src/prime`, where prime is as above. Then inside the file `prime_short.c` where prime_short 
could be `p751, p765, p434` or `p443` you can find the variable `strat_Alice` and the options (commented) `PARALLEL n cores`
and `NEWSTRAT n cores`. By default the option `NEWSTRAT 2 cores` is uncommented.
Now for Bob there are differences,

- for SIKE primes  there is only the variable `strat_Bob` and the options  `PARALLEL n cores`.  
- For eSIDH primes there are two variables, namely `strat_Bob1` and `strat_Bob2`, again only the option `PARALLEL n cores` is available.

In both cases, the option `PARALLEL 2 cores` is uncommented by default.

#### Examples
We present an example of how to compile parallel eSIDH using two cores

To compile the parallel eSIDH  using the prime p_{765}, you must run the following 
> `make PRIME=P765_3_5_1 ARCH=x64 CC=clang OPT_LEVEL=FAST USE_MULX=TRUE USE_ADX=TRUE SET=EXTENDED MODE=NEWSTRAT`

Before running SIDH or SIKE protocols set the variable
> `export OMP_NUM_THREADS=2`

Then you can run
> `./sidhp765/test_SIDH`
***

## Ongoing work for the next version

We are working on: 

- To have the full eSIDH primes roster.
- Add the CRT version of eSIDH
- To explore alternatives to Posix threads and Open MP.
- To include all SIKE primes.
- To configure strategies via compiler flag instead of manually select one.
***
