/*************************** RANDOMC.H ***************** 2007-09-22 Agner Fog *
 *
 * This file contains class declarations and other definitions for the C++
 * library of uniform random number generators.
 *
 * Overview of classes:
 * ====================
 *
 * class CRandomMersenne:
 * Random number generator of type Mersenne twister.
 * Source file mersenne.cpp
 *
 *
 *
 * Member functions (methods):
 * ===========================
 *
 * All these classes have identical member functions:
 *
 * Constructor(uint32 seed):
 * The seed can be any integer. Usually the time is used as seed.
 * Executing a program twice with the same seed will give the same sequence of
 * random numbers. A different seed will give a different sequence.
 *
 * void RandomInit(uint32 seed);
 * Re-initializes the random number generator with a new seed.
 *
 * void RandomInitByArray(uint32 seeds[], int length);
 * In CRandomMersenne only: Use this function if you want to initialize with
 * a seed with more than 32 bits. All bits in the seeds[] array will influence
 * the sequence of random numbers generated. length is the number of entries
 * in the seeds[] array.
 *
 * double Random();
 * Gives a floating point random number in the interval 0 <= x < 1.
 * The resolution is 32 bits in CRandomMother and CRandomMersenne.
 *
 * int IRandom(int min, int max);
 * Gives an integer random number in the interval min <= x <= max.
 * (max-min < MAXINT).
 * The precision is 2^-32 (defined as the difference in frequency between
 * possible output values). The frequencies are exact if max-min+1 is a
 * power of 2.
 *
 * int IRandomX(int min, int max);
 * Same as IRandom, but exact. In CRandomMersenne only.
 * The frequencies of all output values are exactly the same for an
 * infinitely long sequence. (Only relevant for extremely long sequences).
 *
 * uint32 BRandom();
 * Gives 32 random bits.
 *
 *
 * Example:
 * ========
 * #include "mersenne.h"
 *
 * int32 seed = (int32)time(0);        // random seed
 *
 *
 * CRandomMersenne rg(seed);           // make instance of random number generat
 *
 *
 *
 * make random integers in interval from 0 to 99, inclusive:
 *     ir = rg.IRandom(0,99);
 *
 * make random floating point numbers in interval from 0 to 1:
 *     r = rg.Random();
 *
 *  make random bits (hexadecimal):\n");
 *     ir = rg.BRandom();
 *
 *
 *
 *
 *
 * Optimized version:
 * ==================
 * Faster versions of these random number generators are provided as function
 * libraries in asmlib.zip. These function libraries are coded in assembly
 * language and support only x86 platforms, including 32-bit and 64-bit
 * Windows, Linux, BSD, Mac OS-X (Intel based). Use asmlibran.h from asmlib.zip
 *
 *
 * Non-uniform random number generators:
 * =====================================
 * Random number generators with various non-uniform distributions are available
 * in stocc.zip (www.agner.org/random).
 *
 *
 * Further documentation:
 * ======================
 * The file randomc.htm contains further documentation on these random number
 * generators.
 *
 *
 * Copyright:
============
* � 1997 - 2007 Agner Fog. All software in this library is published under the
* GNU General Public License with the further restriction that it cannot be
* used for gambling applications. See licence.htm
*******************************************************************************/


#include <stdio.h>                     // define printf() function
#include <stdlib.h>                    // define exit() function


//#ifndef RANDOMC_H
#define RANDOMC_H


// Define 32 bit signed and unsigned integers.
// Change these definitions, if necessary, to match a particular platform
#if defined(_WIN16) || defined(__MSDOS__) || defined(_MSDOS)
// 16 bit systems use long int for 32 bit integer
typedef long int           int32;   // 32 bit signed integer
typedef unsigned long int  uint32;  // 32 bit unsigned integer
#else
// Most other systems use int for 32 bit integer
typedef int                int32;   // 32 bit signed integer
typedef unsigned int       uint32;  // 32 bit unsigned integer
#endif

// Define 64 bit signed and unsigned integers, if possible
#if (defined(__WINDOWS__) || defined(_WIN32)) && (defined(_MSC_VER) || defined(__INTEL_COMPILER))
// Microsoft and other compilers under Windows use __int64
typedef __int64            int64;   // 64 bit signed integer
typedef unsigned __int64   uint64;  // 64 bit unsigned integer
#define INT64_DEFINED               // Remember that int64 is defined
#elif defined(__unix__) && (defined(_M_IX86) || defined(_M_X64))
// Gnu and other compilers under Linux etc. use long long
typedef long long          int64;   // 64 bit signed integer
typedef unsigned long long uint64;  // 64 bit unsigned integer
#define INT64_DEFINED               // Remember that int64 is defined
#else
// 64 bit integers not defined
// You may include definitions for other platforms here
#endif


/***********************************************************************
Define random number generator classes
***********************************************************************/

class CRandomMersenne {                // Encapsulate random number generator
#if 0
  // Define constants for type MT11213A:
#define MERS_N   351
#define MERS_M   175
#define MERS_R   19
#define MERS_U   11
#define MERS_S   7
#define MERS_T   15
#define MERS_L   17
#define MERS_A   0xE4BD75F5
#define MERS_B   0x655E5280
#define MERS_C   0xFFD58000
#else
  // or constants for type MT19937:
#define MERS_N   624
#define MERS_M   397
#define MERS_R   31
#define MERS_U   11
#define MERS_S   7
#define MERS_T   15
#define MERS_L   18
#define MERS_A   0x9908B0DF
#define MERS_B   0x9D2C5680
#define MERS_C   0xEFC60000
#endif
 public:

  CRandomMersenne(uint32 seed) {      // Constructor
    RandomInit(seed); LastInterval = 0;}

  CRandomMersenne() {      // Constructor
     LastInterval = 0;}



  void RandomInit(uint32 seed) {
    // Initialize and seed
    Init0(seed);

    // Randomize some more
    for (int i = 0; i < 37; i++) BRandom();
  }


  // Seed by more than 32 bits

  void RandomInitByArray(uint32 seeds[], int length) {
    // Seed by more than 32 bits
    int i, j, k;

    // Initialize
    Init0(19650218);

    if (length <= 0) return;

    // Randomize mt[] using whole seeds[] array
    i = 1;  j = 0;
    k = (MERS_N > length ? MERS_N : length);
    for (; k; k--) {
      mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL)) + seeds[j] + j;
      i++; j++;
      if (i >= MERS_N) {mt[0] = mt[MERS_N-1]; i=1;}
      if (j >= length) j=0;}
    for (k = MERS_N-1; k; k--) {
      mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL)) - i;
      if (++i >= MERS_N) {mt[0] = mt[MERS_N-1]; i=1;}}
    mt[0] = 0x80000000UL;  // MSB is 1; assuring non-zero initial array

    // Randomize some more
    mti = 0;
    for (int i = 0; i <= MERS_N; i++) BRandom();
  }



  //int IRandom (int min, int max);     // Output random integer
  //int IRandomX(int min, int max);     // Output random integer, exact
  //double Random();                    // Output random float

  int IRandom(int min, int max) {
    // Output random integer in the interval min <= x <= max
    // Relative error on frequencies < 2^-32
    if (max <= min) {
      if (max == min) return min; else return 0x80000000;
    }
    // Multiply interval with random and truncate
    int r = int((max - min + 1) * Random()) + min;
    if (r > max) r = max;
    return r;
  }


  int IRandomX(int min, int max) {
    // Output random integer in the interval min <= x <= max
    // Each output value has exactly the same probability.
    // This is obtained by rejecting certain bit values so that the number
    // of possible bit values is divisible by the interval length
    if (max <= min) {
      if (max == min) return min; else return 0x80000000;
    }
#ifdef  INT64_DEFINED
    // 64 bit integers available. Use multiply and shift method
    uint32 interval;                    // Length of interval
    uint64 longran;                     // Random bits * interval
    uint32 iran;                        // Longran / 2^32
    uint32 remainder;                   // Longran % 2^32

    interval = uint32(max - min + 1);
    if (interval != LastInterval) {
      // Interval length has changed. Must calculate rejection limit
      // Reject when remainder = 2^32 / interval * interval
      // RLimit will be 0 if interval is a power of 2. No rejection then
      RLimit = uint32(((uint64)1 << 32) / interval) * interval - 1;
      LastInterval = interval;
    }
    do { // Rejection loop
      longran  = (uint64)BRandom() * interval;
      iran = (uint32)(longran >> 32);
      remainder = (uint32)longran;
    } while (remainder > RLimit);
    // Convert back to signed and return result
    return (int32)iran + min;

#else
    // 64 bit integers not available. Use modulo method
    uint32 interval;                    // Length of interval
    uint32 bran;                        // Random bits
    uint32 iran;                        // bran / interval
    uint32 remainder;                   // bran % interval

    interval = uint32(max - min + 1);
    if (interval != LastInterval) {
      // Interval length has changed. Must calculate rejection limit
      // Reject when iran = 2^32 / interval
      // We can't make 2^32 so we use 2^32-1 and correct afterwards
      RLimit = (uint32)0xFFFFFFFF / interval;
      if ((uint32)0xFFFFFFFF % interval == interval - 1) RLimit++;
    }
    do { // Rejection loop
      bran = BRandom();
      iran = bran / interval;
      remainder = bran % interval;
    } while (iran >= RLimit);
    // Convert back to signed and return result
    return (int32)remainder + min;

#endif
  }




  double Random() {
    // Output random float number in the interval 0 <= x < 1
    union {double f; uint32 i[2];} convert;
    uint32 r = BRandom();               // Get 32 random bits
    // The fastest way to convert random bits to floating point is as follows:
    // Set the binary exponent of a floating point number to 1+bias and set
    // the mantissa to random bits. This will give a random number in the
    // interval [1,2). Then subtract 1.0 to get a random number in the interval
    // [0,1). This procedure requires that we know how floating point numbers
    // are stored. The storing method is tested in function RandomInit and saved
    // in the variable Architecture.

    // This shortcut allows the compiler to optimize away the following switch
    // statement for the most common architectures:
#if defined(_M_IX86) || defined(_M_X64) || defined(__LITTLE_ENDIAN__)
    Architecture = LITTLE_ENDIAN1;
#elif defined(__BIG_ENDIAN__)
    Architecture = BIG_ENDIAN1;
#endif

    switch (Architecture) {
    case LITTLE_ENDIAN1:
      convert.i[0] =  r << 20;
      convert.i[1] = (r >> 12) | 0x3FF00000;
      return convert.f - 1.0;
    case BIG_ENDIAN1:
      convert.i[1] =  r << 20;
      convert.i[0] = (r >> 12) | 0x3FF00000;
      return convert.f - 1.0;
    case NONIEEE: default: ;
    }
    // This somewhat slower method works for all architectures, including
    // non-IEEE floating point representation:
    return (double)r * (1./((double)(uint32)(-1L)+1.));
  }


  // Output random bits
  uint32 BRandom() {
    // Generate 32 random bits
    uint32 y;

    if (mti >= MERS_N) {
      // Generate MERS_N words at one time
      const uint32 LOWER_MASK = (1LU << MERS_R) - 1;       // Lower MERS_R bits
      const uint32 UPPER_MASK = 0xFFFFFFFF << MERS_R;      // Upper (32 - MERS_R) bits
      static const uint32 mag01[2] = {0, MERS_A};

      int kk;
      for (kk=0; kk < MERS_N-MERS_M; kk++) {
        y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
        mt[kk] = mt[kk+MERS_M] ^ (y >> 1) ^ mag01[y & 1];}

      for (; kk < MERS_N-1; kk++) {
        y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
        mt[kk] = mt[kk+(MERS_M-MERS_N)] ^ (y >> 1) ^ mag01[y & 1];}

      y = (mt[MERS_N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
      mt[MERS_N-1] = mt[MERS_M-1] ^ (y >> 1) ^ mag01[y & 1];
      mti = 0;
    }

    y = mt[mti++];

#if 1
    // Tempering (May be omitted):
    y ^=  y >> MERS_U;
    y ^= (y << MERS_S) & MERS_B;
    y ^= (y << MERS_T) & MERS_C;
    y ^=  y >> MERS_L;
#endif

    return y;
  }



 private:

  void Init0(uint32 seed) {
    // Detect computer architecture
    union {double f; uint32 i[2];} convert;
    convert.f = 1.0;
    if (convert.i[1] == 0x3FF00000) Architecture = LITTLE_ENDIAN1;
    else if (convert.i[0] == 0x3FF00000) Architecture = BIG_ENDIAN1;
    else Architecture = NONIEEE;

    // Seed generator
    mt[0]= seed;
    for (mti=1; mti < MERS_N; mti++) {
      mt[mti] = (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
    }
  }



  uint32 mt[MERS_N];                  // State vector
  int mti;                            // Index into mt
  uint32 LastInterval;                // Last interval length for IRandomX
  uint32 RLimit;                      // Rejection limit used by IRandomX
  enum TArch {LITTLE_ENDIAN1, BIG_ENDIAN1, NONIEEE}; // Definition of architecture
  TArch Architecture;                 // Conversion to float depends on architecture
};


//#endif
