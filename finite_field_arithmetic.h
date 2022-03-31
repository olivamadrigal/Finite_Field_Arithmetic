#ifndef finite_field_arithmetic_h
#define finite_field_arithmetic_h

/*

In GF(2) we have only two elements: 0 and 1.
Addition and subtration is the equivalent of bitwise XOR
Multiplication is a bitwise AND

Polynomials with coefficients defined over GF(2) can be
represented by binary strings with powers matching the
bit position

x^3 + x^0 = 1 0 0 1 <--- binary string
            3 2 1 0 <--- bit position

Polynomial Arithmetic with coefficients defined over GF(2):

Let GF(p^n) define a Finite Field with p prime defined by
an irreducible polynomial F(t) of degree n.

For security, when F(t) is the irreducible polynomial, we only use
odd prime curves, meaning that n = prime odd.

F(t) has degree(n) represented by n+1 bits
the elemts in F(t) have degree n-1 and are represented by n bits

For elements  and b in the field:

a + b = a - b = XOR

a x b = modular multiplication reducing with m(x)
Step1: c = a x b (multiply out as usual) can implement as shift and add
Step2: if the resulting degree of c is => n, reduce the product using m(x)
via Montgomery, Blakely, or BMM (with threads to exploit parallelism)

Can implement Step1 and Step2 as a single interleaved algorithm.

a / b = a x MI(b) where MI(b) the multiplicative inverse of b mod m(x).

We find MI using Bezout identity and EEA.

For C: max bit length variable may be uint64_t or uint128_t depending on your
platform and compiler.

To work with larger fields, you must use multiprecision and compute accordingly.
For example:

uint64_t E[2]; //can be used to represent a uint128_t as: [63:0][63:0]

Because C does not provide exact bit lengths... you must ensure n fits.
and that the max product of any two elments will fit.

If n = 9, => 10bits to represent
One element has degree d = 8 => 9 bits to represent.
Product of any two elements has max degree Cd = 16 => 17 bits to represent..

For this example: we consider the finite field GF(2^31)
with m(x) = x^31 + x + 1

n = 31 => 32-bits
d = 30 => 31-bits
Cd max = 60 => 61 bits.

we are good with a uint64_t.

https://www.wolframalpha.com/input/?i=IrreduciblePolynomialQ%5Bx%5E31+%2B+x+%2B+1%5Dmod2

this library uses a tiny field of GF(2^3) so that you can see everything and can easily
verify by hand.

it can easily extend to larger fields within small degree and to large fields via multiprecision
in such cases, we must be careful with the word or register size for some computations, such
as when we shift a value.. to have space for the max possible shifts.

static const uint8_t n = 31;
static const uint32_t f = 0x80000003;   //irreducible polynomial

*/

#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

static const uint8_t n = 3; //degree of f
static const uint8_t f = 0x0b;//irreducible polynomial x^3 + x + 1 or 1011 :)

/*-----------------------------------------------------------------------
                            Cardinality
                returns number of elements in the field
 -------------------------------------------------------------------------*/
uint8_t cardinality(void);

/*-----------------------------------------------------------------------
                                Cd Max
        position of the most significant bit set (degree of c)
 -------------------------------------------------------------------------*/
uint8_t degree(uint8_t *c);

/*------------------------------------------------------------------------
                 Generate Elements in the Field
 -------------------------------------------------------------------------*/
uint8_t *GenerateElements(void);

/*------------------------------------------------------------------------
                        GF Addition = GF Subtraction
 -------------------------------------------------------------------------*/
uint8_t GF2nAdd(uint8_t *a, uint8_t *b);

/*------------------------------------------------------------------------
                    GF Multiplication: c = a x b
    normal straighout multiplication example: (x^2 + x)x => x^3 + x^2.
    using shift and add method
 -------------------------------------------------------------------------*/
 uint8_t GF2nMul(uint8_t *a, uint8_t *b);
/*------------------------------------------------------------------------
                    Blakely Reduction: residue = c mod f
                where c = a x b  of c degree >= n.
 -------------------------------------------------------------------------*/
uint8_t GF2nBlakelyReduction(uint8_t *c);

/*------------------------------------------------------------------------
        Interleaved Modular Multiplication and Reduction in GF(2^n)
                        Using Blakely Approach
                        returns a*b MOD F(t)
 -------------------------------------------------------------------------*/
uint8_t InterleavedBlakely(uint8_t *a, uint8_t *b);

/*------------------------------------------------------------------------
                    computes c(t) MOD f(t)
        return quotient q and remainder r through pointers
 -------------------------------------------------------------------------*/
void GF2nDivide(uint8_t *c, uint8_t *F, uint8_t *q, uint8_t *r);

/*------------------------------------------------------------------------
                    Multiplicative Inverse mod F
 
Use of this depends. If we are working with a relatively small curve like
163. We may just precompute all the MI... similar to Sbox and SboxInverse
in AES. Different algorithms depnding on the application.

a = element in GF(2^n) and a1 is it's multiplicative inverse

Example: MI(x^2 + 1) =>  EEA(a,b) = EEA(x^2 + 1, F(t)):
x^2 + 1 ) F(t) => Q = x and R = 1 using normal polynomial division
and keeping +/- as XOR (GF(2) arithmetic)
Verify:  a*a^-1 = 1 mod F(t)
x^2 + 1 (x) = x^3 + x ==> F(t) ) x^3 + x ==> R=1,Q=1 or simple reduction:

1010 a*a^-1
1011 F(t)
----
0010 = x
 -------------------------------------------------------------------------*/
uint8_t MI(uint8_t *a);

/*------------------------------------------------------------------------
                a(x)/b(x) MOD F(t) => a(x)b(x)^-1 MOD F(t)
 -------------------------------------------------------------------------*/
void GF2nDivision(uint8_t *a, uint8_t *b, uint8_t *q, uint8_t *r);

#endif /* finite_field_arithmetic_h */
