#include "finite_field_arithmetic.h"

/*-----------------------------------------------------------------------
                                Cardinality
                returns number of elements in the field
 -------------------------------------------------------------------------*/
uint8_t cardinality(void)
{
    return (uint8_t)pow(2,n);
}

/*------------------------------------------------------------------------
                                Cd Max
        position of the most significant bit set (degree of c)
 -------------------------------------------------------------------------*/
uint8_t degree(uint8_t *c)
{
    return log2(*c);
}

/*------------------------------------------------------------------------
                 Generate Elements in the Field
 -------------------------------------------------------------------------*/
uint8_t *GenerateElements(void)
{
    uint8_t *elements;
    int8_t d;

    d = 7; //max value represented by n bits = 2^n - 1
    elements = calloc(d, sizeof(uint8_t));
    for(int i = d; i >= 0; i--)
        elements[i] = i;
    return elements;
}

/*------------------------------------------------------------------------
                        GF Addition = GF Subtraction
 -------------------------------------------------------------------------*/
uint8_t GF2nAdd(uint8_t *a, uint8_t *b)
{
    return *a ^ *b;
}

/*------------------------------------------------------------------------
                    GF Multiplication: c = a x b
    normal straighout multiplication example: (x^2 + x)x => x^3 + x^2.
    using shift and add method
 -------------------------------------------------------------------------*/
 uint8_t GF2nMul(uint8_t *a, uint8_t *b)
{
    uint8_t c, mask;
    
    mask = 0x01;
    for(int32_t i = 2; i >= 0; i--)
        c = (*a & (mask << i)) ? (c ^ (*b << i)): c;
    return c;
}

/*------------------------------------------------------------------------
                    Blakely Reduction: residue = c mod f
                where c = a x b  of c degree >= n.
 -------------------------------------------------------------------------*/
uint8_t GF2nBlakelyReduction(uint8_t *c)
{
  uint8_t R, mask, d;
  
  R = *c;
  d = degree(c);
  mask = 0x01;
  for(uint8_t i = d; i <= n; i = i + 1)
    R = (R & (mask << i))? R ^ (f << (i-n)) : R;
    
  return R;
}

/*------------------------------------------------------------------------
        Interleaved Modular Multiplication and Reduction in GF(2^n)
        Using Blakely Approach; returns a*b MOD F(t)
 -------------------------------------------------------------------------*/
uint8_t InterleavedBlakely(uint8_t *a, uint8_t *b)
{
    uint8_t r, msb;

    r = 0x00;
    msb = 0x01 << n;
    for(int8_t i = 2; i >= 0; i--)
    {
        r = r << 1;
        r = (*a & (0x01 << i))? r ^ *b : r;
        r = (r & msb)? r ^ f : r;
    }
    return r;
}

/*------------------------------------------------------------------------
                    computes c(t) MOD f(t)
        return quotient q and remainder r through pointers
 -------------------------------------------------------------------------*/
void GF2nDivide(uint8_t *c, uint8_t *F, uint8_t *q, uint8_t *r)
{
    uint8_t deg_c, deg_f, ft;
    
    deg_c = degree(c);
    deg_f = degree(F);
    if(deg_f == 0)
    {
            *q = *c;
            *r = 0x00;
            return;
    }
    *r = *c;
    *q = 0x00;
    ft = *F;
    for(uint8_t i = deg_c, j=deg_c - deg_f; i >= deg_f; i--, j--)
    {
      
       if(*r & (0x01 << i))
       {
            *q = *q | (0x01 << j);
            *r = *r ^ (ft << (i-deg_f));
        }
    }
}

/*------------------------------------------------------------------------
                    Multiplicative Inverse mod F
 
Use of this depends. If we are working with a relatively small curve like
163. We may just precompute all the MI... similar to Sbox and SboxInverse
in AES. Different algorithms depnding on the application...

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

Note: because this is Galois and GF(2^p) is a prime finite field,
it is the case that every element except the identity element of addition
has a multiplicative inverse.
 -------------------------------------------------------------------------*/
uint8_t MI(uint8_t *a)
{
    uint8_t ap, fp, x, xp, y, yp, q, r, t1, t2;
    
    if(*a == 0)return 0;
    if(*a == 1)return 1;
    fp = f;
    ap = *a;
    x = 0x00;
    xp = 0x01;
    y = 0x01;
    yp = 0x00;
    //solve EEA [gdc(a, F) = 1 = ax + Fy] MOD F where x would be the MI.
    //and since any multplice of F mod F is zero.
    while(fp)
    {
        GF2nDivide(&ap, &fp, &q, &r);
        ap = fp;
        fp = r;
        t1 = x;
        x = xp ^ GF2nMul(&q,&x);
        xp = t1;
        t2 = y;
        y = yp ^ GF2nMul(&q,&y);
        yp = t2;
    }
    if(ap != 1)return ap; //in this case, that would be the gcd
    t1 = xp ^ f;
    GF2nDivide(&t1, &f, &q, &r);
    return r;
}

/*------------------------------------------------------------------------
                a(x)/b(x) MOD F(t) => a(x)b(x)^-1 MOD F(t)
 -------------------------------------------------------------------------*/
void GF2nDivision(uint8_t *a, uint8_t *b, uint8_t *q, uint8_t *r)
{
    uint8_t mi, c;
    
    mi = MI(b);
    c = GF2nMul(a, &mi);
    GF2nDivide(&c, &f, q, r );
}

