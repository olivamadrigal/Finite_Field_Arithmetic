#include "finite_field_arithmetic.h"

int main(void)
{
    uint8_t *elements, count, i, j, c, q, r;
    
    count = cardinality();
    printf("\nThere are %d elements in GF(2^3).\n", count);
    printf("\nElements in GF(2^3) in hex format:\n");
    elements = GenerateElements();
    for(uint8_t i =0; i <= 7; i++)
        printf("%02x ", i);
    printf("\n");
    printf("\nThe corresponding degree of each element in GF(2^3):\n");
    for(uint8_t i =0; i <= 7; i++)
        printf("%d ", degree(&i));
    printf("\n");
    
    printf("\n\nGF(2^n) addition: \n");
    count /= 2;
    for(i = 0, j = count; i < count; i++, j++)
    {
        c = GF2nAdd(&elements[i], &elements[j]);
        printf("%02x + %02x = %02x\n", elements[i], elements[j], c);
    }
    
    printf("\n\nGF(2^n) multiplication (a x b) = c: \n");
    for(i = 0, j = count; i < count; i++, j++)
    {
        c = GF2nMul(&elements[i], &elements[j]);
        printf("%02x x %02x = %02x\n", elements[i], elements[j], c);
        r = GF2nBlakelyReduction(&c);
        printf("Blakely Residue: %02x \n", r);//very easy to verify by hand.. c xor f (from the left)
        
    }
    printf("\n\nGF(2^n) modular multiplication (a x b) MOD F(t) = r : \n");
    for(i = 0, j = count; i < count; i++, j++)
    {
        c = InterleavedBlakely(&elements[i], &elements[j]);
        printf("%02x + %02x = %02x\n", elements[i], elements[j], c);
    }
    
    printf("\n\nGF(2^n) modular multiplication (a x b) MOD F(t) = r : \n");
    for(i = 0, j = count; i < count; i++, j++)
    {
        c = InterleavedBlakely(&elements[i], &elements[j]);
        printf("%02x + %02x = %02x\n", elements[i], elements[j], c);
    }
        
    printf("\n\nGF(2^n) division quotient, residue of c MOD F(t) where c = a x b. \n");
    for(i = 0, j = count; i < count; i++, j++)
    {
        c = GF2nMul(&elements[i], &elements[j]);
        GF2nDivide(&c, &f, &q, &r);
        printf("%02x has q= %02x r = %02x\n", c, q, r);
    }
    
    c = 0xa1;
    GF2nDivide(&c, &f, &q, &r);
    printf("~%02x has q= %02x r = %02x\n", c, q, r);
    
    printf("\nFor each element, compute the MI:\n");
    for(uint8_t i = 2; i <= 7; i++)
    {
        printf("a = %02x a^-1 = %02x \n", i, MI(&i));
    }
      
   return 0;
}
