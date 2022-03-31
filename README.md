# Finite_Field_Arithmetic
Finite Field Arithmetic

Finite Field Arithmetic in GF(2^n), uses GF(2^3) as example.

Elements are polynomials of degree at most 2, e.g. x^2 + x + 1 with coefficients defined in GF(2) which allow for direct mapping
in binary. For example: x^2 + x + 1 = 111 (bit positions map to the degree). Elements are printed in hexdecimal representation.


- Addition / Subtraction: a + b = a - b = a ^ b

- Multiplication: a x b (straight out polynomial multiplication

- Modular Reduction: given a product c = a x b, compute residue e = c MOD F

- Blakely or Standard Interleaved Modular Multiplication: multiply and reduce a product

- Divide: a x b = c MOD F, returns quotient and remainder

- Multiplicative Inverse: MI(a) s.t. a*MI(a) = 1 MOD F

- Division: a / b = a x MI(b) MOD F
