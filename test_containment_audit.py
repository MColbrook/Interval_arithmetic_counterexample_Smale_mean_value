#!/usr/bin/env python3
"""Regression tests for the exact unshifted-radius case; not a root certificate."""
from fractions import Fraction as F
from decimal import Decimal,localcontext
from audit_exact_disks import exact_unshifted_equality

def main():
    count=0
    with localcontext() as c:
        c.prec=120
        for k in range(20,31):
            r=F(1,2**k);s=str(Decimal(r.numerator)/Decimal(r.denominator))
            assert exact_unshifted_equality(s,s,r,64);count+=1
            assert not exact_unshifted_equality(s,s,r,113);count+=1
            assert not exact_unshifted_equality(s,s,r*3,64);count+=1
            altered=str(Decimal(r.numerator)/Decimal(r.denominator)+Decimal('1e-18'))
            assert not exact_unshifted_equality(s,altered,r,64);count+=1
    assert not exact_unshifted_equality('0','0',F(0),64);count+=1
    assert not exact_unshifted_equality('-1','-1',F(-1),64);count+=1
    print(f'PASS: {count} exact-radius regression cases; ambiguous or altered cases rejected.')
if __name__=='__main__':main()
