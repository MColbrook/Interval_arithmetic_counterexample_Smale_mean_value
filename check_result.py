#!/usr/bin/env python3
"""Check the theorem against the complete stored interval certificate."""
from pathlib import Path
from decimal import Decimal
from fractions import Fraction
import hashlib, json, re
ROOT=Path(__file__).resolve().parent

def main():
    from audit_certificate import audit
    a=audit(ROOT)
    result=json.loads((ROOT/'result.json').read_text())
    exact=json.loads((ROOT/'independent_exact_audit.json').read_text())
    n=result['degree']
    if result['status']!='PASS' or n!=29984:
        raise ValueError('Incorrect theorem summary')
    if a['degree']!=n or a['critical_points_counted_with_multiplicity']!=n-1 or not a['all_roots_simple']:
        raise ValueError('Critical-set statement is not established')
    if exact['status']!='PASS' or exact['strict_rouche_and_quotient_records_checked']!=n-1:
        raise ValueError('Exact rational check is incomplete')
    h=hashlib.sha256((ROOT/'polynomial.json').read_bytes()).hexdigest()
    if result['specification_sha256']!=h or exact['input_sha256']['polynomial.json']!=h:
        raise ValueError('Certificate and theorem refer to different polynomials')
    lo,hi=Decimal(result['simple_safe_lower']),Decimal(result['simple_safe_upper'])
    if not Decimal(1)<lo<Decimal(a['minimum_safe_lower']) or not hi>Decimal(a['minimum_safe_upper']):
        raise ValueError('Current enclosures do not imply the theorem bounds')
    normalization=json.loads((ROOT/'normalization_audit.json').read_text())
    if normalization['status']!='PASS' or normalization['polynomial_sha256']!=h or normalization['degree']!=n:
        raise ValueError('Normalization check is incomplete')
    if sorted(x['partition'] for x in normalization['records'])!=list(range(16)):
        raise ValueError('Normalization check does not cover all partitions')
    if not result['unique_minimizer_proved'] or not a['unique_minimizer_proved']:
        raise ValueError('Unique minimum is not proved')
    if result['minimizing_disk']!=a['minimizing_disk']:
        raise ValueError('Minimizing disk does not match the certificate')
    if exact.get('refinement_log_sha256')!=hashlib.sha256((ROOT/'certificate_refined.log').read_bytes()).hexdigest():
        raise ValueError('Exact check refers to a different refinement')
    text=(ROOT/'base_comparison.log').read_text()
    pattern=r'^OLD_BASE_FAILURE root=27158 quotient_lower=\[([^,]+),([^\]]+)\] quotient_upper=\[([^,]+),([^\]]+)\] root_count=1 comparison=STRICT_PASS$'
    rows=re.findall(pattern,text,re.M)
    if len(rows)!=1 or not re.search(r'^PASS BASE_COMPARISON degree=29984$',text,re.M):
        raise ValueError('Model-anchor comparison is incomplete')
    low=Fraction(rows[0][0]);high=Fraction(rows[0][3]);guard=Fraction(1,10**20)
    if not low-abs(low)*guard>Fraction('0.99982426177940') or not high+abs(high)*guard<Fraction('0.99982426177944'):
        raise ValueError('Model-anchor comparison does not imply the stated bounds')
    print('PASS: degree 29984; all 29983 critical points simple; published bounds and unique minimum established.')
    print('PASS: the stored local comparison establishes a quotient below one at the model anchor.')

if __name__=='__main__':
    main()
