"""Validate the exact specification and generate rational verifier input."""
from fractions import Fraction
from pathlib import Path
import json
EXPECTED = {'degree': 29984, 'm': 16, 'N': 1874, 'L': 1872, 'correction_blocks': 64, 'lambda': '57607457624961133635359/50000000000000000000', 'correction_scale': '6/5', 'inner_radius': '49/50', 'outer_radius': '1999/2000', 'target_margin': '1/10000000000', 'annular_margin': '1/5000000', 'hermite': True, 'tail_blocks': 8, 'extra_correction_blocks': 128}

BASE_POINT = {'re': '1.04015231339331770677e+00', 'im': '3.02102329270430736141e-01'}

def parameter_lines(spec):
    if spec.get('primitive',False):raise ValueError('This certificate uses q^N, not a primitive.')
    for key in ['degree','m','N','L','correction_blocks']:
        if spec.get(key)!=EXPECTED[key]:
            raise ValueError(f"This certificate requires {key}={EXPECTED[key]}")
    for key in ['lambda','correction_scale','inner_radius','outer_radius','target_margin','annular_margin']:
        if Fraction(spec[key])!=Fraction(EXPECTED[key]):
            raise ValueError(f"Unexpected {key}")
    if spec['S_definition']!='A+1/A':raise ValueError('Unexpected S definition')
    if ('hermite_coefficient' in spec)!=EXPECTED['hermite']:
        raise ValueError('Hermite coefficient presence does not match the verifier')
    b=spec['b']
    out=[spec['A'],spec['r']['re'],spec['r']['im'],spec['a']['re'],spec['a']['im'],
         b['re_numerator']+'/'+b['denominator'],b['im_numerator']+'/'+b['denominator']]
    if len(spec['correction_centers'])!=2:raise ValueError('Expected two centers')
    for x in spec['correction_centers']:out += [x['re'],x['im']]
    rows=sorted(spec['correction_coefficients'],key=lambda x:(x['ell'],x['k']))
    if [(x['ell'],x['k']) for x in rows]!=[(j,k) for j in range(EXPECTED['correction_blocks']) for k in range(16)]:
        raise ValueError('Incomplete or duplicated correction coefficient table')
    for x in rows:out += [x['re'],x['im']]
    if EXPECTED['hermite']:
        x=spec['hermite_coefficient'];out += [x['re'],x['im']]
    if Fraction(spec['lambda']) != EXPECTED['N']*Fraction(spec['high_power_coefficient']):raise ValueError('Amplitude mismatch')
    out += [spec['high_power_coefficient']]
    if len(spec.get('free_coefficients', [])) != 4: raise ValueError('Expected a complete cubic remainder')
    if spec.get('free_center_definition') != '(correction_centers[0]+correction_centers[1])/2, exactly': raise ValueError('Unexpected free center')
    for x in spec['free_coefficients']: out += [x['re'], x['im']]
    if spec.get('tail_blocks')!=EXPECTED['tail_blocks']:raise ValueError('Wrong top-band width')
    if spec.get('tail_base_exponent')!=EXPECTED['N']-EXPECTED['tail_blocks']:raise ValueError('Wrong top-band exponent')
    tr=sorted(spec.get('tail_coefficients',[]),key=lambda x:(x['ell'],x['k']))
    if [(x['ell'],x['k']) for x in tr]!=[(j,k) for j in range(EXPECTED['tail_blocks']) for k in range(16)]:raise ValueError('Incomplete top band')
    for x in tr:out += [x['re'],x['im']]
    if spec.get('extra_correction_blocks')!=EXPECTED['extra_correction_blocks']:raise ValueError('Wrong extra-correction block count')
    if spec.get('extra_correction_definition')!='(z-a)^2*sum_{ell=0}^{extra_correction_blocks-1}sum_{k=0}^{15}x[ell,k]*(z/(6/5))^k*q(z)^ell':raise ValueError('Wrong extra-correction definition')
    er=sorted(spec.get('extra_correction_coefficients',[]),key=lambda x:(x['ell'],x['k']))
    if [(x['ell'],x['k']) for x in er]!=[(j,k) for j in range(EXPECTED['extra_correction_blocks']) for k in range(16)]:raise ValueError('Incomplete extra-correction table')
    for x in er:out += [x['re'],x['im']]
    phase = spec.get('tail_phase', {})
    if Fraction(phase['re']) != Fraction('-77141935446226489841583961398938702277129/77191935446226489841583961398938702277129') or Fraction(phase['im']) != Fraction('2777893724501110028770000000000000000000/77191935446226489841583961398938702277129'):
        raise ValueError('Unexpected tail phase')
    if Fraction(phase['re'])**2+Fraction(phase['im'])**2 != 1:
        raise ValueError('Tail phase must be exactly unimodular')
    out += [phase['re'], phase['im']]
    nb=spec['normalization_base']
    expected_base=BASE_POINT
    if any(Fraction(nb[k])!=Fraction(expected_base[k]) for k in ['re','im']):
        raise ValueError('Normalization base mismatch')
    if spec['normalization_definition']!='p(x)=(P(normalization_base+x)-P(normalization_base))/Pprime(normalization_base), exactly':
        raise ValueError('Wrong normalization definition')
    out += [nb['re'],nb['im']]
    return [str(Fraction(x)) for x in out]

if __name__=='__main__':
    root=Path(__file__).resolve().parent
    exact='\n'.join(parameter_lines(json.loads((root/'polynomial.json').read_text())))+'\n'
    if (root/'parameters.txt').read_text()!=exact:raise SystemExit('FAIL: exact input mismatch')
    print('PASS: polynomial.json and parameters.txt agree as exact rational numbers')
