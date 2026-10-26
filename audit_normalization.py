#!/usr/bin/env python3
"""Audit exact normalization data and recorded nonzero derivative enclosures.
This checks completed evidence; it does not replace the interval computation.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,re
root=Path(__file__).resolve().parent
spec=json.loads((root/'polynomial.json').read_text())
nb=spec['normalization_base'];rows=[]
for k in range(16):
 text=(root/'logs'/f'part_{k:02d}.log').read_text()
 lines=[s for s in text.splitlines() if s.startswith('NORMALIZATION ')]
 if len(lines)!=1:raise SystemExit('Missing or ambiguous normalization record')
 pairs=re.findall(r'\[([^,\]]+),([^\]]+)\]',lines[0])
 if len(pairs)!=6:raise SystemExit('Unexpected normalization format')
 iv=[]
 for lo,hi in pairs:
  a,b=F(lo),F(hi);guard=(1+max(abs(a),abs(b)))*F(1,10**18);iv.append((a-guard,b+guard))
 for j,key in enumerate(['re','im']):
  if not iv[j][0]<=F(nb[key])<=iv[j][1]:raise SystemExit('Wrong normalization point')
 if not any(a>0 or b<0 for a,b in iv[4:]):raise SystemExit('Derivative rectangle includes zero')
 if not re.search(r'^PASS mode=allroots sheet='+str(k)+r' .*degree='+str(spec['degree'])+' ',text,re.M):raise SystemExit('Incomplete partition')
 rows.append({'partition':k,'record':lines[0]})
assert spec['normalization_definition']=='p(x)=(P(normalization_base+x)-P(normalization_base))/Pprime(normalization_base), exactly'
assert spec['normalization_base']!=spec['a']
phase=spec['tail_phase'];assert F(phase['re'])**2+F(phase['im'])**2==1
out={'status':'PASS','degree':spec['degree'],'normalization_base':nb,'distinct_from_model_anchor':True,'nonzero_normalization_derivative_in_every_partition':True,'phase_exactly_unimodular':True,'records':rows,'polynomial_sha256':hashlib.sha256((root/'polynomial.json').read_bytes()).hexdigest(),'warning':'Audit of recorded evidence, not an independent recomputation of analytic enclosures.'}
(root/'normalization_audit.json').write_text(json.dumps(out,indent=2)+'\n');print('PASS normalization and nonzero derivative audit')
