#!/usr/bin/env python3
"""Audit exact inputs, complete index coverage, bounds, and stored hashes.

This is a consistency/integrity audit of an already completed interval run.
It does not replace compiling and rerunning the mathematical verifier.
"""
from __future__ import annotations
import argparse, hashlib, importlib.util, json, re
from decimal import Decimal, getcontext
from fractions import Fraction
from pathlib import Path
getcontext().prec=90

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def guard(v):return (1+abs(v))*Decimal('1e-18')
def records(text,prefix):
    out={}
    for row in re.finditer(r'^'+prefix+r' (\d+) .* quotient_lower=([^ ]+) quotient_upper=([^\n]+)$',text,re.M):
        j=int(row[1]);lo=Decimal(row[2]);hi=Decimal(row[3]);lo-=guard(lo);hi+=guard(hi)
        if j in out or lo<=1 or hi<lo:raise ValueError('Invalid or repeated interval record')
        out[j]=(lo,hi)
    return out

def audit(root:Path):
    s=json.loads((root/'polynomial.json').read_text()); n=s['degree']
    mi=importlib.util.spec_from_file_location('exact_inputs',root/'make_parameters.py');m=importlib.util.module_from_spec(mi);mi.loader.exec_module(m)
    if '\n'.join(m.parameter_lines(s))+'\n'!=(root/'parameters.txt').read_text():raise ValueError('Exact input mismatch')
    if n!=s['m']*s['N'] or s['m']!=16 or s['L']!=s['N']-2:raise ValueError('Wrong degree identity')
    if Fraction(s['high_power_coefficient'])<=0:raise ValueError('Leading coefficient may vanish')
    b=s.get('extra_correction_blocks',0)
    degree_bounds={'model_term':16*s['L']+16,'original_correction':16*s['correction_blocks']+5,'extra_correction':16*b+1 if b else 0,'nonleading_tail':16*s['N']-1,'leading_term':16*s['N']}
    if any(v>=n for k,v in degree_bounds.items() if k!='leading_term'):raise ValueError('Leading-degree cancellation not excluded')
    summary=json.loads((root/'certificate_summary.json').read_text())
    if summary['status']!='PASS':raise ValueError('No completed base run')
    for name,h in summary['inputs']['sha256'].items():
        if digest(root/name)!=h:raise ValueError('Input changed: '+name)
    if digest(root/'certificate.log')!=summary['certificate_log_sha256']:raise ValueError('Main log hash mismatch')
    rows={}; partitions=[]
    for k in range(16):
        log=root/'logs'/f'part_{k:02d}.log'; text=log.read_text()
        if not re.search(r'^PASS mode=allroots sheet='+str(k)+r' .*degree='+str(n)+r' ',text,re.M):raise ValueError('Incomplete partition')
        if not re.search(r'^DISJOINT count='+str(n-1)+r' ',text,re.M):raise ValueError('No disjointness check')
        part=records(text,'CRITICAL_DISK')
        if sorted(part)!=list(range(k,n-1,16)):raise ValueError('Partition index mismatch')
        rows.update(part);partitions.append({'partition':k,'roots':len(part),'sha256':digest(log)})
    if sorted(rows)!=list(range(n-1)):raise ValueError('Critical-set coverage mismatch')
    base_lower=min(x[0] for x in rows.values());refined=False
    refsummary=root/'certificate_refined_summary.json'
    if refsummary.exists():
        ref=json.loads(refsummary.read_text())
        if ref['status']!='PASS':raise ValueError('Invalid refinement status')
        for name,h in ref['sha256'].items():
            if digest(root/name)!=h:raise ValueError('Refinement input changed: '+name)
        for k in range(16):
            if digest(root/'logs'/f'part_{k:02d}.log')!=ref['base_log_sha256'][str(k)]:raise ValueError('Refinement base-log mismatch')
        log=(root/'certificate_refined.log').read_text(); rr=records(log,'REFINED_DISK')
        wanted=sorted(int(x) for x in (root/'refinement_indices.txt').read_text().split())
        if sorted(rr)!=wanted or not re.search(r'^PASS REFINEMENT degree='+str(n)+r' count='+str(len(rr))+r' ',log,re.M):raise ValueError('Refinement coverage failure')
        for j,(lo,hi) in rr.items():
            old=rows[j];rows[j]=(max(lo,old[0]),min(hi,old[1]))
            if rows[j][0]>rows[j][1]:raise ValueError('Contradictory enclosures')
        refined=True
    lower=min(v[0] for v in rows.values());winner=min(rows,key=lambda j:rows[j][1]);upper=rows[winner][1]
    other=min(v[0] for j,v in rows.items() if j!=winner)
    out={'status':'PASS: exact-input and certificate-coverage audit','degree':n,'critical_points_counted_with_multiplicity':n-1,'partitions':partitions,'all_roots_simple':True,'degree_bounds':degree_bounds,'base_safe_lower':str(base_lower),'refinement_included':refined,'minimum_safe_lower':str(lower),'minimum_safe_upper':str(upper),'minimizing_disk':winner,'unique_minimizer_proved':upper<other,'all_other_quotients_greater_than':str(other),'specification_sha256':digest(root/'polynomial.json'),'arithmetic_source_sha256':digest(root/'rigorous_interval.hpp'),'verifier_source_sha256':digest(root/'verify.cpp'),'warning':'This audit checks exact input consistency, result coverage, and artifact integrity; it does not rerun the interval calculations.'}
    (root/'certificate_audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='partitions'},indent=2))
    return out
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--directory',type=Path,default=Path(__file__).resolve().parent);args=p.parse_args();audit(args.directory.resolve())
