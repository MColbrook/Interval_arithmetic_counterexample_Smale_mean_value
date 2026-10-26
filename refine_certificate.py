#!/usr/bin/env python3
"""Refine only critical disks that could contain the smallest quotient.

This requires a completed exhaustive certificate. It does not replace the full
root count, and it never infers an inequality from a numerical root list.
"""
from pathlib import Path
from fractions import Fraction
from decimal import Decimal,getcontext
import argparse,subprocess,json,re,hashlib,importlib.util,time,shutil
getcontext().prec=80

def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def parse(text,prefix):
 out={}
 for x in re.finditer(r'^'+prefix+r' (\d+) .* quotient_lower=([^ ]+) quotient_upper=([^\n]+)$',text,re.M):
  k=int(x[1]);lo=Decimal(x[2]);hi=Decimal(x[3]);lo-=(abs(lo)+1)*Decimal('1e-18');hi+=(abs(hi)+1)*Decimal('1e-18')
  if k in out:raise RuntimeError('Repeated disk index')
  if not lo>Decimal('1.0000000001') or hi<lo:raise RuntimeError('Inadequate/invalid printed enclosure')
  out[k]=(lo,hi)
 return out

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--directory',type=Path,default=Path(__file__).resolve().parent);ap.add_argument('--compiler',default='g++');args=ap.parse_args();root=args.directory.resolve();start=time.time()
 spec=json.loads((root/'polynomial.json').read_text());degree=spec['degree']
 im=importlib.util.spec_from_file_location('exact_inputs',root/'make_parameters.py');mod=importlib.util.module_from_spec(im);im.loader.exec_module(mod)
 if '\n'.join(mod.parameter_lines(spec))+'\n'!=(root/'parameters.txt').read_text():raise RuntimeError('Exact parameter mismatch')
 # Bind the refinement to the input/source identity of the base certificate.
 if not (root/'certificate_summary.json').exists():
  raise RuntimeError('No completed base certificate')
 identity=json.loads((root/'certificate_summary.json').read_text())['inputs']['sha256']
 for name in ['verify.cpp','rigorous_interval.hpp','parameters.txt','critical_disks.txt']:
  if digest(root/name)!=identity[name]:raise RuntimeError('Base input/source changed: '+name)
 rows={};loghashes={}
 for k in range(16):
  path=root/'logs'/f'part_{k:02d}.log';text=path.read_text()
  if not re.search(r'^PASS mode=allroots sheet='+str(k)+r' .*degree='+str(degree)+r' ',text,re.M):raise RuntimeError('Base partition incomplete')
  part=parse(text,'CRITICAL_DISK')
  if sorted(part)!=list(range(k,degree-1,16)):raise RuntimeError('Base index coverage mismatch')
  rows.update(part);loghashes[str(k)]=digest(path)
 if sorted(rows)!=list(range(degree-1)):raise RuntimeError('Incomplete base root count')
 upper=min(x[1] for x in rows.values())
 indices=sorted(k for k,(lo,hi) in rows.items() if lo<=upper+Decimal('1e-12'))
 (root/'refinement_indices.txt').write_text('\n'.join(map(str,indices))+'\n')
 source=root/'refine_disks.cpp'
 if not source.exists():shutil.copy2(Path(__file__).resolve().parent/'refine_disks.cpp',source)
 flags=['-O2','-std=c++17','-fno-fast-math','-ffp-contract=off','-frounding-math',f'-DN_PARAM={spec["N"]}',f'-DBLOCKS_PARAM={spec["correction_blocks"]}',f'-DTAIL_BLOCKS_PARAM={spec["tail_blocks"]}',f'-DEXTRA_BLOCKS_PARAM={spec["extra_correction_blocks"]}']
 variants=[['-lmpfr','-lgmpxx','-lgmp'],['-l:libmpfr.so.6','-lgmpxx','-lgmp']]
 for links in variants:
  result=subprocess.run([args.compiler,*flags,'refine_disks.cpp','-o','refine_disks',*links],cwd=root,capture_output=True,text=True)
  if result.returncode==0:break
 else:raise RuntimeError(result.stdout+result.stderr)
 print('REFINING',len(indices),'of',degree-1,'counted disks',flush=True)
 log=root/'certificate_refined.log'
 with log.open('w') as f:
  result=subprocess.run([str(root/'refine_disks'),'parameters.txt','critical_disks.txt','refinement_indices.txt'],cwd=root,stdout=f,stderr=subprocess.STDOUT)
 text=log.read_text()
 if result.returncode or not re.search(r'^PASS REFINEMENT degree='+str(degree)+r' count='+str(len(indices))+r' ',text,re.M):raise RuntimeError('Refinement incomplete')
 refined=parse(text,'REFINED_DISK')
 if sorted(refined)!=indices:raise RuntimeError('Refinement coverage mismatch')
 for k,(lo,hi) in refined.items():
  oldlo,oldhi=rows[k];rows[k]=(max(lo,oldlo),min(hi,oldhi))
  if rows[k][0]>rows[k][1]:raise RuntimeError('Disjoint enclosures for same root')
 kmin=min(rows,key=lambda k:rows[k][1]);minimum_upper=rows[kmin][1];other_lower=min(rows[k][0] for k in rows if k!=kmin)
 summary={'status':'PASS','degree':degree,'complete_base_critical_points':degree-1,'refined_disks':len(indices),'minimum_safe_lower':str(min(x[0] for x in rows.values())),'minimum_safe_upper':str(minimum_upper),'minimizing_disk':kmin,'unique_minimizer_proved':minimum_upper<other_lower,'all_other_critical_quotients_greater_than':str(other_lower),'printed_endpoint_guard':'(1+abs(endpoint))*1e-18','seconds':time.time()-start,'compiler_flags':flags,'base_log_sha256':loghashes,'sha256':{name:digest(root/name) for name in ['polynomial.json','parameters.txt','critical_disks.txt','verify.cpp','rigorous_interval.hpp','refine_disks.cpp','refinement_indices.txt','certificate_refined.log']}}
 (root/'certificate_refined_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2),flush=True)
if __name__=='__main__':main()
