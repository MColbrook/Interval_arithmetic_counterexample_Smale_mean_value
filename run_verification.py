#!/usr/bin/env python3
"""Compile and verify all disjoint critical-point disks, with optional parallelism."""
from pathlib import Path
from fractions import Fraction
from decimal import Decimal, getcontext
from concurrent.futures import ThreadPoolExecutor,as_completed
import argparse,hashlib,json,os,platform,re,subprocess,sys,time
from make_parameters import EXPECTED,parameter_lines

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 pa=argparse.ArgumentParser(description=__doc__);pa.add_argument('--compiler',default=os.environ.get('CXX','g++'));pa.add_argument('--jobs',type=int,default=1)
 args=pa.parse_args();root=Path(__file__).resolve().parent;start=time.time();spec=json.loads((root/'polynomial.json').read_text())
 exact='\n'.join(parameter_lines(spec))+'\n'
 if (root/'parameters.txt').read_text()!=exact:raise SystemExit('FAIL: polynomial input mismatch')
 lines=(root/'critical_disks.txt').read_text().splitlines()
 if len(lines)!=EXPECTED['degree']-1:raise SystemExit('FAIL: wrong disk count')
 for line in lines:
  v=line.split()
  if len(v)!=3 or Fraction(v[2])<=0:raise SystemExit('FAIL: malformed disk data')
  for x in v:Fraction(x)
 flags=['-O2','-std=c++17','-fno-fast-math','-ffp-contract=off','-frounding-math','-DN_PARAM='+str(EXPECTED['N']),'-DBLOCKS_PARAM='+str(EXPECTED['correction_blocks']),'-DTAIL_BLOCKS_PARAM='+str(EXPECTED['tail_blocks']),'-DEXTRA_BLOCKS_PARAM='+str(EXPECTED['extra_correction_blocks'])]
 variants=[['-lmpfr','-lgmpxx','-lgmp']]
 if sys.platform.startswith('linux'):variants.append(['-l:libmpfr.so.6','-lgmpxx','-lgmp'])
 for links in variants:
  cp=subprocess.run([args.compiler,*flags,'verify.cpp','-o','verify',*links],cwd=root,text=True,capture_output=True)
  if cp.returncode==0:break
 else:raise SystemExit('Compilation failed:\n'+cp.stdout+cp.stderr)
 names=['polynomial.json','parameters.txt','critical_disks.txt','verify.cpp','rigorous_interval.hpp','make_parameters.py','run_verification.py']
 metadata={'platform':platform.platform(),'compiler':subprocess.check_output([args.compiler,'--version'],text=True).splitlines()[0],'flags':flags,'sha256':{x:digest(root/x) for x in names}}
 outdir=root/'logs';outdir.mkdir(exist_ok=True)
 def work(k):
  log=outdir/f'part_{k:02d}.log'
  with log.open('w') as f:
   p=subprocess.run([str(root/'verify'),'allroots',str(k),str(root/'parameters.txt'),str(root/'critical_disks.txt')],cwd=root,text=True,stdout=f,stderr=subprocess.STDOUT)
  text=log.read_text()
  if p.returncode or not re.search(r'^PASS mode=allroots sheet='+str(k)+r' .*degree='+str(EXPECTED['degree'])+r' ',text,re.M):raise RuntimeError(f'Part {k} did not pass; see {log.name}')
  count=len(re.findall(r'^CRITICAL_DISK ',text,re.M));print(f'PART {k:02d} PASS: {count} disks',flush=True)
  return k,text
 results={}
 with ThreadPoolExecutor(max_workers=max(1,min(16,args.jobs))) as pool:
  futures=[pool.submit(work,k) for k in range(16)]
  for future in as_completed(futures):
   k,text=future.result();results[k]=text
 alltext='CERTIFICATE_INPUTS '+json.dumps(metadata,sort_keys=True)+'\n'+''.join(results[k] for k in range(16))
 records=re.findall(r'^CRITICAL_DISK (\d+) .* quotient_lower=([^ ]+) quotient_upper=([^\n]+)$',alltext,re.M)
 ids=sorted(int(r[0]) for r in records)
 if ids!=list(range(EXPECTED['degree']-1)):raise SystemExit('FAIL: missing or repeated critical disk')
 getcontext().prec=60;lower=min(Decimal(r[1]) for r in records);upper=min(Decimal(r[2]) for r in records)
 target=Decimal(1)+Decimal(Fraction(EXPECTED['target_margin']).numerator)/Decimal(Fraction(EXPECTED['target_margin']).denominator)
 if lower<=target:raise SystemExit('FAIL: recorded lower bound is inadequate')
 alltext+=f'COMPLETE degree={EXPECTED["degree"]} checked={len(records)} minimum_printed_lower={lower} minimum_printed_upper={upper}\n'
 (root/'certificate.log').write_text(alltext)
 summary={'status':'PASS','degree':EXPECTED['degree'],'method':'One disjoint Rouche disk per critical point; count equals derivative degree','critical_points_counted_with_multiplicity':len(records),'all_critical_points_simple':True,'strict_smale_lower_bound':str(1+Fraction(EXPECTED['target_margin'])),'minimum_printed_lower':str(lower),'minimum_printed_upper':str(upper),'numeric_printing_guard':'0.000000000000001','minimum_safe_lower':str(lower-Decimal('1e-15')),'minimum_safe_upper':str(upper+Decimal('1e-15')),'seconds':time.time()-start,'inputs':metadata,'certificate_log_sha256':digest(root/'certificate.log')}
 (root/'certificate_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print('COMPLETE',json.dumps(summary,indent=2),flush=True)
if __name__=='__main__':main()
