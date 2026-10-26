#!/usr/bin/env python3
"""Rebuild and check disk arithmetic, the minimizing critical disk, and the model-anchor comparison."""
from pathlib import Path
import argparse, os, subprocess, sys
from make_parameters import EXPECTED
ROOT=Path(__file__).resolve().parent

def compile_program(compiler,source,binary):
    flags=['-O2','-std=c++17','-fno-fast-math','-ffp-contract=off','-frounding-math',
           '-DN_PARAM='+str(EXPECTED['N']),'-DBLOCKS_PARAM='+str(EXPECTED['correction_blocks']),
           '-DTAIL_BLOCKS_PARAM='+str(EXPECTED['tail_blocks']),
           '-DEXTRA_BLOCKS_PARAM='+str(EXPECTED['extra_correction_blocks'])]
    variants=[['-lmpfr','-lgmpxx','-lgmp']]
    if sys.platform.startswith('linux'):
        variants.append(['-l:libmpfr.so.6','-lgmpxx','-lgmp'])
    for libraries in variants:
        result=subprocess.run([compiler,*flags,source,'-o',binary,*libraries],
                              cwd=ROOT,text=True,capture_output=True)
        if result.returncode==0:
            return
    raise RuntimeError('Compilation failed:\n'+result.stdout+result.stderr)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler',default=os.environ.get('CXX','g++'))
    args=parser.parse_args()
    out=ROOT/'local_checks';out.mkdir(exist_ok=True)
    tasks=[('test_ball_arithmetic.cpp','test_ball_arithmetic',[],out/'arithmetic.log'),
           ('verify.cpp','verify',['disk','28891','parameters.txt','critical_disks.txt'],out/'minimizing_disk.log'),
           ('check_base_comparison.cpp','check_base_comparison',
            ['parameters.txt','critical_disks.txt','27158'],ROOT/'base_comparison.log')]
    for source,binary,arguments,log in tasks:
        print('BUILD:',source,flush=True)
        compile_program(args.compiler,source,binary)
        with log.open('w',encoding='utf-8') as stream:
            subprocess.run([str(ROOT/binary),*arguments],cwd=ROOT,
                           stdout=stream,stderr=subprocess.STDOUT,check=True)
        print('PASS:',binary,'(see '+str(log.relative_to(ROOT))+')',flush=True)
    print('PASS: local rebuild checks completed.',flush=True)

if __name__=='__main__':
    main()
