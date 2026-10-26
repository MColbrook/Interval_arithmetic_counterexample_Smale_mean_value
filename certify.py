#!/usr/bin/env python3
"""Reproduce the complete degree-29984 interval proof."""
from pathlib import Path
import argparse, os, subprocess, sys
ROOT=Path(__file__).resolve().parent

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--jobs',type=int,default=4)
    parser.add_argument('--compiler',default=os.environ.get('CXX','g++'))
    args=parser.parse_args()
    if not 1<=args.jobs<=16:
        parser.error('--jobs must be between 1 and 16')
    commands=[
        ['make_parameters.py'],['test_containment_audit.py'],
        ['check_local.py','--compiler',args.compiler],
        ['run_verification.py','--jobs',str(args.jobs),'--compiler',args.compiler],
        ['refine_certificate.py','--compiler',args.compiler],
        ['audit_certificate.py'],['audit_exact_disks.py'],
        ['audit_normalization.py'],['check_result.py']]
    for arguments in commands:
        print('RUN:',' '.join(arguments),flush=True)
        subprocess.run([sys.executable,str(ROOT/arguments[0]),*arguments[1:]],cwd=ROOT,check=True)
    print('PASS: the complete interval computation establishes the theorem.',flush=True)

if __name__=='__main__':
    main()
