#!/usr/bin/env python3
"""Check exact inputs and the supplied complete certificate without recomputing its interval calculations."""
from pathlib import Path
import subprocess, sys
ROOT=Path(__file__).resolve().parent

def main():
    for script in ('make_parameters.py','test_containment_audit.py','audit_certificate.py',
                   'audit_exact_disks.py','audit_normalization.py','check_result.py'):
        print('CHECK:',script,flush=True)
        subprocess.run([sys.executable,str(ROOT/script)],cwd=ROOT,check=True)
    print('PASS: stored certificate records are consistent and imply the published bounds.',flush=True)
    print('The native interval calculations were not recomputed; use certify.py to rerun them.',flush=True)

if __name__=='__main__':
    main()
