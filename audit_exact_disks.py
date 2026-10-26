#!/usr/bin/env python3
"""Exact-rational audit of disk disjointness and logged strict inequalities.

This audit independently proves disjointness of the rational containing disks.
It also checks all printed Rouche comparisons and containment statements with
an outward allowance for decimal formatting. It does not recompute the analytic
or interval enclosures; the native verifier supplies those.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import re
import time

PAIR = r'\[([^,\]]+),([^\]]+)\]'
ROW = re.compile(r'^(CRITICAL_DISK|REFINED_DISK) (\d+) radius=' + PAIR +
                 r' lhs=' + PAIR + r' rhs=' + PAIR +
                 r' quotient_lower=([^ ]+) quotient_upper=([^\n]+)$', re.M)
CENTER = re.compile(r'^FINITE_CENTER id=(\d+) center=' + PAIR +
                    r' \+ i' + PAIR + r' shift=' + PAIR +
                    r' radius=' + PAIR + r'$', re.M)
GUARD = F(1, 10**20)  # More than the maximum 22-significant-digit print loss.


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def bounds(lo: str, hi: str) -> tuple[F, F]:
    a, b = F(lo), F(hi)
    if a > b:
        raise ValueError('Reversed printed interval')
    return a - abs(a)*GUARD, b + abs(b)*GUARD


def exact_disjoint(disks: list[tuple[F, F, F]]) -> int:
    order = sorted(range(len(disks)), key=lambda j: disks[j][0])
    maxradius = max(d[2] for d in disks)
    pairs = 0
    for pos, j in enumerate(order):
        x, y, r = disks[j]
        for nextpos in range(pos+1, len(order)):
            xx, yy, rr = disks[order[nextpos]]
            dx = xx-x
            if dx > 2*maxradius:
                break
            pairs += 1
            if dx*dx+(yy-y)**2 <= (r+rr)**2:
                raise ValueError(f'Containing disks overlap: {j}, {order[nextpos]}')
    return pairs


def exact_unshifted_equality(rawlo: str, rawhi: str, radius: F, precision: int) -> bool:
    # Log endpoints are binary floating-point values with at most `precision`
    # significand bits. At a power of two r the closest other such value is
    # at least r*2^(-precision) away. The guarded decimal print interval is
    # narrower than this spacing, so matching endpoints recover r EXACTLY.
    # This is not a floating comparison or an exception to the Rouche test.
    if radius <= 0 or precision <= 0 or GUARD >= F(1, 2**precision):
        return False
    num, den = radius.numerator, radius.denominator
    if num & (num-1) or den & (den-1):
        return False
    return F(rawlo) == radius == F(rawhi)


def check_log(text: str, disks: list[tuple[F, F, F]], prefix: str,
              equal_cases: list[int] | None = None) -> list[int]:
    fmt = re.findall(r'^ARITHMETIC .*long_double_mantissa=(\d+) ', text, re.M)
    if len(fmt) != 1:
        raise ValueError('Missing or ambiguous endpoint precision')
    precision = int(fmt[0])
    centers = {}
    for match in CENTER.finditer(text):
        j = int(match[1])
        centers[j] = (bounds(match[2], match[3]), bounds(match[4], match[5]))
    ids = []
    for match in ROW.finditer(text):
        if match[1] != prefix:
            continue
        j = int(match[2])
        if not 0 <= j < len(disks):
            raise ValueError('Disk index out of range')
        rlo, rhi = bounds(match[3], match[4])
        _, lhs_hi = bounds(match[5], match[6])
        rhs_lo, _ = bounds(match[7], match[8])
        qlo = F(match[9]) - abs(F(match[9]))*GUARD
        qhi = F(match[10]) + abs(F(match[10]))*GUARD
        if not (rlo > 0 and lhs_hi < rhs_lo and qhi >= qlo > F(10000000001,10000000000)):
            raise ValueError(f'Unproved printed strict inequality: {j}')
        x,y,r = disks[j]
        if j in centers:
            rebox, imbox = centers[j]
            dx = max(abs(v-x) for v in rebox)
            dy = max(abs(v-y) for v in imbox)
        else:
            dx = dy = F(0)
        strict_containment = rhi < r and dx*dx+dy*dy < (r-rhi)**2
        equality = (j not in centers and
                    exact_unshifted_equality(match[3], match[4], r, precision))
        if not strict_containment and not equality:
            raise ValueError(f'Accepted disk containment not certified from log: {j}')
        if equality and equal_cases is not None:
            equal_cases.append(j)
        ids.append(j)
    if len(ids) != len(set(ids)):
        raise ValueError('Repeated accepted disk')
    return sorted(ids)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument('--disjoint-only', action='store_true')
    args=parser.parse_args(); root=args.directory.resolve(); start=time.monotonic()
    spec=json.loads((root/'polynomial.json').read_text())
    disks=[tuple(F(x) for x in row.split()) for row in (root/'critical_disks.txt').read_text().splitlines()]
    if len(disks) != spec['degree']-1 or any(len(d)!=3 or d[2]<=0 for d in disks):
        raise ValueError('Malformed containing disks')
    pairs=exact_disjoint(disks)
    out={'status':'PASS', 'degree':spec['degree'], 'rational_containing_disks':len(disks),
         'exact_rational_pairs_checked':pairs,'method':'Exact rational squared distances and outward-guarded printed inequalities',
         'decimal_print_guard_relative':'1e-20','source_sha256':digest(Path(__file__)),
         'input_sha256':{name:digest(root/name) for name in ['polynomial.json','critical_disks.txt']}}
    if not args.disjoint_only:
        allids=[]; equal_cases=[]
        for k in range(16):
            path=root/'logs'/f'part_{k:02d}.log'; text=path.read_text()
            if not re.search(r'^PASS mode=allroots sheet='+str(k)+r' .*degree='+str(spec['degree'])+r' ',text,re.M):
                raise ValueError('Partition not complete')
            ids=check_log(text,disks,'CRITICAL_DISK',equal_cases)
            if ids != list(range(k,spec['degree']-1,16)):
                raise ValueError('Partition coverage failure')
            allids.extend(ids)
        if sorted(allids)!=list(range(spec['degree']-1)):
            raise ValueError('Incomplete critical-set coverage')
        out['strict_rouche_and_quotient_records_checked']=len(allids)
        refine=root/'certificate_refined.log'
        if refine.exists():
            ids=check_log(refine.read_text(),disks,'REFINED_DISK',equal_cases)
            wanted=sorted(int(x) for x in (root/'refinement_indices.txt').read_text().split())
            if ids!=wanted:
                raise ValueError('Refinement coverage failure')
            out['refined_strict_records_checked']=len(ids)
            out['refinement_log_sha256']=digest(refine)
        out['exact_equal_radius_containments']=len(set(equal_cases))
        out['equal_radius_indices']=sorted(set(equal_cases))
        out['equal_radius_method']='Unshifted centers and uniquely reconstructed point powers of two; containment by exact equality. All Rouche comparisons remain strict.'
        out['warning']='This audit proves exact containing-disk separation and checks recorded inequalities; analytic enclosure validity still depends on the native interval calculation.'
    out['seconds']=time.monotonic()-start
    name='exact_disjointness_audit.json' if args.disjoint_only else 'independent_exact_audit.json'
    (root/name).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2),flush=True)

if __name__=='__main__':
    main()
