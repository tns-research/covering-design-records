#!/usr/bin/env python3
"""Check that a file of blocks is a (v,k,t) covering design. usage: verify.py C_89_8_3_2606.txt 89 8 3"""
import itertools, math, sys
path, v, k, t = sys.argv[1], *map(int, sys.argv[2:5])
blocks = [sorted(map(int, l.split())) for l in open(path) if l.strip()]
assert all(len(set(b)) == k and b[0] >= 1 and b[-1] <= v for b in blocks), "invalid block"
seen = set()
for b in blocks:
    seen.update(itertools.combinations(b, t))
total = math.comb(v, t)
print(f"{len(blocks)} blocks, {len(seen)} of {total} {t}-subsets covered")
sys.exit(0 if len(seen) == total else 1)
