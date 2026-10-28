#!/usr/bin/env python3
"""Verify the files of the delivered snapshot before regenerating its outputs."""
from pathlib import Path
import hashlib

ROOT = Path(__file__).resolve().parent

def main() -> int:
    errors = []
    rows = (ROOT/'SHA256SUMS.txt').read_text(encoding='utf-8').splitlines()
    for row in rows:
        digest, name = row.split('  ', 1)
        p = Path(name)
        if p.is_absolute() or '..' in p.parts:
            errors.append(name + ': unsafe path')
        elif not (ROOT/p).is_file():
            errors.append(name + ': missing')
        elif hashlib.sha256((ROOT/p).read_bytes()).hexdigest() != digest:
            errors.append(name + ': changed')
    if errors:
        print('\n'.join(errors))
        return 1
    print(f'PASS: {len(rows)} delivered files match SHA256SUMS.txt')
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
