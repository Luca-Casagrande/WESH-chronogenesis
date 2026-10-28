#!/usr/bin/env python3
"""Execute the explicit WESH verification inventory and record the results."""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent


def versions() -> dict:
    packages = {}
    for name in ('numpy', 'scipy', 'sympy', 'mpmath'):
        try:
            packages[name] = version(name)
        except PackageNotFoundError:
            packages[name] = 'not installed'
    return {'python': platform.python_version(), 'platform': platform.platform(),
            'packages': packages}


def run_one(relative: str, timeout: int) -> dict:
    path = ROOT / relative
    start = time.monotonic()
    env = os.environ.copy()
    env.pop('PYTHONOPTIMIZE', None)
    env.update({'PYTHONUTF8': '1', 'PYTHONIOENCODING': 'utf-8',
                'PYTHONDONTWRITEBYTECODE': '1',
                'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1'})
    result = {'script': relative}
    try:
        p = subprocess.run([sys.executable, '-B', str(path)], cwd=path.parent,
                           env=env, stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT, encoding='utf-8',
                           errors='replace', timeout=timeout)
        output = p.stdout
        result['returncode'] = p.returncode
    except subprocess.TimeoutExpired as exc:
        output = exc.stdout or ''
        if isinstance(output, bytes):
            output = output.decode('utf-8', errors='replace')
        result.update(returncode=-1, error='timeout')
    except OSError as exc:
        output = str(exc)
        result.update(returncode=-2, error=str(exc))
    result['seconds'] = round(time.monotonic() - start, 3)
    result['passed'] = result['returncode'] == 0
    path.with_suffix('.log').write_text(output, encoding='utf-8')
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--jobs', type=int, default=1,
                        help='Parallel workers (default: 1).')
    parser.add_argument('--timeout', type=int, default=300,
                        help='Timeout per program, in seconds (default: 300).')
    args = parser.parse_args()
    if args.jobs < 1 or args.timeout < 1:
        parser.error('jobs and timeout must be positive')
    if not __debug__:
        parser.error('Use Python without -O or -OO: assertions are part of the tests.')
    manifest = json.loads((ROOT/'verification_manifest.json').read_text(encoding='utf-8'))
    paths = manifest['scripts']
    if not paths or len(paths) != len(set(paths)):
        parser.error('Verification inventory is empty or has duplicates.')
    for rel in paths:
        p = Path(rel)
        if p.is_absolute() or '..' in p.parts or not (ROOT/p).is_file():
            parser.error(f'Missing or invalid verification path: {rel}')
    started = datetime.now(timezone.utc).isoformat()
    start = time.monotonic()
    results = []
    if args.jobs == 1:
        for rel in paths:
            item = run_one(rel, args.timeout)
            results.append(item)
            print(f"{'PASS' if item['passed'] else 'FAIL'}  {rel} ({item['seconds']:.2f}s)", flush=True)
    else:
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = {pool.submit(run_one, p, args.timeout): p for p in paths}
            for future in as_completed(futures):
                item = future.result(); results.append(item)
                print(f"{'PASS' if item['passed'] else 'FAIL'}  {item['script']} ({item['seconds']:.2f}s)", flush=True)
    environment = versions()
    results.sort(key=lambda x: x['script'])
    passed = sum(item['passed'] for item in results)
    summary = {'package_version': manifest['package_version'], 'started_utc': started,
               'completed_utc': datetime.now(timezone.utc).isoformat(),
               'duration_seconds': round(time.monotonic()-start, 3),
               'environment': environment, 'script_count': len(paths),
               'passed_count': passed, 'all_passed': passed == len(paths),
               'scripts': results,
               'optional_published_pdf_check': 'not run by this command'}
    (ROOT/'run_summary.json').write_text(json.dumps(summary, indent=2)+'\n', encoding='utf-8')
    (ROOT/'tested_environment.json').write_text(json.dumps(environment, indent=2)+'\n', encoding='utf-8')
    print(f"\n{passed}/{len(paths)} programs passed; see run_summary.json.")
    return 0 if summary['all_passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
