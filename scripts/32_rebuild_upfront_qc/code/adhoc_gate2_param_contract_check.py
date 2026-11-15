"""READ-ONLY DIAGNOSTIC -- NOT A TEST, and it does not execute Scrublet.

Verifies the Scrublet parameter contract against the installed scanpy.

Checks three things `inspect.signature` alone cannot:
  1. every SCRUBLET_PARAMS key is a real wrapper argument;
  2. each value equals the wrapper's own default, so pinning them is equivalent
     to the original call rather than a silent parameter change;
  3. where each parameter is forwarded inside the wrapper body -- a name can
     exist in the signature and still not reach the underlying Scrublet.

Writes nothing, runs no scoring.
"""
import inspect
import os
import re
import sys

import scanpy as sc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gate2_score_core import SCRUBLET_PARAMS  # noqa: E402

sig = inspect.signature(sc.pp.scrublet).parameters
src = inspect.getsource(sc.pp.scrublet)

print(f'scanpy {sc.__version__} | SCRUBLET_PARAMS has {len(SCRUBLET_PARAMS)} keys\n')
print(f'{"key":<36}{"ours":<16}{"scanpy default":<18}{"equal":<7}in-sig')
print('-' * 84)
mismatched, absent = [], []
for k, v in SCRUBLET_PARAMS.items():
    in_sig = k in sig
    d = sig[k].default if in_sig else '<ABSENT>'
    same = (v == d) if in_sig else False
    if not in_sig:
        absent.append(k)
    elif not same:
        mismatched.append((k, v, d))
    print(f'{k:<36}{str(v):<16}{str(d):<18}{str(same):<7}{in_sig}')

print('\nexpected_doublet_rate present in SCRUBLET_PARAMS:',
      'expected_doublet_rate' in SCRUBLET_PARAMS, '(must be False -- passed separately)')
print('random_state present in SCRUBLET_PARAMS         :',
      'random_state' in SCRUBLET_PARAMS, '(must be False -- passed separately)')
print('batch_key present in SCRUBLET_PARAMS            :',
      'batch_key' in SCRUBLET_PARAMS,
      '(must be False -- this layer scores one capture at a time)')

print('\n--- forwarding inside the wrapper body ---')
for k in list(SCRUBLET_PARAMS) + ['expected_doublet_rate', 'random_state', 'batch_key']:
    lines = [ln.strip() for ln in src.splitlines() if re.search(r'\b' + re.escape(k) + r'\b', ln)]
    forwarded = [ln for ln in lines if '=' in ln and not ln.startswith(k + ':')]
    sample = forwarded[-1] if forwarded else (lines[-1] if lines else '-')
    print(f'  {k:<36}{len(lines):>3} mentions | {sample[:64]}')

print('\nVERDICT')
print('  keys absent from the wrapper signature :', absent or 'none')
print('  values differing from scanpy defaults  :',
      [(k, o, d) for k, o, d in mismatched] or 'none')
