"""Audit direct imports and trap indirect solver calls during native tests."""
import os
from pathlib import Path
import re
import subprocess
import sys

library, shim, probe, *tests = map(lambda x: str(Path(x).resolve()),sys.argv[1:])
# Each binary's limit, sized for builds without a sanitizer; a sanitizer build multiplies it by the factor CMake
# puts in the environment (CERTPOLYSOLVE_INSTRUMENTED_SLOWDOWN in CMakeLists.txt, with its derivation).
LIMIT_SECONDS = 180
limit = LIMIT_SECONDS * int(os.environ.get('CERTPOLYSOLVE_TEST_TIME_FACTOR', '1'))
imports = subprocess.check_output(['nm','-D','--undefined-only',library],text=True)
forbidden = re.compile(r'find_roots|isolate_real_roots|arb_fmpz_poly_(?:real|complex)_roots')
assert not forbidden.search(imports), 'forbidden direct backend import'
for executable in [probe,*tests]:
    env = dict(os.environ)
    dependencies = subprocess.check_output(['ldd',executable],text=True)
    asan = re.search(r'libasan[^ ]* => (\S+)',dependencies)
    env['LD_PRELOAD'] = ':'.join(([asan.group(1)] if asan else [])+[shim])
    process = subprocess.run([executable],env=env,text=True,capture_output=True,timeout=limit)
    expected = 90 if executable==probe else 0
    assert process.returncode==expected, executable+'\n'+process.stdout+'\n'+process.stderr
    if executable==probe: assert 'forbidden backend root solver' in process.stderr
print('direct imports and runtime backend traps passed')
