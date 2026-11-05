"""Standalone rigorous proof verifier. No optimizer or numerical tolerances.
Rebuilds every intermediate in a new temporary directory. A C compiler is
required for the supplied arithmetic-only FLINT accelerator.
"""
from pathlib import Path
import sys,os,shutil,tempfile,subprocess,json,hashlib
import flint
BASE=Path(__file__).resolve().parent
INPUTS=['source-basis.json','source-kernels.json','holomorphic-basis.json','covariants.json','sos-certificate.json','state-parameters.json','exact_field.py','exact_checks.py','integer_evaluations.py','basis_values.py','interval_evaluations.py','arb_helper.c','verify_final.py','check_couplings.py','verify_certificate.py','requirements.txt','value-definition.json']
assert flint.__version__=='0.9.0','Install the pinned dependencies in requirements.txt.'
package=Path(flint.__file__).resolve().parent
libraries=list((package/'.dylibs').glob('libflint*.dylib'))+list((package.parent/'python_flint.libs').glob('libflint*.so*'))+list((package.parent/'flint.libs').glob('libflint*.so*'))
assert libraries,'Cannot locate the FLINT shared library bundled with python-flint.'
compiler=shutil.which('cc')or shutil.which('clang')or shutil.which('gcc');assert compiler,'A C compiler is required.'
env=dict(os.environ);env.pop('BEGIN_POINT',None);env.pop('STOP_POINT',None);env['OPENBLAS_NUM_THREADS']='1'
with tempfile.TemporaryDirectory(prefix='n11-s6-exact-proof-')as name:
 root=Path(name)
 for file in INPUTS:shutil.copy2(BASE/file,root/file)
 suffix='.dylib'if sys.platform=='darwin'else'.so';lib=libraries[0]
 subprocess.run([compiler,'-O3','-shared','-fPIC',str(root/'arb_helper.c'),str(lib),'-Wl,-rpath,'+str(lib.parent),'-o',str(root/('arb_helper'+suffix))],check=True,env=env)
 for args in [['check_couplings.py'],['exact_checks.py'],['integer_evaluations.py','source'],['integer_evaluations.py','holo'],['basis_values.py'],['interval_evaluations.py'],['verify_final.py']]:
  print('VERIFY',*args,flush=True);subprocess.run([sys.executable,str(root/args[0]),*args[1:]],check=True,cwd=root,env=env)
 report=json.loads((root/'independent-verification.json').read_text());assert report['globalOptimalityProven']and report['fullPolynomialIdentityVerified']
 for file,digest in report['inputSha256'].items():assert hashlib.sha256((BASE/file).read_bytes()).hexdigest()==digest
 shutil.copy2(root/'independent-verification.json',BASE/'independent-verification.json');shutil.copy2(root/'pivot-enclosures.json',BASE/'pivot-enclosures.json')
print('PASS: exact global minimum of eleven-qubit S6, with exact analytical attainment.',flush=True)
