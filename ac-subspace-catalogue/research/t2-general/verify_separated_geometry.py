from pathlib import Path
import importlib.util,json,sys
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('separated_geometry',HERE.parent/'algebraic-k3/verify_geometry.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
def verify(write=False):
 spaces=json.loads((HERE/'separated-bases.json').read_text())['subspaces']
 return module.run(write=write,cases=(35,42),source_fn=lambda n:next(st for st in spaces if st['n']==n),output_dir=HERE/'separated-geometry')
if __name__=='__main__':
 (HERE/'separated-geometry').mkdir(exist_ok=True)
 print(json.dumps(verify('--write' in sys.argv),indent=2))
