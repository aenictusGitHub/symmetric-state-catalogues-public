#!/usr/bin/env python3
"""Certify algebraic Majorana points with exact rational root disks."""
from pathlib import Path
import importlib.util,json,sys
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('nested_exact_majorana_root_disks',HERE.parents[1]/'algebraic-k3/verify_geometry.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
def verify(write=False):
 return module.run(write=write,cases=(42,),source_fn=lambda n:json.loads((HERE/'basis.json').read_text()),output_dir=HERE)
if __name__=='__main__':print(json.dumps(verify('--write' in sys.argv),indent=2))
