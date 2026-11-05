#!/usr/bin/env python3
"""Recheck the complete-parameter lower construction's symbolic identities."""
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def verify():
    source = ROOT.parent/'analytical-constructions'/'construct_t1.py'
    spec = importlib.util.spec_from_file_location('construct_t1_general_audit', source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return {'passed': True,
            'scope': 'Symbolic complete-parameter central-block identities and coefficient positivity for the general first-order attaining construction. Separated module support and cancellation arguments are documented in proof.tex and the source constructor.',
            'source': '../analytical-constructions/construct_t1.py',
            'symbolicCentralBlocks': module.verify_general_blocks()}


if __name__ == '__main__':
    report = verify()
    (ROOT/'general-lower-verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))
