#!/usr/bin/env python3
"""Verify the exact five-dimensional N=13 first-order AC lower basis."""
import json
from pathlib import Path
from verify import verify_lower

if __name__ == '__main__':
    result=verify_lower()
    Path(__file__).with_name('lower-verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
