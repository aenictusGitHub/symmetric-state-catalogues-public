#!/usr/bin/env python3
"""Reproduce the complete finite search and independent numerical verification."""
from search import search
from summarize import main
from verify_independently import verify
JOBS=[(24,7,1,450),(21,4,1,700),(22,5,1,700),(23,5,1,700),(21,5,500,700),(22,6,500,700),(23,6,500,700),(21,6,160,700),(21,7,160,700),(22,7,160,700),(23,7,160,700)]
if __name__=='__main__':
 for n,t,restarts,steps in JOBS:search(n,t,restarts,steps,211103+100*n+t)
 main();print(verify()['passed'])
