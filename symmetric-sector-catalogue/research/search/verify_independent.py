"""Independent check by symmetric multiplicities of full-Hilbert-space Paulis."""
import math
import numpy as np

def pauli_sector(c,k):
 n=len(c)-1
 idx=np.arange(2**n)
 pop=np.array([int(i).bit_count() for i in idx])
 psi=np.asarray([c[w]/math.sqrt(math.comb(n,w)) for w in pop])
 total=0.
 for x in range(k+1):
  for y in range(k-x+1):
   z=k-x-y
   flip=(1<<(x+y))-1
   yz=((1<<y)-1)<<x | ((1<<z)-1)<<(x+y)
   phase=(1j)**y * (-1.)**np.array([int(i&yz).bit_count() for i in idx])
   expect=np.vdot(psi[idx^flip],phase*psi)
   if abs(expect.imag)>1e-9:raise ValueError('Non-real Hermitian expectation')
   multiplicity=math.comb(n,k)*math.comb(k,x)*math.comb(k-x,y)
   total+=multiplicity*expect.real**2
 return total

if __name__=='__main__':
 import json,glob
 from pathlib import Path
 for file in glob.glob(str(Path(__file__).resolve().parent/'runs'/'n*-s*-*.json')):
  if any(t in file for t in ['sparse','patterns','refined']):continue
  d=json.load(open(file));c=np.array([complex(*v) for v in d['coefficients']])
  independent=pauli_sector(c,d['k']);difference=independent-d['value']
  print(d['n'],d['k'],d['mode'],format(independent,'.12f'),'difference',format(difference,'.2e'))
  assert abs(difference)<2e-7
