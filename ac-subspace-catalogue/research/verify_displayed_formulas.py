#!/usr/bin/env python3
"""Parse the actual compact LaTeX and compare its coefficients exactly.
The strict parser supports only the mathematical notation used in these displays;
unrecognized notation is rejected rather than silently ignored. Requires SymPy.
"""
from pathlib import Path
import json,re,hashlib,importlib.util
import sympy as s
ROOT=Path(__file__).resolve().parents[1]

def clean(text):
 for old in [r'\begin{aligned}',r'\end{aligned}',r'\qquad',r'\quad',r'\left',r'\right',r'\bigl',r'\bigr',r'\cdot',r'\,',r'\\','&','{}']:
  text=text.replace(old,'')
 return re.sub(r'\s+','',text)

class ScalarParser:
 def __init__(self,text,env,n=None):
  self.text=text;self.env=env;self.pos=0;self.n=n;self.kets=[]
  self.dicke=[s.Symbol('D_'+str(r)) for r in range(n+1)] if n is not None else []
 def parse(self):
  value=self.sum()
  if self.pos!=len(self.text):raise ValueError(('Unexpected notation',self.text[self.pos:]))
  return s.simplify(value)
 def sum(self):
  value=self.product()
  while self.pos<len(self.text) and self.text[self.pos] in '+-':
   sign=self.text[self.pos];self.pos+=1;term=self.product();value=value+term if sign=='+' else value-term
  return value
 def product(self):
  sign=1
  while self.pos<len(self.text) and self.text[self.pos] in '+-':
   if self.text[self.pos]=='-':sign=-sign
   self.pos+=1
  factors=[]
  while self.pos<len(self.text) and self.text[self.pos] not in '+-})':
   if self.text[self.pos]=='/':
    if not factors:raise ValueError('Division without numerator')
    self.pos+=1;factors[-1]/=self.atom()
   else:factors.append(self.atom())
  if not factors:raise ValueError(('Missing scalar',self.text,self.pos))
  return sign*s.prod(factors)
 def argument(self):
  if self.pos>=len(self.text):raise ValueError('Missing command argument')
  if self.text[self.pos] in '{\\':return self.atom()
  # An unbraced TeX argument contains a single character/token: \frac14.
  token=self.text[self.pos];self.pos+=1
  if token.isdigit():return s.Integer(token)
  if token=='i':return s.I
  if token in self.env:return self.env[token]
  raise ValueError(('Unknown command argument',token))
 def atom(self):
  text=self.text
  if text[self.pos] in '{(':
   closer='}' if text[self.pos]=='{' else ')'
   self.pos+=1;value=self.sum()
   if self.pos>=len(text) or text[self.pos]!=closer:raise ValueError('Unclosed group')
   self.pos+=1;return value
  if text[self.pos]=='|' and self.n is not None:
   ket=re.match(r'\|D_(?:\{(\d+)\}|(\d))\^\{\(?(\d+)\)?\}\\rangle',text[self.pos:])
   if not ket:raise ValueError(('Unknown ket',text[self.pos:]))
   n=int(ket[1] or ket[2]);r=int(ket[3]);assert n==self.n and 0<=r<=n
   self.pos+=ket.end();self.kets.append(r);return self.dicke[r]
  for command in [r'\sqrt',r'\frac',r'\overline']:
   if text.startswith(command,self.pos):
    self.pos+=len(command)
    if command==r'\sqrt':return s.sqrt(self.argument())
    if command==r'\overline':return s.conjugate(self.argument())
    numerator=self.argument();denominator=self.argument();return numerator/denominator
  for name in sorted(self.env,key=len,reverse=True):
   if text.startswith(name,self.pos):
    self.pos+=len(name);return self.env[name]
  char=text[self.pos]
  if char.isdigit():
   found=re.match(r'\d+',text[self.pos:])[0];self.pos+=len(found);return s.Integer(found)
  self.pos+=1
  if char=='i':return s.I
  if char in self.env:return self.env[char]
  raise ValueError(('Unsupported token',char,text))

def scalar(text,env):return ScalarParser(text,env).parse()

def phases(latex):
 env={};modulus_assertions=[]
 for assignment in clean(latex).rstrip('.').split(',') if latex else []:
  lhs,rhs=assignment.split('=')
  if lhs.startswith('|'):
   name=lhs.strip('|');assert rhs=='1' and s.simplify(env[name]*s.conjugate(env[name])-1)==0
   modulus_assertions.append(name)
  else:
   assert re.fullmatch(r'[a-z]|\\(?:alpha|beta|gamma)',lhs);env[lhs]=scalar(rhs,env)
 return env,modulus_assertions

def dicke_vector(latex,n,env):
 text=clean(latex);lhs,rhs=text.split('=')
 lhs_match=re.fullmatch(r'(.*?)\|(?:u_\{?\d+\}?|[ABC])\\rangle',lhs);assert lhs_match,('lhs',lhs)
 scale=scalar(lhs_match[1],env) if lhs_match[1] else s.Integer(1)
 parser=ScalarParser(rhs,env,n);expression=s.expand(parser.parse()/scale)
 coeff=[expression.coeff(ket) for ket in parser.dicke]
 assert parser.kets and s.simplify(expression-s.Add(*(c*k for c,k in zip(coeff,parser.dicke))))==0,('nonlinear or non-ket term',expression)
 return coeff,scale,len(parser.kets)

def definitions(latex,env=None):
 env={} if env is None else env.copy()
 for assignment in clean(latex).rstrip('.').split(','):
  lhs,rhs=assignment.split('=')
  assert re.fullmatch(r'[a-z]|\\mathcalN|[pq]_(?:\{[0246]\}|[0246])',lhs),lhs
  env[lhs]=scalar(rhs,env)
 return env

def verify_new_record(rec):
 n=rec['n'];folder=ROOT/'research'/'algebraic-k3'/('n'+str(n))
 raw=(folder/'basis.json').read_bytes();source=json.loads(raw)
 if 'subspaces' in source:source=source['subspaces'][0]
 assert source['id']==rec['id'] and len(source['basis'])==len(rec['basis'])==3
 for key in ['derivationLatex','algebraicDefinition']:
  if key in source:assert rec[key]==source[key],(rec['id'],key)
 for actual,saved in zip(rec['basis'],source['basis']):
  for key in ['stateLatex','coefficientLatex','exactCoefficients']:
   if key in saved:
    target='coefficientDefinitionLatex' if key=='coefficientLatex' and isinstance(saved[key],str) else key
    assert actual[target]==saved[key],(actual['id'],target)
 # The new display cases have scalar amplitudes and algebraic variables, so
 # they must not enter the older phase-only parser.
 if n==14:
  env=definitions(rec['derivationLatex']);N=env[r'\mathcalN']
  parsed=[dicke_vector(rec['basis'][0]['stateLatex'],14,env)[0]]
  for index,(name,sign,prefactor) in enumerate([('p',1,s.I),('q',-1,s.S.One)],1):
   expected=(r'|u_'+str(index+1)+r'\rangle=\frac{'+('i' if index==1 else '1')+r'}{32\sqrt{\mathcal N}}\sum_{r\in\{0,2,4,6\}}\sqrt{\binom{14}{r}}\,'+name+r'_r\bigl(|D_{14}^{(r)}\rangle'+('+' if sign==1 else '-')+r'|D_{14}^{(14-r)}\rangle\bigr)')
   assert clean(rec['basis'][index]['stateLatex'])==clean(expected)
   terms=definitions(rec['basis'][index]['coefficientDefinitionLatex'],env);row=[s.S.Zero]*15
   for r in[0,2,4,6]:
    amplitude=prefactor*s.sqrt(s.binomial(14,r))*terms[name+'_'+str(r)]/(32*s.sqrt(N))
    row[r]=amplitude;row[14-r]=sign*amplitude
   parsed.append(row)
 elif n==15:
  symbols=s.symbols(' '.join(source['algebraicDefinition']['variables']),real=True)
  env={clean(s.latex(x)):x for x in symbols}
  parsed=[dicke_vector(b['stateLatex'],15,env)[0] for b in rec['basis']]
 else:
  assert n==17
  variables=s.symbols('x0:36',real=True);fixed={int(i):s.Rational(v)for i,v in source['algebraicDefinition']['fixedCoordinates'].items()}
  env={r'\xi_{'+str(i)+'}':fixed.get(i,x)for i,x in enumerate(variables)}
  parsed=[]
  for index,offset in enumerate([2,2,0]):
   expected=r'|u_{'+str(index+1)+r'}\rangle=\sum_{\ell=0}^{5}\frac{\sqrt{\binom{17}{3\ell+'+str(offset)+r'}}}{s_{3\ell+'+str(offset)+r'}}(\xi_{'+str(12*index)+r'+2\ell}+i\xi_{'+str(12*index+1)+r'+2\ell})|D_{17}^{(3\ell+'+str(offset)+r')}\rangle'
   assert rec['basis'][index]['stateLatex']==expected
   row=[s.S.Zero]*18
   for ell in range(6):
    r=3*ell+offset;re_index=12*index+2*ell;scale=source['algebraicDefinition']['scales'][r]
    row[r]=s.sqrt(s.binomial(17,r))*(env[r'\xi_{'+str(re_index)+'}']+s.I*env[r'\xi_{'+str(re_index+1)+'}'])/scale
   parsed.append(row)
 reports=[];count=0
 locals_={str(v):v for value in env.values() for v in value.free_symbols}
 for actual,row in zip(rec['basis'],parsed):
  exact=[s.sympify(z,locals=locals_)for z in actual['exactCoefficients']]
  assert len(row)==len(exact)
  for r,(a,b)in enumerate(zip(row,exact)):
   assert s.simplify(s.expand_complex(a-b))==0,(actual['id'],r,'displayed coefficient differs')
   if isinstance(actual.get('coefficientLatex'),list):
    assert s.simplify(scalar(clean(actual['coefficientLatex'][r]),env)-b)==0,(actual['id'],r,'component display differs')
   count+=1
  reports.append({'id':actual['id'],'coefficientsChecked':len(exact),'allDisplayedAndOmittedCoefficientsMatchExactly':True})
 # Establish that the source formulas just matched belong to a fully certified
 # exact subspace, rather than only to a matching set of decimal previews.
 spec=importlib.util.spec_from_file_location('display_source_n'+str(n),folder/'verify.py')
 module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
 proof=module.verify(write_report=False) if n==17 else module.verify()
 assert proof['passed']
 return {'id':rec['id'],'sourceSha256':hashlib.sha256(raw).hexdigest(),'displayMethod':'Exact parsing of scalar definitions and coefficient expressions; explicit finite-sum templates are expanded coefficient by coefficient. Display fields must also match the separately verified source record.','sourceExactVerifierPassed':True,'basis':reports},count

def main():
 raw=(ROOT/'data/catalogue.json').read_bytes();data=json.loads(raw);reports=[];total=0
 for rec in data['subspaces']:
  targets={'radical-n12-k2-t2','radical-n13-k2-t3','radical-n18-k2-t3','cyclic-radical-n7-k3-t1','mod3-time-reversal-n11-d2-t2','mod4-reflection-n16-d3-t2'}
  if rec['id'] not in targets and not rec['id'].startswith('molien-'):continue
  env,explicit=phases(rec.get('derivationLatex',''))
  phase_names=[name for name in env if not name.startswith('\\')]
  for name in phase_names:
   value=env[name];assert s.simplify(s.expand_complex(value*s.conjugate(value)-1))==0,(rec['id'],name,'not unit phase')
  basis_reports=[]
  for basis in rec['basis']:
   parsed,scale,terms=dicke_vector(basis['stateLatex'],rec['n'],env)
   exact=[s.sympify(c) for c in basis['exactCoefficients']]
   assert len(parsed)==len(exact)
   for r,(a,b) in enumerate(zip(parsed,exact)):
    assert s.simplify(s.expand_complex(a-b))==0,(basis['id'],r,a,b)
    total+=1
   basis_reports.append({'id':basis['id'],'normalizationFactorOnLeft':str(scale),'displayedTerms':terms,'allDisplayedAndOmittedCoefficientsMatchExactly':True,'coefficientsChecked':len(exact)})
  reports.append({'id':rec['id'],'scalarDefinitionsParsedFromDisplay':{name:str(value) for name,value in env.items()},'unitPhasesVerified':phase_names,'everyPhaseUnitModulusVerified':True,'basis':basis_reports})
 assert len(reports)==6+sum(r['id'].startswith('molien-') for r in data['subspaces']),('Missing compact-formula records',len(reports))
 for identifier in ['radical-n14-k3-t2','algebraic-n15-k3-t2','algebraic-n17-k3-t2']:
  rec=next(r for r in data['subspaces']if r['id']==identifier)
  report,count=verify_new_record(rec);reports.append(report);total+=count
 out={'passed':True,'catalogueSha256':hashlib.sha256(raw).hexdigest(),'method':'Strict independent parsing of actual displayed LaTeX, including scalar and phase definitions, Greek variables, complex conjugation, grouped ket sums, fractions and left-hand normalization factors; exact coefficient-wise comparison with saved exactCoefficients. Builder functions are not imported.','subspacesChecked':len(reports),'basisFormulasChecked':sum(len(r['basis']) for r in reports),'complexCoefficientsChecked':total,'subspaces':reports}
 Path(__file__).with_name('displayed-formulas-audit.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
