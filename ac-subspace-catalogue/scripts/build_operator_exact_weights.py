#!/usr/bin/env python3
"""Derive exact weights of rho=P_V/d from explicitly exact Dicke bases.

No decimal fitting is used. Expressions are parsed with a restricted AST;
normalization and every Gram entry are checked in exact arithmetic. Tensor
coefficients use the exact integer Hahn form of orthonormal T_lm, equivalent
to Serrano-Ensastiga/Braun (2020), Eq.36. Squared norms include both signs of m.
Requires SymPy; run with Python3 and `--ids ID...` for focused checks.
"""
from __future__ import annotations
import ast, argparse, hashlib, json, math, signal, time
from functools import lru_cache
from pathlib import Path
import sympy as sp
ROOT=Path(__file__).resolve().parents[1]

class Unsupported(ValueError): pass
class TimeLimit(TimeoutError): pass

def parse_exact(source):
    if not isinstance(source,str): raise Unsupported('Exact coefficient is not a mathematical string')
    tree=ast.parse(source,mode='eval')
    if len(list(ast.walk(tree)))>500: raise Unsupported('Expression too long')
    def walk(n):
        if isinstance(n,ast.Expression):return walk(n.body)
        if isinstance(n,ast.Constant) and isinstance(n.value,int) and not isinstance(n.value,bool):return sp.Integer(n.value)
        if isinstance(n,ast.Name) and n.id=='I':return sp.I
        if isinstance(n,ast.UnaryOp) and isinstance(n.op,(ast.UAdd,ast.USub)):return walk(n.operand)*(1 if isinstance(n.op,ast.UAdd) else -1)
        if isinstance(n,ast.BinOp):
            a,b=walk(n.left),walk(n.right)
            if isinstance(n.op,ast.Add):return a+b
            if isinstance(n.op,ast.Sub):return a-b
            if isinstance(n.op,ast.Mult):return a*b
            if isinstance(n.op,ast.Div):return a/b
            if isinstance(n.op,ast.Pow) and b.is_Integer and abs(b)<=100:return a**b
        if isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id=='sqrt' and len(n.args)==1 and not n.keywords:return sp.sqrt(walk(n.args[0]))
        raise Unsupported('Only exact integers, I, sqrt, and arithmetic are accepted')
    return walk(tree)

def simplify(x):
    x=sp.expand(x)
    if x.is_Rational:return x
    return sp.simplify(x)

@lru_cache(None)
def tensor(n,l,m,r):
    c=r+m
    S=sum((-1)**k*math.comb(n-l,k)*math.comb(l,r-k)*math.comb(l,c-k)
          for k in range(max(0,r-l,c-l),min(n-l,r,c)+1))
    if not S:return sp.S.Zero
    f=math.factorial
    ratio=sp.Rational((2*l+1)*f(l+m)*f(l-m)*f(n-r)*f(r)*f(c)*f(n-c),f(n+l+1)*f(n-l)*f(l)**2)
    return (-1)**c*S*sp.sqrt(ratio)

def exact_basis(record):
    n=record['n'];basis=[]
    if record.get('evidence')!='exact':raise Unsupported('Numerical source')
    for vector in record['basis']:
        saved=vector.get('exactCoefficients')
        if not saved or len(saved)!=n+1:raise Unsupported('No complete explicit exact basis')
        basis.append({i:parse_exact(s) for i,s in enumerate(saved) if s!='0'})
    for a,u in enumerate(basis):
        for b,v in enumerate(basis[:a+1]):
            gram=simplify(sum(sp.conjugate(x)*v.get(i,0) for i,x in u.items()))
            if gram!=int(a==b):raise Unsupported('Exact source basis is not orthonormal')
    return basis

def rho_from_basis(basis):
    rho={};d=len(basis)
    for u in basis:
        for i,x in u.items():
            for j,y in u.items():rho[i,j]=rho.get((i,j),0)+x*sp.conjugate(y)/d
    return {ij:simplify(v) for ij,v in rho.items() if v!=0}

def compact_latex(square):
    value=sp.sqrt(square)
    if square==0:return '0'
    options=[sp.latex(value)]
    if square.is_Rational:
        p,q=map(int,sp.fraction(square))
        options.append('\\sqrt{'+(str(p) if q==1 else '\\frac{'+str(p)+'}{'+str(q)+'}')+'}')
    tex=min(options,key=len)
    if len(tex)>80 or sp.count_ops(value)>18:return None
    if square.is_Rational and any(len(str(abs(int(v))))>12 for v in sp.fraction(square)):return None
    return tex

def state_weights(record):
    basis=exact_basis(record);rho=rho_from_basis(basis);n=record['n'];weights={};all_squares=[]
    bands={m:[(r,rho[r,r+m]) for r in range(n-m+1) if rho.get((r,r+m),0)!=0] for m in range(n+1)}
    for l in range(1,n+1):
        square=0
        for m in range(l+1):
            if not bands[m]:continue
            coefficient=sp.expand(sum(tensor(n,l,m,r)*z for r,z in bands[m]))
            square+=(1 if m==0 else 2)*coefficient*sp.conjugate(coefficient)
        square=simplify(square);all_squares.append(square)
        tex=compact_latex(square)
        if tex is not None:
            weights[str(l)]={'latex':tex,'squared':str(square),'numericValue':float(sp.N(sp.sqrt(square),17)),'verified':True}
    purity=simplify(sum(all_squares)+sp.Rational(1,n+1)-sp.Rational(1,len(basis)))
    if purity!=0:raise ArithmeticError('Exact multipole purity check failed')
    return {'n':n,'dimension':len(basis),'normalization':'rho=P_V/d','weights':weights,'basisOrthonormality':'exact','multipolePurity':'exact','completeRanksComputed':n}

def atomic_write(path,text):
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(text)
    temporary.replace(path)

def write(states,audit):
    # Checkpoint runs generated under an older display filter are normalized here.
    for state in states.values():
        short={}
        for rank,weight in state['weights'].items():
            tex=compact_latex(parse_exact(weight['squared']))
            if tex is not None:short[rank]={**weight,'latex':tex}
        state['weights']=short
        state['omittedLongRankCount']=state['n']-len(short)
    audit['summary']={'exactStates':len(states),'displayedWeights':sum(len(state['weights']) for state in states.values()),'nonzeroDisplayedWeights':sum(w['numericValue']>0 for state in states.values() for w in state['weights'].values()),'skippedRecords':len(audit.get('skipped',[]))}
    payload={'schemaVersion':1,'normalization':'rho=P_V/d','method':'Exact source basis, exact Gram identity, exact integer-Hahn orthonormal spin tensors, exact squared multipole norms; no fitting to decimal data.','reference':'E. Serrano-Ensastiga and D. Braun, Phys. Rev. A 101, 022332 (2020), Eqs.36-43','sourceSha256':hashlib.sha256((ROOT/'data/catalogue.json').read_bytes()).hexdigest(),'states':states}
    text=json.dumps(payload,ensure_ascii=False,indent=2)+'\n'
    atomic_write(ROOT/'data/operator-exact-weights.json',text)
    atomic_write(ROOT/'data/operator-exact-weights.js','/* Exact multipole weights derived by scripts/build_operator_exact_weights.py. */\nwindow.OperatorExactWeights='+json.dumps(payload,ensure_ascii=False,separators=(',',':'))+';\n')
    atomic_write(ROOT/'scripts/operator-exact-weights-audit.json',json.dumps(audit,ensure_ascii=False,indent=2)+'\n')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--ids',nargs='*');ap.add_argument('--resume',action='store_true');ap.add_argument('--finalize-only',action='store_true');ap.add_argument('--seconds-per-state',type=int,default=25);args=ap.parse_args()
    records=json.loads((ROOT/'data/catalogue.json').read_text())['subspaces'];states={};audit={'computed':[],'skipped':[]}
    if args.resume and (ROOT/'data/operator-exact-weights.json').exists():
        previous=json.loads((ROOT/'data/operator-exact-weights.json').read_text())
        current_sha=hashlib.sha256((ROOT/'data/catalogue.json').read_bytes()).hexdigest()
        if previous.get('sourceSha256')!=current_sha:raise SystemExit('Catalogue changed: exact checkpoints must be regenerated without --resume')
        states=previous['states']
        audit_path=ROOT/'scripts/operator-exact-weights-audit.json'
        if audit_path.exists():audit=json.loads(audit_path.read_text())
        computed={entry['id']:entry for entry in audit.get('computed',[])}
        audit['computed']=[computed.get(id,{'id':id,'ranks':len(state['weights']),'status':'verified checkpoint'}) for id,state in states.items()]
        audit.pop('summary',None)
    if args.finalize_only:
        audit['summary']={'exactStates':len(states),'displayedWeights':sum(len(state['weights']) for state in states.values()),'skippedThisRun':len(audit.get('skipped',[]))}
        write(states,audit);print(json.dumps(audit['summary']),flush=True);return
    priority={'radical-n13-k2-t3':-3,'acliege-Π_7_2_2':-2,'line-ac-n8-t3':-1}
    records=sorted(records,key=lambda r:(priority.get(r['id'],r['n']),r['id']))
    def alarm(*_):raise TimeLimit('Exact symbolic calculation exceeded per-record budget')
    signal.signal(signal.SIGALRM,alarm)
    for record in records:
        id=record['id']
        if args.ids and id not in args.ids:continue
        if id in states:continue
        started=time.monotonic();signal.alarm(args.seconds_per_state)
        try:
            result=state_weights(record);states[id]=result;audit['computed'].append({'id':id,'ranks':len(result['weights']),'seconds':round(time.monotonic()-started,3)})
            print(json.dumps({'id':id,'exactStates':len(states),'displayWeights':len(result['weights']),'seconds':round(time.monotonic()-started,2)}),flush=True)
        except (Unsupported,TimeLimit) as e:audit['skipped'].append({'id':id,'reason':str(e)})
        except Exception as e:audit['skipped'].append({'id':id,'reason':type(e).__name__+': '+str(e)});print('ERROR',id,str(e),flush=True)
        finally:signal.alarm(0)
        write(states,audit)
    audit['summary']={'exactStates':len(states),'displayedWeights':sum(len(s['weights']) for s in states.values()),'skippedThisRun':len(audit['skipped'])};write(states,audit);print(json.dumps(audit['summary']),flush=True)
if __name__=='__main__':main()
