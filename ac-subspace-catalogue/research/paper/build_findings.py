from pathlib import Path
import sympy as s
import json, math, hashlib
import numpy as np
ROOT=Path(__file__).resolve().parent
I=s.I

def sparse(n,d): return [s.sympify(d.get(r,0)) for r in range(n+1)]
def record(id,n,t,rows,page,equation,label,kind='paper-exact'):
 rows=[sparse(n,row) for row in rows]
 return {'id':id,'n':n,'j':str(s.Rational(n,2)),'dimension':len(rows),'order':t,'label':label,'kind':kind,'basisConvention':'Normalized Dicke basis, r=0,...,N; the paper spin basis is ordered m=j,...,-j, so r=j-m.','exactBasis':[[str(s.simplify(z)) for z in row] for row in rows],'basis':[[[float(s.re(z).evalf()),float(s.im(z).evalf())] for z in row] for row in rows],'reference':{'pdfPage':page,'printedPage':f'022435-{page}','equation':equation}}

records=[]
records.append(record('paper-eq42',4,1,[{0:s.Rational(1,2),2:I/s.sqrt(2),4:s.Rational(1,2)},{0:s.Rational(1,2),2:-I/s.sqrt(2),4:s.Rational(1,2)}],6,'42','The unique spin-2 two-dimensional 1-AC subspace, up to rotation'))
records.append(record('paper-d1',5,1,[{0:s.sqrt(s.Rational(3,8)),4:s.sqrt(s.Rational(5,8))},{1:-s.sqrt(s.Rational(5,8)),5:s.sqrt(s.Rational(3,8))}],11,'D1','Spin-5/2 subspace V1'))
a=s.sqrt(20-5*s.sqrt(7));b=s.sqrt(10*s.sqrt(7)-20);c=s.sqrt(16-5*s.sqrt(7))
records.append(record('paper-d5',5,1,[{1:a/4,3:-b/4,5:c/4},{0:c/4,2:b/4,4:a/4}],11,'D5','Spin-5/2 subspace V2'))
records.append(record('paper-d7',6,1,[{0:s.sqrt(s.Rational(2,5)),5:s.sqrt(s.Rational(3,5))},{1:-s.sqrt(s.Rational(3,5)),6:s.sqrt(s.Rational(2,5))},{3:1}],11,'D7','Three-dimensional spin-3 subspace'))
phase=(s.sqrt(5)+2*I*s.sqrt(7))/s.sqrt(33)
records.append(record('paper-d10',9,1,[{3:s.sqrt(s.Rational(5,8)),7:-s.sqrt(s.Rational(3,8))},{2:-s.sqrt(s.Rational(3,8)),6:-s.sqrt(s.Rational(5,8))},{0:s.sqrt(s.Rational(11,2))*phase/4,4:s.sqrt(3)/4,8:s.sqrt(s.Rational(15,2))/4},{1:s.sqrt(s.Rational(15,2))/4,5:-s.sqrt(3)/4,9:s.sqrt(s.Rational(11,2))*phase/4}],12,'D10 (begins on PDF p. 11)','Four-dimensional spin-9/2 subspace'))
records.append(record('paper-e1',7,2,[{0:s.sqrt(s.Rational(3,10)),5:s.sqrt(s.Rational(7,10))},{2:s.sqrt(s.Rational(7,10)),7:-s.sqrt(s.Rational(3,10))}],12,'E1','Two-dimensional spin-7/2 subspace'))
records.append(record('paper-e2',10,2,[{0:s.sqrt(s.Rational(2,7)),7:s.sqrt(s.Rational(5,7))},{3:-s.sqrt(s.Rational(5,7)),10:s.sqrt(s.Rational(2,7))}],12,'E2','Two-dimensional spin-5 subspace'))
# Equivalent published bases used in the uniqueness/Schmidt-decomposition proofs.
records.append(record('paper-c7',4,1,[{2:1},{0:1/s.sqrt(2),4:1/s.sqrt(2)}],11,'C7','Dicke / spin-cat basis of the spin-2 subspace'))
records.append(record('paper-f7',4,1,[{0:1/s.sqrt(3),3:s.sqrt(s.Rational(2,3))},{1:s.sqrt(s.Rational(2,3)),4:-1/s.sqrt(3)}],12,'F7','Rotated spin-2 subspace basis'))

# General family (46), retaining the explicit admissibility condition for (47).
def family1(n):
 j=s.Rational(n,2);rows=[]
 for a in range(max(0,int(s.floor((j-1)/2))+1)):
  rows.append({2*a:1/s.sqrt(2),n-2*a:1/s.sqrt(2)})
 include_center=n%2==0 and 2*j-3-4*s.floor((j-1)/2)>0
 if include_center:rows.append({n//2:1})
 return rows
for n in range(2,21):
 records.append(record(f'family46-n{n}',n,1,family1(n),8,'46–47 (Eq. 46 starts on PDF p. 7)','Separated opposite-spin pairs','paper-family-verified'))
# Exact support-based interpretation of (51), avoiding a noninteger dimension.
def family2(n):
 j=s.Rational(n,2);target=j*(j+1)/3;kappa=int(s.ceiling(j-s.sqrt(target)))
 rows=[];occupied=[]
 for a in range(n):
  hi=j-3*a;lo=j-3*a-1-kappa
  indices=[3*a,3*a+kappa+1,n-3*a-kappa-1,n-3*a]
  if min(indices)<0 or max(indices)>n or len(set(indices))<4:break
  merged=sorted(occupied+indices)
  if any(y-x<3 for x,y in zip(merged,merged[1:])):break
  alpha2=s.simplify((target-lo*lo)/(2*(hi*hi-lo*lo)))
  beta2=s.simplify((hi*hi-target)/(2*(hi*hi-lo*lo)))
  if not(alpha2>0 and beta2>0):break
  rows.append({indices[0]:s.sqrt(alpha2),indices[1]:s.sqrt(beta2),indices[2]:s.sqrt(beta2),indices[3]:s.sqrt(alpha2)})
  occupied=merged
 return rows,kappa
for n in list(range(10,21))+[22,35]:
 rows,kappa=family2(n)
 if rows:
  r=record(f'family51-n{n}',n,2,rows,8,'50–51; Table III for N=10,22,35','Separated four-support second-order family','paper-family-verified')
  if len(rows)==1 and 12<=n<=20:
   r['paperConstructionOrder']=2;r['order']=3;r['independentOrderUpgrade']='Symmetric probabilities cancel odd diagonal moments, and all occupied Dicke indices are separated by at least four. Independently verified through order three.'
  r['kappa']=kappa;r['dimensionRule']='Count the explicitly generated nonnegative integer a values with positive weights and all global Dicke-support separations at least three.'
  records.append(r)
# An immediate family extending the patterns in E1 and E2.
for n in [13,16,19]:
 j=s.Rational(n,2);q=(j+1)/3;split=int(j+q)
 rows=[{0:s.sqrt(q/(j+q)),split:s.sqrt(j/(j+q))},{n-split:-s.sqrt(j/(j+q)),n:s.sqrt(q/(j+q))}]
 r=record(f'pair-extension-n{n}',n,2,rows,12,'Independent extension of the pattern in E2','Two-state second-order family for N≡1 mod 3','independent-corollary')
 r['proof']='Each basis vector has mean Jz zero and second moment j*q=j(j+1)/3, q=(j+1)/3. For N≥10 and N≡1 mod3, all occupied Dicke indices in the pair are separated by at least three, so every rank≤2 off-diagonal matrix element vanishes.'
 records.append(r)

# Explicit sparse pure 3-AC family requested for the full N=12,...,20 range.
for n in range(12,21):
 j=s.Rational(n,2);a=max(4,math.ceil(n/4));inner=j-a
 p=s.simplify((j*(j+1)/3-inner*inner)/(2*(j*j-inner*inner)))
 q=s.Rational(1,2)-p
 assert p>0 and q>0 and min(a,n-2*a)>=4
 r=record(f'sparse-third-n{n}',n,3,[{0:s.sqrt(p),a:s.sqrt(q),n-a:s.sqrt(q),n:s.sqrt(p)}],None,'Independent sparse-moment construction','Sparse four-support pure 3-AC state','independent-corollary')
 r['reference']={'construction':'Independent sparse-moment construction; uses the paper definition and Eq. 8 tensor structure.'}
 r['parameters']={'a':a,'endpointProbabilityPerPoint':str(p),'innerProbabilityPerPoint':str(q)}
 r['proof']='All occupied Dicke indices are at least four apart, eliminating every off-diagonal tensor component through rank3. Reflection-symmetric probabilities cancel Jz and Jz³. The chosen p,q enforce 2p*j²+2q*(j−a)²=j(j+1)/3, so every diagonal rank1–3 multipole vanishes.'
 records.append(r)

# Rounded E3/E4 inputs are not exact states and are not assumed normalized.
values={'a':[-.0929878947,.6014799921],'b':[.7239719068,.1184090590],'c':[.3285593938,-.1799267606],'d':[-.0455578468,-.4480124508],'e':[.5237199686,.0902146965],'f':[.6915079166,.0880343260],'g':[-.0083833249,.7694216889],'h':[-.1396861199,-.5193356365]}
rows=[]
for support in [{2:'a',5:'b',8:'c',11:'d',14:'e'},{0:'e',3:'d',6:'c',9:'b',12:'a'},{1:'f',4:'g',7:'h',10:'g',13:'f'}]:
 vector=np.array([complex(*values[support[r]]) if r in support else 0j for r in range(15)]);rows.append(vector)
numerical={'id':'paper-e3-rounded','n':14,'j':'7','dimension':3,'order':2,'kind':'paper-numerical','label':'Rounded spin-7 three-dimensional second-order example','reference':{'pdfPage':12,'equations':['E3','E4']},'rawBasis':[[[z.real,z.imag] for z in row] for row in rows],'rawSquaredNorms':[float(np.vdot(row,row).real) for row in rows],'basis':[[[z.real,z.imag] for z in row/np.linalg.norm(row)] for row in rows],'claim':'The printed decimals support a numerical witness only. Normalize before checking; no exact existence or maximality proof follows from rounding.'}

findings={'source':{'title':'Quantum metrology of rotations with mixed spin states','authors':['Eduardo Serrano-Ensástiga','Chryssomalis Chryssomalakos','John Martin'],'journal':'Physical Review A','volume':111,'article':'022435','year':2025,'doi':'10.1103/PhysRevA.111.022435','url':'https://doi.org/10.1103/PhysRevA.111.022435','published':'2025-02-20','received':'2024-07-09','accepted':'2025-01-31','pdfPages':14,'localSource':'/Users/jmartin/Documents/Latex/Mypapers/PRA111_022435(2025).pdf'},
'definition':{'n':'N=2j qubits','k':'Complex dimension of the subspace, not number of displayed states with individual AC properties.','t':'Every normalized vector in the whole subspace is at least t-anticoherent.','paperReference':{'pdfPage':3,'definition':2,'lemma':1,'equations':['10','11','12']},'compressionCondition':'B^† T_{LM} B=0 for all 1≤L≤t and −L≤M≤L.','projectorCondition':'P_V T_{LM} P_V=0.','reducedCrossCondition':'Tr_{N−t}(|ψ_a><ψ_b|)=δ_ab I_{Sym^t}/(t+1) for an orthonormal basis.','generalBasisCondition':'For a nonorthonormal full-rank basis, Tr_{N−t}(|v_a><v_b|)=<v_b|v_a> I/(t+1).','warning':'Checking only the diagonal reduced states is insufficient; all coherent superpositions require the cross conditions.'},
'bounds':[{'id':'paper-schmidt-bound','formulaLatex':r'k(t+1)\le N-t+1,\qquad k\le\left\lfloor\frac{N-t+1}{t+1}\right\rfloor','status':'proved','references':[{'pdfPage':9,'equation':'52'},{'pdfPage':12,'lemma':'F1','equation':'F3'},{'pdfPage':13,'location':'Paragraph before proof of Theorem 1'}],'proof':'Choose a common symmetric t-qubit basis in the maximally entangled Schmidt decomposition of every orthonormal code vector. The cross-reduction condition makes all k(t+1) partner vectors mutually orthogonal in Sym^(N−t)(C²), which has dimension N−t+1.'},
{'id':'tensor-gram-bound','formulaLatex':r'k\left(\left\lfloor\frac t2\right\rfloor+1\right)^2\le N+1','status':'proved-independent-corollary','proof':'Let s=floor(t/2), and choose Hilbert–Schmidt orthonormal spin tensor operators A_α spanning ranks 0,...,s. There are (s+1)² such operators. Every A_α†A_β is a linear combination of tensor ranks ≤2s≤t, so its compression to V is Tr(A_α†A_β) I_V/(N+1)=δ_αβ I_V/(N+1). Therefore the k(s+1)² vectors sqrt(N+1) A_α|ψ_a> are orthonormal in H^(j), of dimension N+1.','attribution':'An independent consequence of the paper’s compression characterization and standard spin-tensor product decomposition; no novelty claim.','consequence':'For t=2, k≤floor((N+1)/4). With E2 this proves the maximum dimension is exactly2 at N10,t2; the paper bound alone is3.'}],
'figure2':{'pdfPage':8,'status':'numerically found dimensions, not rigorous maxima','n':list(range(1,18)),'t1':[0,1,1,2,2,3,3,3,4,4,4,5,5,5,6,6,6],'t2':[0,0,0,1,0,1,2,1,1,2,2,2,2,3,3,3,3],'notes':'Read from the complete figure. A reported zero is not a nonexistence proof; use independent bounds/exclusions when available.'},
'families':[{'id':'paper-pair-family','equations':['46','47','48'],'pdfPages':[7,8],'basisLatex':r'|\psi_a\rangle=(|D_N^{(2a)}\rangle+|D_N^{(N-2a)}\rangle)/\sqrt2','aRange':'a=0,...,floor((j−1)/2)','centerCondition':'Append |D_N^(N/2)> only for integer j with 2j−3−4floor((j−1)/2)>0, equivalently positive even integer j.','verifiedDimension':'floor((j−1)/2)+1+1_{j positive even integer}. For integer j≥1, this equals floor(j/2)+1.','proof':'Different supports are separated by at least two; every basis vector balances opposite m, so Jz compression vanishes and all J± compression entries vanish.'},
{'id':'paper-four-support-family','equations':['49','50','51'],'pdfPage':8,'kappa':'ceil(j−sqrt(j(j+1)/3))','support':'r=3a,3a+κ+1,N−3a−κ−1,N−3a','weights':'Let u=j−3a, v=j−3a−κ−1, μ=j(j+1)/3. α²=(μ−v²)/(2(u²−v²)), β²=(u²−μ)/(2(u²−v²)); amplitudes α,β,β,α.','admissibility':'Use only nonnegative integer a with all four indices in [0,N], positive weights, and global support separations at least3. Direct verification determines the emitted dimension.','proof':'All supports are at least three apart, eliminating every rank≤2 ladder term, including cross terms. Opposite m pairs eliminate the first diagonal moment, and α,β solve the second diagonal moment exactly.'}],
'sourceCaveats':[{'id':'eq48-center-condition','pdfPages':[7,8],'equations':['47','48'],'issue':'The integer-j dimension in printed Eq48 omits the explicit condition for adding the middle vector in the preceding text. At j3, the pair with m=±1 has nonzero J+ cross matrix element with m0.','handling':'Retain the Eq47 support condition and verify all matrix elements.'},{'id':'eq49-noninteger','pdfPage':8,'equation':'49','issue':'The printed minimum contains a floor on its first argument only. At j10, κ4 gives literal k2=5/3.','handling':'Do not assign a fractional dimension or claim an intended correction; count admissible integer basis vectors and verify the subspace condition directly.'},{'id':'eq13-upper-triangle','pdfPage':3,'equations':['13','14'],'issue':'Eq13 prints the sum m1≤m2, while the Hilbert–Schmidt trace in Eq14 sums all ordered matrix entries. Their numerical values are not generally equal. The simultaneous zero condition after summing all M remains the relevant characterization.','handling':'Use the complete matrix compression Frobenius norm, including both off-diagonal directions.'},{'id':'e3-rounded-unnormalized','pdfPage':12,'equations':['E3','E4'],'issue':'Basis coefficients are approximate and the listed vectors are not normalized. Lemma1 permits a nonorthonormal basis, but projector/purity calculations need normalization or orthonormalization.','handling':'Keep the raw decimals with provenance and verify the normalized entire span numerically.'}],
'exactWitnesses':records,'numericalWitnesses':[numerical]}
(ROOT/'findings.json').write_text(json.dumps(findings,ensure_ascii=False,indent=2)+'\n')
print('exact witnesses',len(records),'numeric',len(findings['numericalWitnesses']))
