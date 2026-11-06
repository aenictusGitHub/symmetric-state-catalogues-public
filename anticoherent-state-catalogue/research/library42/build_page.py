"""Render the library audit page from current imported states and verification.

The historical library42 directory name is retained for existing links.
"""
from pathlib import Path
from decimal import Decimal
import html,json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
states=json.loads((HERE/'states.json').read_text())['states']
exact=json.loads((HERE/'exact-lower-states.json').read_text())
report=json.loads((HERE/'verification.json').read_text())
catalogue=json.loads((ROOT/'data/catalogue.json').read_text())
config_path=HERE/'config.json'
config=json.loads(config_path.read_text()) if config_path.exists() else {'maxQubits':100}
max_n=max(st['n'] for st in states)
configured_n=config.get('maxQubits',100)
catalogue_numbers=catalogue.get('coreNumbers',[])+catalogue.get('extendedNumbers',[])
first_excluded=next(t for t in range(1,configured_n+2)
                    if ((t//2+1)**2 if t%2==0 else (t//2+1)*(t//2+2))>configured_n)
max_order=max(catalogue.get('orders',[first_excluded])) if max(catalogue_numbers,default=0)>=configured_n else first_excluded

def decimal(value):
    return Decimal(str(value))

def display(value):
    return format(value,'.4E').replace('E-','e−').replace('E+','e+')

def source_id(state):
    return state.get('numerical',{}).get('sourceRecordId') or state.get('sourceRecordId') or state['id'].removeprefix('liege-')

audits={r['id']:r for r in report.get('records',[])}

def reductions(state):
    # Additional refined records carry their independent audit with the state.
    return (audits.get(source_id(state)) or state['verification']).get('reductions',[])

def through_residual(state):
    values=[decimal(r['frobeniusDistance']) for r in reductions(state) if r['t']<=state['achievedOrder']]
    if not values:values=[decimal(state['verification']['maximumResidual'])]
    return max(values)

through_max=max(through_residual(state) for state in states)
next_values=[decimal(r['frobeniusDistance']) for state in states for r in reductions(state) if r['t']==state['achievedOrder']+1]
next_min=min(next_values) if next_values else None
precision=report.get('workingPrecisionDigits',100)
rows=[]
for n in range(2,max_n+1):
    available=[st for st in states if st['n']==n]
    if not available:
        rows.append(f'<tr><th scope="row">{n}</th><td>Not recorded</td><td>—</td></tr>')
        continue
    top=max(st['achievedOrder'] for st in available)
    best=min((st for st in available if st['achievedOrder']==top),key=through_residual)
    rows.append(f'<tr><th scope="row">{n}</th><td>{top}</td><td>{display(through_residual(best))}</td></tr>')

next_text=f'; the smallest residual at the next order is {display(next_min)}' if next_min is not None else ''
additional_n60='<p>The numerical N=60, order-10 representative is included with the same independent reduced-state checks. Its reported order is not a proof of maximality.</p>' if any(st['n']==60 and st['achievedOrder']==10 for st in states) else ''
extra100_path=HERE/'extra-states-100.json'
extra100=json.loads(extra100_path.read_text())['states'] if extra100_path.exists() else []
additional_n100=(f'<p>{len(extra100)} additional original local records were independently refined and checked at their reported orders. They provide numerical states at N=65,70,72,75,84,96, including higher observed orders at N=65,70,72,75. These are numerical observations, not exact existence or maximum-order proofs. <a href="extra-source-audit-100.json">Read the source audit</a>.</p>' if extra100 else '')
extra100_downloads=('<a href="extra-source-audit-100.json">Additional N=61–100 source checks ↓</a><a href="extra-states-100.json">Additional numerical components ↓</a><a href="import_extra100.py">Additional-source importer ↓</a>' if extra100_path.exists() else '')
page=f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>ACstates library · Verified import through N={max_n}</title><link rel="stylesheet" href="../../styles.css"><style>main{{max-width:1000px;padding-block:40px}}h1{{font-size:clamp(32px,5vw,48px)}}h2{{margin-block:28px 14px}}.scroll{{overflow-x:auto}}table{{width:100%;border-collapse:collapse}}th,td{{padding:9px 16px;text-align:left;border-bottom:1px solid var(--line)}}td:last-child{{font:12px monospace;max-width:500px;overflow-wrap:anywhere}}.downloads{{display:flex;flex-wrap:wrap;gap:18px;margin-block:25px}}</style></head><body><main><a href="../../index.html">← AC-state catalogue</a><p class="eyebrow">Original library · Import and verification</p><h1>Every qubit number through N={max_n}</h1><p>The catalogue includes every row from N=1 to N={max_n}, with orders 1–{max_order}. Numerical library states are shown in blue; exact constructions and proved exclusions remain distinct. An unfilled higher-order cell records missing evidence rather than a failed-search exclusion.</p><h2>Numerical library representatives</h2><p>{len(states)} saved library representatives through N={max_n} were imported from the audited refinement of the original ACstates folder and additional source checks. Their reported orders were independently checked from saved high-precision Dicke components using {precision}-digit partial traces. The largest Frobenius residual through a reported order is {display(through_max)}{next_text}. These are numerical checks, not exact existence or maximum-order proofs.</p><p>The original coordinate and Dicke files disagree in some records. The catalogue follows the audited coordinate-derived representatives and their refinements, rather than substituting an inconsistent raw Dicke file. Original files and the refinement audit remain linked in the <a href="../../../majorana-degeneracy-catalogue/index.html">Majorana-degeneracy catalogue</a>.</p><p>The unlisted candidates acs.29.6_*.txt and acs.13.3.4.txt fail their claimed orders and are not counted as witnesses. The file named acs.48.8.txt contains 50 Majorana coordinates and fails its claimed eighth order; it is not used as an eighth-order N=48 witness. Additional accepted and rejected source checks are recorded in the audit downloads.</p>{additional_n60}{additional_n100}<h2>Exact lower-order examples</h2><p>{len(exact['states'])} additional analytical states were constructed to match the selected lower order exactly. Normalization, all required reductions and failure at the next order are checked with exact symbolic arithmetic. Together with existing representatives, these provide exact-order defaults whenever available. Where no matching-order construction is recorded, a higher-order representative is labelled by its own order.</p><div class="downloads"><a href="states.json">Saved numerical components ↓</a><a href="verification.json">{precision}-digit verification report ↓</a><a href="verify_numeric.py">Independent numerical verifier ↓</a><a href="exact-lower-states.json">New exact examples ↓</a><a href="verify_exact.py">Independent exact verifier ↓</a><a href="rejected-n29-t6.json">N=29 rejected candidate ↓</a><a href="rejected-n48-t8.json">N=48 rejected file ↓</a><a href="catalogue-audit.json">Catalogue consistency audit ↓</a><a href="additional-source-audit.json">Additional source checks ↓</a><a href="extra-state.json">N=37 refinement audit ↓</a><a href="extra-states-60.json">N=51 and N=60 refinement audits ↓</a>{extra100_downloads}</div><h2>Highest order observed in the imported library</h2><p>This column reports the highest audited numerical order available at each size. It does not claim a proved maximum.</p><div class="scroll"><table><thead><tr><th scope="col">N</th><th scope="col">Reported order</th><th scope="col">Largest residual through that order</th></tr></thead><tbody>{''.join(rows)}</tbody></table></div></main></body></html>'''
(HERE/'index.html').write_text(page)
