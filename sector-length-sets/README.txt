SECTOR-LENGTH GEOMETRY — PURE SYMMETRIC STATES, N=3–6

The site needs no build step, accounts, API keys, or external network services.
A browser supporting WebGL is required. Local HTTP preview has been tested.
Opening index.html directly depends on your browser allowing local script/font
assets; GitHub Pages serves these static files normally. For local preview:
  python3 -m http.server 8766
then visit http://localhost:8766/ in a browser.

Publish on GitHub Pages
-----------------------
Copy the COMPLETE contents of this folder into, for example,
  symmetric-sector-catalogue/sector-length-sets/
in the repository already serving the catalogue. The new page is then at
  .../symmetric-sector-catalogue/sector-length-sets/index.html
All local asset paths are relative, so a different subdirectory also works.
No npm installation, compilation, API keys, or deployment service is needed.
Keep vendor/, data/, references/, and .nojekyll with the HTML files.

Suggested link from the existing catalogue:
  <a href="sector-length-sets/index.html">Sector-length sets, N=3–6</a>

What the page shows
-------------------
N=3: exact interval, named-state nodes.
N=4: exact triangle, named-state nodes; GHZ and balanced Dicke remain separate
     choices at their common sector point.
N=5: exact proved range, six named boundary states and 41 sampled boundary
     states, each with an explicit analytical expression and exact parameter.
N=6: current proved necessary outer region, constructive violet patches, and
     the manuscript's 27 nodes: 16 violet positions and 11 analytically excluded
     gray positions, plus 31 explicit attained boundary samples. The shared
     psi6/chi6 position has two state choices. No interior cloud is shown.

Hover a violet point to preview its state and Majorana sphere. Click to select
it; the selection is restored when leaving other previews. The typeset state
selector provides keyboard/touch access, including sampled boundary states.
Use arrows or type a name to navigate; Enter selects and Escape closes.
Dragging either three-dimensional view rotates it. The six-qubit view also has
main/alternative/rear presets; Majorana rotation has keyboard-accessible buttons.
The gray-node readout explains exclusion and never invents a state or sphere.

Known analytical state expressions appear before numerical coefficients.
Open 'Exact parameters & definitions' for longer algebraic definitions.
All mathematical labels use $...$ notation and local KaTeX typesetting.
Text uses the LaTeX Latin Modern Roman fonts, also included locally.

The original manuscripts and figures in the research workspace are unchanged.
The references/ directory contains copies of the supporting N5 and N6 PDFs.

Scientific scope
----------------
The N=6 body is an OUTER REGION, not the full exact achievable set. Gray surfaces
can be excluded or uncertified; every displayed gray NODE is analytically
excluded. The refined four-copy lower cap is excluded except at the coherent
product point. The display uses refined analytical surface meshes, with bisection-clipped
constructive color transitions and smooth normals within each sheet. True
creases stay distinct. Mesh interpolation remains a display approximation,
not a new exact claim for all interpolated triangle interiors. No convex hull of states is declared
achievable. See notes.html and the provenance files for full definitions.

Data and reproduction
----------------------
Normalized Dicke amplitudes are [real,imaginary] pairs, ordered k=0,...,N.
Sector lengths use the sum-of-squared-Pauli-correlations convention, S0=1.
The browser independently recomputes all displayed sectors and reconstructs
each state from its computed Majorana stars to check the fidelity.

data/small-n.json and data/n6.json retain the full research witness datasets,
including witnesses intentionally NOT plotted. The smaller *-display.js files
are the same geometry and selected records, without the hidden sample clouds.
scripts/build_display_data.py regenerates these lightweight browser datasets.

The math and state audits run with Node, with no third-party modules:
  node test_majorana_math.cjs
  node audit-atlas.cjs
  node audit-visible.cjs

Boundary/mesh export scripts record their parent research-workspace inputs;
regenerating the N6 analytical geometry requires those original Julia libraries.
The delivered page and independent numerical audits require no Julia runtime.

Third-party licenses
--------------------
Three.js (MIT): vendor/THREE-LICENSE.txt
KaTeX (MIT): vendor/katex/LICENSE.txt
Latin Modern (GUST Font License): vendor/fonts/GUST-FONT-LICENSE.txt
All third-party resources are local. The scientific sources are linked in notes.

Prepared 6 November 2026. This folder has not been published or pushed to GitHub.

Marker presentation update
--------------------------
Majorana stars and set vertices use resolution-independent circular vector
markers with thin darker outlines, matching the catalogue style. Front and
rear points have distinct fills; rear set nodes remain lighter. The marker
positions, state catalogues, analytical expressions and hover targets are
unchanged. No extra sector-set nodes were introduced.

Hull, symmetry and navigation update
-----------------------------------
Constellation hulls draw only genuine edges, with no face-triangulation diagonals.
Rotational C_n axes are detected from the full multiplicity-weighted constellation.
Use the order selector for highest order, a selected order, or all axes. Continuous
axial cases show C_infinity; balanced Dicke states also have a continuous family
of transverse C_2 axes, described in the summary rather than discretized.

Drag to move the N=3,4,5 sets. For N=6 use Shift-drag or right-drag to move and
ordinary drag to rotate. Reset view restores the initial framing.

Both panels have Export SVG buttons. They preserve the current framing, displayed
axes, true hull edges, markers, labels and leaders. The exports contain vector
geometry with embedded LaTeX fonts; math labels use SVG foreignObject/XHTML, which
is supported in browsers. No screenshot raster is embedded. The state selector typesets both its selected value and its open choices,
including Greek letters, subscripts and superscripts. A hidden native select
retains the state-data interface and provides the no-enhancement fallback.

Align highest C_n with z rotates the Majorana display frame so a highest-order
axis is vertical. It also resets the camera for a clear view. The stored Dicke
coefficients retain their original basis; this is a collectively rotated view of
the same sector-length point. Reset restores the original display frame. SVG
export includes the aligned constellation and its axes.

Exact values and smooth surfaces
-------------------------------
All selectable named states and boundary samples have verified analytical
sector tuples. Rational values are shown as fractions and algebraic values as
radicals. Parameter-dependent tuples include expandable exact definitions.
Numerical Dicke coefficients and independent checks remain available below.
The exact tuples are in data/exact-sector-values.js and its JSON twin.

The N=6 skin uses 156,910 triangles in 17 analytical sheets and full-resolution
true-edge paths. Three conforming refinement passes concentrate triangles near
violet/gray transitions, preserving the same constructive tests. Attainment transitions are clipped by 20 bisection steps,
with new vertices projected onto their analytical graph. Normals are shared
within each sheet, while genuine creases remain. The original research
geometry and all 27 manuscript nodes are unchanged. New geometry is in
data/n6-smooth-geometry.js; scripts/export_smooth_n6.jl reproduces it.

Labels default to key extrema plus the selected or hovered state. Use Labels
to show all known-state labels and family callouts, or hide them. All points
remain selectable in every mode. Export SVG follows the chosen label mode.

Sector-plot zoom update
-----------------------
Sector sets allow up to 80x magnification. Scroll over the feature to enlarge
it around the pointer; drag to move the N=3–5 views, or Shift/right-drag for N=6.
Reset view restores the original framing. Majorana zoom keeps its own range.

Six-qubit boundary samples
--------------------------
Eight selectable states sample the filtered-octahedral C3 color seam from octa
to product. Seven sample the proved attained straight C5-junction-to-omega
equality segment. Each has normalized Dicke coefficients, analytical formulas,
exact sector values and its independently reconstructed Majorana constellation.
New points have compact violet markers and a label only when focused. They are
also available under Sampled boundary states. The attained C5-to-omega segment
and the filtered-octahedral C3 curve are traced in violet, with rear portions
lighter and dashed.
The nearby curved rank-four outer crease is unattainable except at its C5 end;
no artificial state is assigned to it. See notes.html for the distinction.

Long state definitions
-----------------------
Long exact parameter definitions are collapsed when entering a different state.
They remain expandable, and a manually opened panel stays open while the same
state is displayed.

Curve and display update
-------------------------
The filtered-octahedral C3 trace uses 721 points, including each displayed F_t
state exactly. Rear portions are lighter and dashed. Both front and rear
lines leave the face-depth buffer unchanged; rear dashes keep 6px strokes and
5px gaps across zoom levels, including SVG export.
The N=5 plot includes a brief nonconvexity note. Violet shaded swatches identify the
Attained region for N=4,5 and the Attained surface for N=6. Full evidence
remains in notes.

Paired-face boundary samples and parameter layout
-------------------------------------------------
The straight and elliptic borders of the attained paired-purity face have
additional explicit pure-state samples. The elliptic arc uses a three-Dicke
C2 family; the straight part uses the general real-Hodge construction.
Long parameter lists keep the selected parameter on its own line above the
definitions, so a fraction cannot drift beside the tall aligned block.
