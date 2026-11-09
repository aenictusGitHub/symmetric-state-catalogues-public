"""Build SVG glyph data from the explorers' vendored fonts (requires fontTools)."""
import json
from pathlib import Path

from fontTools.pens.svgPathPen import SVGPathPen
from fontTools.ttLib import TTFont

root = Path(__file__).resolve().parents[1]
vendor = root / "sector-length-sets/vendor"
sources = [(p.stem, p) for p in sorted((vendor / "katex/fonts").glob("*.ttf"))]
sources += [
    ("Latin Modern Roman-Regular", vendor / "fonts/lmroman10-regular.otf"),
    ("Latin Modern Roman-Bold", vendor / "fonts/lmroman10-bold.otf"),
    ("Latin Modern Roman-Italic", vendor / "fonts/lmroman10-italic.otf"),
]
fonts = {}
for name, source in sources:
    with TTFont(source) as font:
        glyph_set = font.getGlyphSet()
        paths = {}
        for codepoint, glyph in sorted(font.getBestCmap().items()):
            pen = SVGPathPen(glyph_set)
            glyph_set[glyph].draw(pen)
            paths[chr(codepoint)] = pen.getCommands()
        fonts[name] = [font["head"].unitsPerEm, paths]
output = root / "sector-length-sets/data/label-outlines.js"
output.write_text(
    "/* Generated from vendored KaTeX and Latin Modern fonts. See vendor font licenses. */\n"
    "globalThis.ATLAS_LABEL_OUTLINES = " + json.dumps(fonts, separators=(",", ":")) + ";\n"
)
