"""Embed the existing local math/text fonts for standalone SVG downloads."""
from pathlib import Path
import base64
import json
import re

ROOT = Path(__file__).resolve().parents[1]
css = (ROOT / "vendor/katex/katex.min.css").read_text()


def woff(match):
    data = (ROOT / "vendor/katex" / match.group(1)).read_bytes()
    return 'src:url(data:font/woff2;base64,' + base64.b64encode(data).decode() + ') format("woff2")'


css = re.sub(r"src:url\((fonts/[^)]+\.woff2)\)[^}]*", woff, css)
assert "url(fonts/" not in css
for filename, style, weight in [
    ("lmroman10-regular.otf", "normal", 400),
    ("lmroman10-italic.otf", "italic", 400),
    ("lmroman10-bold.otf", "normal", 700),
]:
    data = base64.b64encode((ROOT / "vendor/fonts" / filename).read_bytes()).decode()
    css += ('@font-face{font-family:"Latin Modern Roman";font-style:' + style
            + ';font-weight:' + str(weight) + ';src:url(data:font/otf;base64,'
            + data + ') format("opentype")}')

(ROOT / "data/export-fonts.js").write_text(
    "/* Local font licenses are preserved under vendor/. */\n"
    + "globalThis.ATLAS_EXPORT_FONT_CSS=" + json.dumps(css, ensure_ascii=False) + ";\n"
)
print("Embedded local fonts for portable SVG labels.")
