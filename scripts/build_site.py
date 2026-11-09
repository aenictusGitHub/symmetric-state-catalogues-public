"""Build the public catalogues and sector-length atlas, retaining large certificate downloads on GitHub."""
from pathlib import Path
from urllib.parse import quote, unquote, urlsplit
import os
import re
import shutil

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / "_site"
CATALOGUES = ("symmetric-sector-catalogue", "anticoherent-state-catalogue", "ac-subspace-catalogue", "majorana-degeneracy-catalogue", "sector-length-sets", "reduced-purity-catalogue")
REPOSITORY = os.environ.get("GITHUB_REPOSITORY", "aenictusGitHub/symmetric-state-catalogues-public")
REVISION = os.environ.get("GITHUB_SHA", "main")
LIMIT = 10 * 1024**2

if DEST.exists():
    shutil.rmtree(DEST)
DEST.mkdir()
files = [ROOT / name for name in ("index.html", "styles.css", ".nojekyll")]
for name in CATALOGUES:
    files.extend(f for f in (ROOT / name).rglob("*") if f.is_file() and not any(x in f.parts for x in ("__pycache__", "node_modules", ".git")) and f.name != ".DS_Store")
# The analytical atlas mesh is a runtime script; keep it on Pages.
large = {f.resolve() for f in files if f.stat().st_size > LIMIT
         and not (f.is_relative_to(ROOT / "sector-length-sets") and f.suffix in {".js", ".css"})}

def rewrite_link(match, source):
    url = match.group(2)
    parsed = urlsplit(url)
    if parsed.scheme or parsed.netloc or not parsed.path:
        return match.group(0)
    target = (source.parent / unquote(parsed.path)).resolve()
    if target not in large:
        return match.group(0)
    remote = "https://raw.githubusercontent.com/" + REPOSITORY + "/" + REVISION + "/" + quote(target.relative_to(ROOT).as_posix())
    return match.group(1) + remote + match.group(3)

for source in files:
    if source.resolve() in large:
        continue
    target = DEST / source.relative_to(ROOT)
    target.parent.mkdir(parents=True, exist_ok=True)
    if source.suffix == ".html":
        html = re.sub(r'(\b(?:href|src)=["\'])([^"\']+)(["\'])', lambda m: rewrite_link(m, source), source.read_text())
        target.write_text(html)
    else:
        shutil.copy2(source, target)
size = sum(f.stat().st_size for f in DEST.rglob("*") if f.is_file())
assert size < 1024**3, "Published site exceeds GitHub Pages size limit"
assert not any("slocc-state-catalogue" in str(f.relative_to(DEST)) for f in DEST.rglob("*"))
print(f"Built catalogues and sector-length atlas: {size / 1024**2:.1f} MiB; {len(large)} large certificate files linked to GitHub.")
