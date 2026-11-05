#!/usr/bin/env python3
"""Regenerate the browser data and CSV from the canonical catalogue JSON."""
import csv
import json
from pathlib import Path

root = Path(__file__).resolve().parents[1]
data = json.loads((root / "data/catalogue.json").read_text())
keys = set()
sources = {s["id"] for s in data["sources"]}
for entry in data["entries"]:
    key = (entry["n"], entry["k"], entry["mode"])
    assert key not in keys, f"Duplicate result {key}"
    keys.add(key)
    assert 1 <= entry["k"] <= entry["n"]
    assert entry["mode"] in {"min", "max"}
    assert entry["status"] in {"proven", "reported", "numerical", "bound"}
    assert entry["references"] and all(ref["sourceId"] in sources for ref in entry["references"])
    if entry["status"] == "bound":
        assert entry.get("boundDirection") in {"lower", "upper"}
for source in data["sources"]:
    if source.get("pdf") or source.get("href"):
        local_path = source.get("pdf", source.get("href"))
        assert (root / local_path).is_file(), f"Missing source: {local_path}"
(root / "data/catalogue.js").write_text(
    "window.SECTOR_CATALOGUE = " + json.dumps(data, ensure_ascii=False, separators=(",", ":")) + ";\n"
)
fields = ["n", "k", "mode", "value", "decimal", "status", "stateName", "sourcePage", "sourceLocation", "notes"]
with (root / "data/catalogue.csv").open("w", newline="", encoding="utf-8") as stream:
    writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(data["entries"])
print(f"Validated {len(keys)} unique entries and {len(sources)} source records; regenerated JS and CSV.")
