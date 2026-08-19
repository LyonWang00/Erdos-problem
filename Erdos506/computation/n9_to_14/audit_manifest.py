#!/usr/bin/env python3
"""Audit terminal certificates and write their source manifest."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
CERTIFICATES = HERE / "certificates"
EXPECTED_N10 = [10, 0, 0, 0, 10, 20, 2, 0, 15]


def read(name: str) -> dict[str, object]:
    return json.loads((CERTIFICATES / name).read_text(encoding="utf-8"))


def main() -> None:
    largest_block = read("largest_block_reduction.json")
    n9 = read("n9_04_exact_support.json")
    n10 = read("n10_07_recursive.json")
    n10_geometry = read("n10_08_final_geometry.json")
    n11 = read("n11_06_recursive.json")
    n12 = read("n12_06_recursive.json")
    n13_n14 = read("n13_n14_global_inequalities.json")
    retained_n10 = [
        record["signature"] for record in n10["levels"]["10"]
        if not str(record["status"]).startswith("INFEASIBLE")
    ]
    checks = {
        "largest_block_reduction": largest_block["status"] == "PASS"
            and [record["n"] for record in largest_block["records"]]
            == list(range(7, 16)),
        "n9_closed": n9["status"] == "PASS"
            and n9["retained_signature_count"] == 0,
        "n10_unique_residual": retained_n10 == [EXPECTED_N10],
        "n10_geometry_closed": n10_geometry["status"] == "PASS",
        "n11_closed": n11["status"] == "PASS"
            and n11["retained_signature_count"] == 0,
        "n12_closed": n12["status"] == "PASS"
            and n12["retained_signature_count"] == 0,
        "n13_n14_closed": n13_n14["status"] == "PASS"
            and all(record["retained_signature_count"] == 0
                    for record in n13_n14["records"]),
    }
    files = sorted(HERE.glob("*.py")) + sorted(
        path for path in CERTIFICATES.glob("*.json")
        if path.name != "manifest.json"
    )
    document = {
        "schema_version": 1,
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "files": {
            str(path.relative_to(HERE)).replace("\\", "/"): {
                "bytes": path.stat().st_size,
                "sha256": sha256(path.read_bytes()).hexdigest(),
            }
            for path in files
        },
    }
    output = CERTIFICATES / "manifest.json"
    output.write_text(json.dumps(document, indent=2) + "\n",
                      encoding="utf-8")
    print(json.dumps({
        "status": document["status"],
        "checks": checks,
        "file_count": len(files),
        "output": str(output),
    }, indent=2))
    if document["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
