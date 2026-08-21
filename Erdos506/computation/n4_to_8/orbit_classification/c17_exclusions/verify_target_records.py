from itertools import combinations
import importlib.util
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("c17_types", HERE / "c17_types.py")
targets = importlib.util.module_from_spec(spec)
spec.loader.exec_module(targets)


def main():
    report = {}
    all_triples = {tuple(c) for c in combinations(range(8), 3)}
    for class_id in sorted(targets.CLASSES):
        rec = targets.CLASSES[class_id]
        seen = {}
        duplicate = []
        for kind, blocks in [("L", rec["lines"]), ("C", rec["circle3"] + rec["circle4"])]:
            for block in blocks:
                ids = targets.parse_block(block)
                for tri in combinations(ids, 3):
                    if tri in seen:
                        duplicate.append([tri, seen[tri], kind + block])
                    seen[tri] = kind + block
        missing = sorted(all_triples - set(seen))
        report[class_id] = {
            "line_blocks": rec["lines"],
            "circle3": rec["circle3"],
            "circle4": rec["circle4"],
            "covered_triples": len(seen),
            "duplicate_cover_count": len(duplicate),
            "missing_triples": ["".join(map(str, t)) for t in missing],
            "ok": len(seen) == 56 and not duplicate and not missing,
        }
    print(json.dumps(report, indent=2, sort_keys=True))
    if not all(v["ok"] for v in report.values()):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
