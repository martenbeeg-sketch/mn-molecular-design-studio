"""Scientific validity is independent of process completion and acceptance."""
from __future__ import annotations

import math

# Physical output domains, not task acceptance thresholds.
RANGES = {"ipsae": (0, 1), "ipsae_min": (0, 1), "ipsae_max": (0, 1), "ipsae_avg": (0, 1),
          "iptm": (0, 1), "ptm": (0, 1), "binder_plddt": (0, 100), "ipae": (0, None), "ipde": (0, None),
          "binder_rmsd": (0, None), "monomer_rmsd": (0, None), "qc_sequence_length": (1, None), "qc_sequence_pi": (0, 14)}


def finite_number(value, *, allow_bool=False):
    if value is None or (isinstance(value, bool) and not allow_bool):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError):
        return None
    return number if math.isfinite(number) else None


def measurement_outcome(metrics, required, error="", *, boolean_fields=()):
    def valid(key):
        number = finite_number(metrics[key], allow_bool=key in boolean_fields)
        low, high = RANGES.get(key, (None, None))
        return number is not None and (low is None or number >= low) and (high is None or number <= high)
    fields = {key: "missing" if metrics.get(key) is None else
              "valid" if valid(key) else "invalid"
              for key in required}
    invalid = [f"{key}: {state}" for key, state in fields.items() if state != "valid"]
    return {"status": "error" if error else "invalid" if "invalid" in fields.values() else
            "missing" if invalid else "measured", "fields": fields,
            "error": str(error or ("Required measurements unavailable: " + "; ".join(invalid) if invalid else ""))}


def measurement_summary(outcomes):
    valid = sum(item["status"] == "measured" for item in outcomes)
    return {"candidate_count": len(outcomes), "scored_count": valid,
            "error_count": sum(item["status"] == "error" for item in outcomes),
            "missing_count": sum(item["status"] == "missing" for item in outcomes),
            "invalid_count": sum(item["status"] == "invalid" for item in outcomes),
            "scientific_status": "failed" if not valid else "valid" if valid == len(outcomes) else "partial"}
