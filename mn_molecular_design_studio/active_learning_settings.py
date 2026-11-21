"""Validated study settings, independent of UI, model runtime and scheduler."""
from __future__ import annotations

import math
from dataclasses import asdict, dataclass


def _integer(value, name, minimum):
    if type(value) is not int or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}.")


@dataclass(frozen=True)
class ActiveLearningSettings:
    binder_percent_start: float = 10.0
    binder_percent_target: float = 0.01
    tightening_rounds: int = 10
    max_iterations: int = 0

    def __post_init__(self):
        for name in ("binder_percent_start", "binder_percent_target"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or not 0 < value <= 100:
                raise ValueError(f"{name} must be a finite percentage in (0, 100].")
        if self.binder_percent_target > self.binder_percent_start:
            raise ValueError("The target binder percentile cannot exceed the starting percentile.")
        _integer(self.tightening_rounds, "tightening_rounds", 1)
        _integer(self.max_iterations, "max_iterations", 0)

    def percentile(self, round_number: int) -> float:
        _integer(round_number, "round_number", 1)
        if self.tightening_rounds == 1 or round_number >= self.tightening_rounds:
            return float(self.binder_percent_target)
        fraction = (round_number - 1) / (self.tightening_rounds - 1)
        return self.binder_percent_start * (self.binder_percent_target / self.binder_percent_start) ** fraction

    def iteration_limit_reached(self, completed_rounds: int) -> bool:
        _integer(completed_rounds, "completed_rounds", 0)
        return bool(self.max_iterations and completed_rounds >= self.max_iterations)

    def convergence_count(self, frozen_reference_count: int) -> int:
        """Call once at target percentile; the controller must persist the result."""
        _integer(frozen_reference_count, "frozen_reference_count", 1)
        return max(1, math.ceil(frozen_reference_count * self.binder_percent_target / 100))

    def to_dict(self):
        return asdict(self)


def next_round_selection_counts(*, initial_seed_count: int, remaining_compounds: int,
                                max_batch_size: int, seed_fraction_percent: float = 50.0,
                                predicted_fraction_percent: float = 80.0) -> dict[str, int]:
    """Size a later docking batch relative to the actual initial centroid sample."""
    if initial_seed_count < 1:
        raise ValueError("The initial centroid sample must contain at least one compound.")
    if not 1 <= float(seed_fraction_percent) <= 100:
        raise ValueError("The next-round fraction of the initial centroid sample must be between 1 and 100 percent.")
    if not 0 <= float(predicted_fraction_percent) <= 100:
        raise ValueError("The predicted-compound fraction must be between 0 and 100 percent.")
    if max_batch_size < 0:
        raise ValueError("The maximum docking batch size cannot be negative.")

    target_count = max(1, int(math.floor(initial_seed_count * float(seed_fraction_percent) / 100.0 + 0.5)))
    total = min(target_count, max(0, int(remaining_compounds)))
    if max_batch_size:
        total = min(total, int(max_batch_size))
    predicted = int(math.floor(total * float(predicted_fraction_percent) / 100.0 + 0.5))
    return {"total": total, "predicted": predicted, "random": total - predicted, "target_total": target_count}
