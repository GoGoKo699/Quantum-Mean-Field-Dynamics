#!/usr/bin/env python3
"""Finite checks of the constructive precision and workspace lemmas.

Uses only Python's standard library.  The checks below are deliberately about
the nontrivial constructions, not a simulated quantum speedup:

* exact-rational Gaussian truncation/quantization parameter inequalities;
* dyadic endpoint signs, masks, and partial coarse-increment variances;
* exhaustive small Toeplitz-hash pairwise independence;
* indexed HashPRG blocks versus the recursive output definition;
* clean reversible bridges for noninvertible state updates, including arbitrary
  target registers, inverse schedules, exact leaf counts, and tape movement;
* explicit PRG error and quantum-mean-estimation numerical budgets.

The generic mathematical guarantees are proved in the accompanying chapters.
No amplitude-estimation experiment or fault-tolerant gate benchmark is run.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import json
import math
from pathlib import Path
import random


def power_of_two(exponent: int) -> Fraction:
    return Fraction(1 << exponent) if exponent >= 0 else Fraction(1, 1 << -exponent)


def ceil_log2(value: Fraction) -> int:
    """Exact ceiling of log2 for a positive rational."""
    if value <= 0:
        raise ValueError("logarithm argument must be positive")
    exponent = value.numerator.bit_length() - value.denominator.bit_length()
    if power_of_two(exponent) < value:
        exponent += 1
    assert power_of_two(exponent - 1) < value <= power_of_two(exponent)
    return exponent


def gaussian_parameters(eta: Fraction) -> tuple[int, int]:
    if not 0 < eta < 1:
        raise ValueError("eta must lie in (0,1)")
    cutoff = 1
    while 4 * power_of_two(-cutoff) * (cutoff * cutoff + 2 * cutoff + 2) > (eta / 3) ** 4:
        cutoff += 1
    uniform_bits = max(
        cutoff + 1 + ceil_log2((6 / eta) ** 2),
        ceil_log2(96 * (cutoff + 1) / eta),
    )
    return cutoff, uniform_bits


def check_gaussian_budgets() -> list[dict[str, int | float]]:
    records = []
    for accuracy_bits in (2, 4, 8, 16, 32, 64):
        eta = power_of_two(-accuracy_bits)
        cutoff, uniform_bits = gaussian_parameters(eta)
        tail_fourth = 4 * power_of_two(-cutoff) * (cutoff * cutoff + 2 * cutoff + 2)
        radial_squared = power_of_two(1 + cutoff - uniform_bits)
        angular_squared = 128 * cutoff * power_of_two(-2 * uniform_bits)
        assert tail_fourth <= (eta / 3) ** 4
        assert radial_squared <= (eta / 6) ** 2
        assert angular_squared <= (eta / 6) ** 2
        records.append(
            {
                "requested_error_bits": accuracy_bits,
                "normal_cutoff_parameter_s": cutoff,
                "bits_per_uniform": uniform_bits,
                "random_bits_per_normal": 2 * uniform_bits,
                "fourth_tail_budget_ratio": float(tail_fourth / (eta / 3) ** 4),
                "squared_radial_budget_ratio": float(radial_squared / (eta / 6) ** 2),
                "squared_angular_budget_ratio": float(angular_squared / (eta / 6) ** 2),
            }
        )
    return records


def clipped_duration(endpoint: Fraction, left: Fraction, width: Fraction) -> Fraction:
    return min(width, max(Fraction(0), endpoint - left))


def check_endpoint_alignment() -> dict[str, int]:
    endpoints_checked = 0
    mask_checks = 0
    variance_checks = 0
    for temporal_level in range(4):
        width = power_of_two(-temporal_level)
        for interval_index in range(1 << temporal_level):
            left = interval_index * width
            midpoint = left + width / 2
            for finest_level in range(6):
                endpoint_bits = max(1, finest_level - temporal_level + 1)
                bin_width = width * power_of_two(-endpoint_bits)
                assert temporal_level + endpoint_bits >= finest_level
                for cell in range(1 << endpoint_bits):
                    rounded = left + (Fraction(cell) + Fraction(1, 2)) * bin_width
                    for within_bin in (Fraction(1, 4), Fraction(3, 4)):
                        continuous = left + (cell + within_bin) * bin_width
                        assert (continuous < midpoint) == (rounded < midpoint)
                        for path_level in range(finest_level + 1):
                            step = power_of_two(-path_level)
                            assert continuous // step == rounded // step
                            for index in range(1 << path_level):
                                grid_left = index * step
                                assert (continuous > grid_left) == (rounded > grid_left)
                                mask_checks += 1
                        fine_step = power_of_two(-finest_level)
                        durations = [
                            clipped_duration(rounded, index * fine_step, fine_step)
                            for index in range(1 << finest_level)
                        ]
                        assert sum(durations, Fraction(0)) == rounded
                        if finest_level:
                            for coarse_index in range(1 << (finest_level - 1)):
                                expected = clipped_duration(
                                    rounded, 2 * coarse_index * fine_step, 2 * fine_step
                                )
                                assert durations[2 * coarse_index] + durations[2 * coarse_index + 1] == expected
                                variance_checks += 1
                        endpoints_checked += 1
    return {
        "rational_endpoint_pairs": endpoints_checked,
        "active_mask_comparisons": mask_checks,
        "partial_coarse_variance_identities": variance_checks,
    }


def toeplitz_hash(value: int, diagonal_seed: int, offset: int, block_bits: int) -> int:
    result = offset
    for row in range(block_bits):
        bit = 0
        for column in range(block_bits):
            matrix_bit = (diagonal_seed >> (row - column + block_bits - 1)) & 1
            bit ^= matrix_bit & ((value >> column) & 1)
        result ^= bit << row
    return result


def check_pairwise_independence() -> list[dict[str, int]]:
    records = []
    for block_bits in (1, 2, 3):
        domain = range(1 << block_bits)
        expected_count = 1 << (block_bits - 1)
        pairs = 0
        for first in domain:
            for second in domain:
                if first == second:
                    continue
                counts: Counter[tuple[int, int]] = Counter()
                for diagonal_seed in range(1 << (2 * block_bits - 1)):
                    for offset in domain:
                        counts[
                            (
                                toeplitz_hash(first, diagonal_seed, offset, block_bits),
                                toeplitz_hash(second, diagonal_seed, offset, block_bits),
                            )
                        ] += 1
                assert len(counts) == 1 << (2 * block_bits)
                assert set(counts.values()) == {expected_count}
                pairs += 1
        records.append({"block_bits": block_bits, "ordered_distinct_input_pairs": pairs})
    return records


HashParameters = tuple[int, int]


def indexed_hashprg_block(
    start: int, hashes: list[tuple[HashParameters, HashParameters]], index: int, block_bits: int
) -> int:
    value = start
    for depth in reversed(range(len(hashes))):
        diagonal_seed, offset = hashes[depth][(index >> depth) & 1]
        value = toeplitz_hash(value, diagonal_seed, offset, block_bits)
    return value


def recursive_hashprg(
    start: int, hashes: list[tuple[HashParameters, HashParameters]], block_bits: int
) -> list[int]:
    if not hashes:
        return [start]
    low_hash, high_hash = hashes[-1]
    low_start = toeplitz_hash(start, *low_hash, block_bits)
    high_start = toeplitz_hash(start, *high_hash, block_bits)
    return recursive_hashprg(low_start, hashes[:-1], block_bits) + recursive_hashprg(
        high_start, hashes[:-1], block_bits
    )


def check_block_addressing() -> dict[str, int]:
    rng = random.Random(30017)
    blocks = 0
    configurations = 0
    for block_bits in (1, 2, 3, 5):
        for depth in range(1, 7):
            for _ in range(8):
                start = rng.randrange(1 << block_bits)
                hashes = [
                    tuple(
                        (rng.randrange(1 << (2 * block_bits - 1)), rng.randrange(1 << block_bits))
                        for _ in range(2)
                    )
                    for _ in range(depth)
                ]
                recursive = recursive_hashprg(start, hashes, block_bits)
                indexed = [
                    indexed_hashprg_block(start, hashes, index, block_bits)
                    for index in range(1 << depth)
                ]
                assert recursive == indexed
                configurations += 1
                blocks += len(indexed)
    return {"seed_configurations": configurations, "addressed_blocks": blocks}


# A leaf operation means target ^= f_step(source; seed).  All local Boolean
# scratch has already been uncomputed by the leaf bridge construction.
Leaf = tuple[int, int, int]


def checkpoint_schedule(
    start: int,
    depth: int,
    arity: int,
    source: int = 0,
    target: int = 1,
    scratch_base: int = 2,
) -> list[Leaf]:
    if depth == 0:
        return [(start, source, target)]
    child_length = arity ** (depth - 1)
    intermediate = list(range(scratch_base, scratch_base + arity - 1))
    child_scratch = scratch_base + arity - 1
    boundaries = [source] + intermediate + [target]
    children = [
        checkpoint_schedule(
            start + index * child_length,
            depth - 1,
            arity,
            boundaries[index],
            boundaries[index + 1],
            child_scratch,
        )
        for index in range(arity)
    ]
    forward = [leaf for child in children for leaf in child]
    cleanup = [leaf for child in reversed(children[:-1]) for leaf in reversed(child)]
    return forward + cleanup


def step_function(index: int, value: int, seed: int, state_mask: int) -> int:
    """Intentionally noninjective, seed-dependent arithmetic update."""
    return (2 * value + (value >> 1) + (index + 1) * (seed + 3)) & state_mask


def execute_schedule(schedule: list[Leaf], registers: list[int], seed: int, state_mask: int) -> None:
    for index, source, target in schedule:
        assert source != target
        registers[target] ^= step_function(index, registers[source], seed, state_mask)


def direct_composition(length: int, value: int, seed: int, state_mask: int) -> int:
    for index in range(length):
        value = step_function(index, value, seed, state_mask)
    return value


def check_checkpoint_bridges() -> list[dict[str, int | float]]:
    records = []
    state_bits = 3
    state_mask = (1 << state_bits) - 1
    for arity, maximum_depth in ((2, 5), (3, 4), (4, 3)):
        for depth in range(maximum_depth + 1):
            schedule = checkpoint_schedule(0, depth, arity)
            length = arity**depth
            leaf_count = (2 * arity - 1) ** depth
            register_count = 2 + (arity - 1) * depth
            assert len(schedule) == leaf_count
            assert max(max(source, target) for _, source, target in schedule) < register_count
            indices = [index for index, _, _ in schedule]
            assert indices[0] == indices[-1] == 0
            travel = sum(abs(second - first) for first, second in zip(indices, indices[1:]))
            assert travel == 2 * (leaf_count - length)
            basis_cases = 0
            for seed in range(4):
                for source_value in range(1 << state_bits):
                    expected = direct_composition(length, source_value, seed, state_mask)
                    for target_value in range(1 << state_bits):
                        registers = [source_value, target_value] + [0] * (register_count - 2)
                        execute_schedule(schedule, registers, seed, state_mask)
                        assert registers == [source_value, target_value ^ expected] + [0] * (register_count - 2)
                        execute_schedule(list(reversed(schedule)), registers, seed, state_mask)
                        assert registers == [source_value, target_value] + [0] * (register_count - 2)
                        basis_cases += 1
            records.append(
                {
                    "arity": arity,
                    "recursion_depth": depth,
                    "irreversible_steps": length,
                    "leaf_bridges": leaf_count,
                    "configuration_registers": register_count,
                    "classical_tape_row_travel": travel,
                    "basis_cases_and_inverse_checks": basis_cases,
                    "gate_exponent": math.log(2 * arity - 1, arity),
                }
            )
    return records


def check_prg_error_budgets() -> dict[str, int | float]:
    cases = 0
    maximum_bad_hash_ratio = Fraction(0)
    for state_width in (1, 8, 31, 128):
        for depth in (1, 2, 5, 12):
            for chi_bits in (2, 5, 16, 64):
                chi = power_of_two(-chi_bits)
                block_bits = (
                    3 * state_width
                    + 2 * depth
                    + ceil_log2(Fraction(16 * (depth + 1)))
                    + 3 * ceil_log2(1 / chi)
                    + 2
                )
                theta = chi / 2
                bad_hash_bound = (
                    Fraction(((1 << depth) - 1) ** 2 * depth)
                    * power_of_two(3 * state_width - block_bits)
                    / theta**2
                )
                assert bad_hash_bound <= chi / 16
                assert theta / 2 + bad_hash_bound < chi
                maximum_bad_hash_ratio = max(maximum_bad_hash_ratio, bad_hash_bound / chi)
                cases += 1
    # Include increasing worst-case finite ranges: a pointwise small Hh bound is
    # deliberately not assumed in the first/second-moment transfer.
    for exponent in range(1, 17):
        rho = power_of_two(-2 * exponent)
        sigma = power_of_two(-exponent)
        finite_range = power_of_two(exponent)
        threshold = min(Fraction(1, 4), rho / (2 * finite_range), sigma**2 / finite_range**2)
        chi = power_of_two(-ceil_log2(1 / threshold))
        assert 2 * finite_range * chi <= rho
        assert finite_range**2 * chi <= sigma**2
    return {
        "exact_parameter_cases": cases,
        "largest_bad_hash_bound_divided_by_chi": float(maximum_bad_hash_ratio),
        "polynomial_range_moment_transfer_cases": 16,
    }


def check_mean_estimation_budgets() -> list[dict[str, int | float]]:
    # Use the rational upper bound pi < 22/7 and sqrt(2) < 3/2.  Thus these
    # inequalities do not depend on binary floating-point evaluations of pi.
    pi_upper = Fraction(22, 7)
    sqrt_two_upper = Fraction(3, 2)
    records = []
    for ratio in (Fraction(101, 100), Fraction(2), Fraction(7), Fraction(64), Fraction(10000)):
        sigma = Fraction(1)
        error = 1 / ratio
        bins = ceil_log2(8 * ratio)
        query_parameter = 1 << ceil_log2(128 * pi_upper * ratio * (bins + 1))
        first_sum = 4 * pi_upper * sqrt_two_upper * sigma * (bins + 1) / query_parameter
        second_sum = 4 * pi_upper**2 * sigma * (1 << bins) / query_parameter**2
        tail = sigma / (1 << bins)
        assert first_sum + second_sum < error / 8
        assert tail <= error / 8
        assert first_sum + second_sum + tail < error / 4
        assert first_sum + second_sum + tail + error / 4 < error
        records.append(
            {
                "sigma_over_error": float(ratio),
                "largest_bin_index": bins,
                "power_of_two_query_parameter": query_parameter,
                "amplitude_error_bound_over_error": float((first_sum + second_sum) / error),
                "tail_bound_over_error": float(tail / error),
            }
        )
    # Coefficient ledger includes optional PRG and final coefficient rounding.
    assert Fraction(1, 4) + Fraction(1, 16) + Fraction(1, 4) + Fraction(1, 16) + Fraction(1, 16) == Fraction(11, 16)
    return records


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).with_name("resource_construction_check.json"),
    )
    args = parser.parse_args()
    result = {
        "scope": "Finite construction checks, not quantum execution or a hardware advantage benchmark",
        "all_checks_passed": True,
        "gaussian_budget_checks": check_gaussian_budgets(),
        "endpoint_alignment_checks": check_endpoint_alignment(),
        "toeplitz_pairwise_independence_checks": check_pairwise_independence(),
        "hashprg_block_addressing_checks": check_block_addressing(),
        "checkpoint_bridge_checks": check_checkpoint_bridges(),
        "prg_budget_checks": check_prg_error_budgets(),
        "quantum_mean_budget_checks": check_mean_estimation_budgets(),
        "coefficient_error_budget_fraction": "11/16",
        "ternary_checkpoint_exponent": math.log(5, 3),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"all_checks_passed": True, "output": str(args.output)}, indent=2))


if __name__ == "__main__":
    main()
