#!/usr/bin/env python3
"""Pure generator of the point-pair profiles used by the isolated model."""

from __future__ import annotations

from functools import lru_cache


@lru_cache(maxsize=None)
def profiles(n: int, maximum_block_size: int) -> tuple[tuple[int, ...], ...]:
    """Return all exact off-pair partitions; no solver or file is involved."""
    sizes = tuple(range(3, maximum_block_size + 1))
    answer: set[tuple[int, ...]] = set()
    for line_size in (None, *sizes):
        line = tuple(int(line_size == k) for k in sizes)
        used = 0 if line_size is None else line_size - 2
        work = [0] * len(sizes)

        def rec(index: int, remaining: int) -> None:
            if index == len(sizes):
                if remaining == 0:
                    answer.add(line + tuple(work))
                return
            weight = sizes[index] - 2
            for number in range(remaining // weight + 1):
                work[index] = number
                rec(index + 1, remaining - number * weight)
            work[index] = 0

        rec(0, n - 2 - used)
    return tuple(sorted(answer))
