---
description: "Use when editing C/C++ core memory handling in NeEstimator (malloc/free, ownership, early returns, ASan findings)."
name: "C Core Memory Safety"
applyTo:
  - "**/*.c"
  - "**/*.h"
  - "**/*.cpp"
---

# C/C++ Core Memory Safety

- Preserve scientific behavior and keep C/C++ entrypoints aligned.
- Keep edits localized, especially in Ne2x.c.
- For every allocation, identify owner and exactly one release path.
- When a function has multiple early returns, prefer a shared cleanup section to avoid leaks.
- Initialize pointers to NULL and lengths/counters to safe defaults before first use.
- Guard string and buffer writes with explicit bounds.
- Keep CLI argument contract stable: i:<info_file>, optional o:<option_file>, m:/m+:, c:.

## Required Validation After Memory-Sensitive Changes

1. Run make
2. Run make test
3. Run make asan
4. Run ./Ne2x_asan i:test_info.txt o:test_option.txt
