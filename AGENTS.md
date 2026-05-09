# AGENTS.md

## Scope

- Primary focus is the C/C++ core in this repository.
- Do not modify the Java GUI under NeEstimator2x unless explicitly requested.
- Preserve behavior compatibility between C and C++ entrypoints unless the task explicitly asks for a behavior change.

## Quick Start Commands

Run from repository root.

- Build all core artifacts:
  - make
- Run C vs C++ regression:
  - make test
- Build AddressSanitizer binary:
  - make asan
- Run a representative ASan scenario:
  - ./Ne2x_asan i:test_info.txt o:test_option.txt

## Architecture Map (C/C++)

- Core monolithic implementation and CLI dispatch:
  - [Ne2x.c](Ne2x.c)
- C API boundary used by C++ wrapper:
  - [Ne2x_api.h](Ne2x_api.h)
- Thin C++ forwarding entrypoint:
  - [main.cpp](main.cpp)
- Build orchestration and standard targets:
  - [Makefile](Makefile)
- Regression procedure details and output comparison behavior:
  - [regression_test.sh](regression_test.sh)
- Human-facing project overview and usage:
  - [README.md](README.md)
- Debug workflow details:
  - [DEBUG_GUIDE.md](DEBUG_GUIDE.md)

## Behavioral Constraints for Code Changes

- Keep the CLI argument contract stable:
  - i:<info_file>, optional o:<option_file>, m:/m+: and c: flows.
- Keep C and C++ execution paths functionally aligned:
  - C binary and C++ wrapper must continue to produce matching outputs in regression tests.
- Prefer minimal, localized edits in Ne2x.c and avoid broad refactors unless requested.
- Treat AddressSanitizer findings as first-class failures for touched code paths.

## Validation Checklist After C/C++ Changes

- Run make
- Run make test
- If memory-sensitive code is touched, run:
  - make asan
  - ./Ne2x_asan i:test_info.txt o:test_option.txt

## Known Testing Detail

- Regression comparison intentionally ignores runtime timestamp lines (for example Starting time, Ending time, Time) to avoid false negatives from clock differences.
