---
name: investigate-regression-mismatch
description: "Investigate make test failures or C vs C++ output mismatches in NeEstimator. Use when regression output differs, to classify timestamp-only noise vs scientific result changes and identify smallest fix."
argument-hint: "Optional failing command or file name"
user-invocable: true
---

# Investigate Regression Mismatch

Use this workflow when C and C++ outputs differ in regression testing.

## Standard Procedure

1. Reproduce the failure with make test.
2. If needed, run bash regression_test.sh to inspect step-level failures.
3. Compare files under .regression/baseline and .regression/current:
   - c.stdout vs cpp.stdout
   - c.test_Ne.txt vs cpp.test_Ne.txt
   - c.test_NexLD.txt vs cpp.test_NexLD.txt (when present)
4. Normalize runtime-only lines before comparing:
   - Starting time:
   - Ending time:
   - Time:
5. Classify mismatch type:
   - timestamp-only difference
   - missing/extra file
   - scientific output difference
6. Map the mismatch to the likely code area:
   - CLI dispatch and flow in Ne2x.c
   - C API boundary in Ne2x_api.h
   - C++ forwarding entrypoint in main.cpp
   - comparison behavior in regression_test.sh
7. Propose the smallest behavior-safe fix and list exact files touched.
8. Re-run make test and report final status.

## Output Format

- Reproduction result
- Mismatch classification
- Minimal fix plan
- Validation result

## Guardrails

- Do not modify Java GUI code under NeEstimator2x unless explicitly requested.
- Preserve C and C++ behavior compatibility unless the task explicitly asks for a behavior change.
- Prefer minimal, localized edits.
