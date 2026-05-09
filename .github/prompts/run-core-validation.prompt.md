---
description: "Run standard NeEstimator C/C++ core validation: build, regression parity, and ASan scenario."
name: "Run Core Validation"
argument-hint: "Optional extra validation scenario"
agent: "agent"
---

Validate the C/C++ core from repository root using this sequence:

1. Run make
2. Run make test
3. Run make asan
4. Run ./Ne2x_asan i:test_info.txt o:test_option.txt

Then provide a compact report:

- PASS/FAIL for each step
- first failing command output excerpt (if any)
- likely root cause
- smallest next fix

If the user provided an extra validation scenario, run it after step 4 and include it in the report.
