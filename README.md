
# NeEstimator2.X

This repository contains the source code of NeEstimator v2.x software (see
[NeEstimator software](http://www.molecularfisherieslaboratory.com.au/neestimator-software/)
at [Molecular Fisheries Laboratory](http://www.molecularfisherieslaboratory.com.au/))
received from Jennifer Ovenden and maintained by the NeEstimator development team
until 2019.

This folder contains the source code of NeEstimator v2.x binary software and the
Java application, which wraps the binary software and provides a GUI interface.

## Porting intent

The current work introduces an incremental C-to-C++ migration path focused on
safety and maintainability, while preserving scientific behavior.

In this phase:

- `Ne2x` remains the original C executable.
- `main.cpp` contains a C++ CLI dispatcher equivalent to the legacy C entrypoint,
  and forwards execution to the same C core routines.
- `Ne2x_cpp` is the executable produced from `main.cpp` and linked against the same C core.
- `libne2x.a` exposes the C core so C and C++ binaries execute the same logic.

This allows modernization (tooling, diagnostics, future refactors) without
rewriting the estimation algorithms up front.

## Build (Linux)

From project root:

```bash
make
```

This builds:

- `Ne2x` (original C executable)
- `libne2x.a` (C static library)
- `Ne2x_cpp` (C++ wrapper executable)

To build the GUI application (legacy Java wrapper), enter `NeEstimator2x` and run:

```bash
cd NeEstimator2x
ant
```

This generates `NeEstimator2x/dist/NeEstimator2x.jar`.

## AddressSanitizer check

To build an ASan-instrumented binary:

```bash
make asan
```

This produces `Ne2x_asan`, which is used to detect memory errors such as
out-of-bounds accesses, use-after-free, and leaks.

Run a representative scenario:

```bash
./Ne2x_asan i:test_info.txt o:test_option.txt
```

Expected outcome:

- process exits successfully
- no AddressSanitizer/LeakSanitizer report is printed

## Testing C vs C++ binaries

Use the regression script to verify that C and C++ wrapper outputs are identical
for the reference test case:

```bash
make test
```

or equivalently:

```bash
bash regression_test.sh
```

The script:

1. builds `Ne2x` and `Ne2x_cpp`
2. runs both with `i:test_info.txt o:test_option.txt`
3. compares stdout and generated output files (`test_Ne.txt`, `test_NexLD.txt` when present)

Note: the comparison ignores runtime-only timestamp lines (for example
`Starting time:`, `Ending time:`, `Time:`), so the test focuses on scientific
outputs rather than execution clock differences.

Expected success message:

```text
PASS: C and C++ wrapper outputs are identical for test_info/test_option.
```

## CLI invocation modes

Based on the current dispatcher in `main.cpp`, the program can be invoked in
**5 high-level modes**:

1. No arguments (interactive/direct mode)
2. `i:` mode (single run from info directive file)
3. `m:` mode (multi-file batch, short syntax)
4. `m+:` mode (multi-file batch, extended syntax)
5. `c:` mode (common-settings batch mode)

### Syntax summary

```text
./Ne2x_cpp
./Ne2x_cpp i:<info_file> [o:<option_file>] [rm]
./Ne2x_cpp m:<multi_file> [rm]
./Ne2x_cpp m+:<multi_file> [rm]
./Ne2x_cpp c:<common_file> [rm]
```

Equivalent forms are accepted by `Ne2x` as well.

### Meaning of each mode

1. `./Ne2x_cpp`
	Runs direct/interactive flow (same as legacy behavior when no CLI argument
	is provided).

2. `i:<info_file>`
	Runs one analysis from the info directive file. You may also provide
	`o:<option_file>` to supplement options.

3. `m:<multi_file>`
	Runs batch processing using the short multi-file directive format.

4. `m+:<multi_file>`
	Runs batch processing using the extended multi-file directive format.

5. `c:<common_file>`
	Runs batch processing with common settings followed by input-file entries.

### Optional `rm` argument

`rm` means "remove control file(s) after use":

- In `m:`, `m+:`, and `c:` modes, `rm` requests removing the first control file.
- In `i:` mode, `rm` can appear as the second argument (if no `o:` is used) or
  as the third argument (after `o:`).

Notes:

- The parser accepts only command families starting with `i:`, `m:`, `m+:`, or `c:`.
- Any malformed command prefix prints `Illegal argument!`.
- The C++ dispatcher intentionally preserves legacy behavior for compatibility.

## AI customization for this project (VS Code)

This repository includes project-level customization files for coding agents:

- Instruction: `.github/instructions/c-memory-safety.instructions.md`
- Prompt: `.github/prompts/run-core-validation.prompt.md`
- Skill: `.github/skills/investigate-regression-mismatch/SKILL.md`

How to use them:

1. `c-memory-safety.instructions.md`
	- Auto-applies on C/C++ files (`*.c`, `*.h`, `*.cpp`).
	- Use it when touching memory-sensitive code paths, especially in `Ne2x.c`.
	- It enforces a safety checklist and the expected validation sequence (`make`, `make test`, `make asan`, `./Ne2x_asan ...`).

2. `/run-core-validation`
	- Type `/` in chat and select `Run Core Validation`.
	- Runs the standard validation flow and returns a compact PASS/FAIL summary with the first failing excerpt and likely root cause.
	- You can provide an optional extra scenario argument.

3. `/investigate-regression-mismatch`
	- Type `/` in chat and select `investigate-regression-mismatch`.
	- Use it when `make test` fails or C/C++ outputs diverge.
	- It classifies differences (timestamp-only vs real scientific mismatch), points to likely files, proposes the smallest safe fix, and revalidates with `make test`.

Tip: keep these customizations in sync with `AGENTS.md` when build/test workflows change.

## Run with the GUI (Linux)

To use the compiled binary with the GUI application, the binary software must be
in the same directory as the GUI application. Place
`NeEstimator2x/dist/NeEstimator2x.jar` in the same position as the binary.
Then rename the binary to `Ne2L`.

For example, from project root:

```bash
mv Ne2x Ne2L
cp NeEstimator2x/dist/NeEstimator2x.jar .
java -jar NeEstimator2x.jar
```

## Citation

If you use this code, you should cite the original methods as well as NeEstimator
program note. For example:

> “We estimated Ne using the molecular co-ancestry method of Nomura (2008),
> as implemented in NeEstimator V2.1 (Do et al. 2014.).”
>
> Do, C., Waples, R. S., Peel, D., Macbeth, G. M., Tillett, B. J. & Ovenden,
> J. R. (2014). NeEstimator V2: re-implementation of software for the
> estimation of contemporary effective population size (Ne) from genetic data.
> Molecular Ecology Resources. 14, 209-214.
