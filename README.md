
# NeEstimator2.X

This repository contain the source code of NeEstimator v2.x software (see
[NeEstimator software](http://www.molecularfisherieslaboratory.com.au/neestimator-software/)
at [Molecular Fisheries Laboratory](http://www.molecularfisherieslaboratory.com.au/))
received from Jennifer Ovenden and maintained by NeEstimator development team
until 2019.

This folder contains the source code of NeEstimator v2.x binary software and the
java application which wraps the binary software and provide a GUI interface.

## Porting intent

The current work introduces an incremental C-to-C++ migration path focused on
safety and maintainability, while preserving scientific behavior.

In this phase:

- `Ne2x` remains the original C executable.
- `main.cpp` is the minimal C++ entrypoint forwarding `argc/argv` to the same C core.
- `Ne2x_cpp` is the executable produced from `main.cpp` and linked against the same C core.
- `libne2x.a` exposes the C core so C and C++ binaries execute the same logic.

This allows modernization (tooling, diagnostics, future refactors) without
rewriting the estimation algorithms up front.

## Build (linux)

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

## Run with the GUI (linux)

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
