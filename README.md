
# NeEstimator2.X

This repository contains the source code of NeEstimator v2.x software (see
[NeEstimator software](http://www.molecularfisherieslaboratory.com.au/neestimator-software/)
at [Molecular Fisheries Laboratory](http://www.molecularfisherieslaboratory.com.au/))
received from Jennifer Ovenden and maintained by the NeEstimator development team
until 2019.

This folder contains the source code of NeEstimator v2.x binary software and the
Java application, which wraps the binary software and provides a GUI interface.

## Documentation

The documentation on NeEstimator is available in the [Help.pdf](docs/Help.pdf) file,
inside the [docs](docs) folder of this repository.

## Download pre built binaries

Pre-built binaries for Linux, Windows, and macOS are available in the
[Releases](https://github.com/bunop/NeEstimator2.X/releases) section of this
repository, with the compiled Java GUI application included as well:

| Asset | Description |
|--------|-----------|
| Ne2.exe | Windows executable |
| Ne2L | Linux executable |
| Ne2M | macOS executable |
| NeEstimator2x.jar | Java GUI application |
| Source code (zip) | Source code (zip) |
| Source code (tar.gz) | Source code (tar.gz) |

Simply download the appropriate executable for your platform and place it in the same
folder as `NeEstimator2x.jar` to use the GUI application (java required). For CLI usage,
the executable can be used directly without the GUI.

## Build

### Linux

From project root:

```bash
make
```

This builds `Ne2x` (original C executable) in the current folder.

### Windows

Install MinGW (same tooling used in CI):

```bash
choco install mingw
```

From project root, run:

```bash
mingw32-make
```

If your setup exposes `make` as an alias, this is equivalent:

```bash
make
```

### macOS

Install GCC (same tooling used in CI):

```bash
brew install gcc
```

From project root:

```bash
make
```

### Java GUI

To build the GUI application (legacy Java wrapper), enter `NeEstimator2x` folder,
then compile with `ant`:

```bash
cd NeEstimator2x
ant
```

This generates `NeEstimator2x/dist/NeEstimator2x.jar`.

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

## CLI invocation modes

> [!TIP]
> The easiest way to run it via the CLI is to start from the GUI Java
> application: set all your parameters there and then click on the "Create Parameter Files"
> button. This will generate a single `info` file that can be used in single-run mode
> (`i:`).

The program can be invoked in **5 high-level modes**:

1. No arguments (interactive/direct mode)
2. `i:` mode (single run from info directive file)
3. `m:` mode (multi-file batch, short syntax)
4. `m+:` mode (multi-file batch, extended syntax)
5. `c:` mode (common-settings batch mode)

### Syntax summary

```text
./Ne2x
./Ne2x i:<info_file> [o:<option_file>] [rm]
./Ne2x m:<multi_file> [rm]
./Ne2x m+:<multi_file> [rm]
./Ne2x c:<common_file> [rm]
```

### Meaning of each mode

1. `./Ne2x`
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

## Troubleshooting and known limitations (from GitHub issues)

This section summarizes recurring problems reported in
[Issues](https://github.com/bunop/NeEstimator2.X/issues) and practical
workarounds.

### GUI does not start: `HeadlessException` / missing `$DISPLAY`

Reference: [#7](https://github.com/bunop/NeEstimator2.X/issues/7)

The Java GUI requires an active X server. If you run over SSH or on a headless
node, `java -jar NeEstimator2x.jar` can fail with `HeadlessException`.

Suggestions:

- Run GUI locally on a desktop session.
- If remote access is required, use X forwarding (for example `ssh -X`) and
	verify `$DISPLAY` is set.
- For clusters/HPC, prefer CLI mode with directive files (`i:`, `m:`, `m+:`,
	`c:`).

### CLI/directive files are strict and can fail silently

References: [#6](https://github.com/bunop/NeEstimator2.X/issues/6),
[#7](https://github.com/bunop/NeEstimator2.X/issues/7),
[#8](https://github.com/bunop/NeEstimator2.X/issues/8)

Directive-file syntax is sensitive. Common pitfalls:

- Directory paths in directive files should end with `/`.
- Inline comments in example files use `*` (not shell-style `#`).
- Spacing/formatting in control files must follow examples closely.

Recommended validation command:

```bash
./Ne2x i:test_info.txt o:test_option.txt
```

If this succeeds but your own run fails, the issue is likely in your custom
input/control files.

### Temporary-file pressure on `/tmp` with large runs

Reference: [#1](https://github.com/bunop/NeEstimator2.X/issues/1)

The C executable uses temporary-file based strategies that can create high disk
pressure on temporary storage for large jobs.

Practical advice:

- Monitor free space on the filesystem hosting temporary files (often `/tmp`).
- On HPC, avoid running where temporary storage is small or shared and close to
	quota.
- If temporary storage fills up, output quality may degrade (for example
	unexpected `NaN`/`Inf` values).

### Installation clarity by platform

Reference: [#2](https://github.com/bunop/NeEstimator2.X/issues/2)

Confirmed by user reports, both binary and `.jar` files are needed for GUI-driven runs:

- Windows: `Ne2.exe` and `NeEstimator2x.jar` in the same folder.
- macOS: `Ne2M` and `NeEstimator2x.jar` in the same folder.
- Linux: `Ne2L` and `NeEstimator2x.jar` in the same folder.

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
