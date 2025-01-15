
# NeEstimator2.X

This repository contain the source code of NeEstimator v2.x software (see
[NeEstimator software](http://www.molecularfisherieslaboratory.com.au/neestimator-software/)
at [Molecular Fisheries Laboratory](http://www.molecularfisherieslaboratory.com.au/))
received from Jennifer Ovenden and maintained by NeEstimator development team
until 2019.

This folder contains the source code of NeEstimator v2.x binary software and the
java application which wraps the binary software and provide a GUI interface.

## Build (linux)

To build all the software, you require JDK 1.8 (at least) and Apache Ant. You can
build the binary software simply by:

```bash
make
```

In the project root directory. The binary software will be placed in the same
directory.
To build the GUI application, you need to enter into `NeEstimator2x`:

```bash
cd NeEstimator2x
ant
```

this will generate the `NeEstimator2x/dist/NeEstimator2x.jar` application.

## Run (linux)

To use the compiled binary with the GUI application, the binary software must be
in the same directory as the GUI application. Place the `NeEstimator2x/dist/NeEstimator2x.jar` file in the same position of the binary software. Next, the
binary software shoud be renamed into `Ne2L`. For example in the project root directory, you can do like this:

```bash
mv Ne2x Ne2L
cp NeEstimator2x/dist/NeEstimator2x.jar .
```

Then start the GUI application with:

```bash
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
