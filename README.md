# GADGET Plotting

Script for the generation of figures, GIFs and videos from the data produce by a GAGET3/4 simulation.

- It only works with the traditional output format (binary data) which is the default in GADGET3 (SnapFormat=1) and a compatibility option in GADGET4 (legacy format selected with SnapFormat=1 too).
- It is just an script intended to be imported as is, it is not a module nor a package.
- A small testing data set is provided in testing/test_snapshots/.
- The testing script testing/testing.jl shows examples of every function, how to import the script, and provides a sanity check, as it should run without errors.
- A testing/Manifest.toml and testing/Project.toml files provide the dependencies.

## Functions

There are four tiers of functions:

- AUXILIARY FUNCTIONS: This are only for internal used, and compensate lack of functionality in Base or the used libraries, e.g. Unitful.jl, or isolate a bit a repeated logic.
- DATA ACQUISITION FUNCTIONS: This are only for internal used. They take the raw data, applied some transformation and return it as familiar data estructures.
- PLOTTING FUNCTIONS: This functions take data in the format outputted by the data acquisition functions and return plot objects.
- PIPELINE FUNCTIONS: This are the ones intended to be externally used. They take the location of the snapshot files and configuration parameters to produce a series of figures/gif/videos automatically.

NOTE: Some of the pipeline functions may produce by default a large number of images (but it can be configure to do less). And they may take a long time to run, especially if the function utilices the pgfplotsx backend of Plots.jl.

## Documentation

Each function is documented within the script, where the docstring explains each function, its arguments and its returns.
For examples refer to testing/testing.jl, note that it expect a particular structure of its files, namely:

    .
    ├── gadget_plotting.jl
    └── testing  
        ├── results
        │   ├── ...
        │   └── ...
        ├── testing.jl
        ├── Manifest.toml 
        └── Project.toml

Which is the one present in the repo.

### References

[GADGET2](https://wwwmpa.mpa-garching.mpg.de/gadget/)

[GADGET4](https://wwwmpa.mpa-garching.mpg.de/gadget4/)

### No guarantee

This script may break at any moment, and some functions are intended for data produce by GADGET3 which is not a public code. So no guaranties are given.