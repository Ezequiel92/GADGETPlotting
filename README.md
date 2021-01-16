# GADGET Plotting

Script to produce figures, GIFs and videos from the data produce by a GAGET3/4 simulation.

- It only works with the traditional output format (binary data) which is the default in GADGET3 (SnapFormat=1) and a compatibility option in GADGET4 (legacy format selected with SnapFormat=1 too).
- It is just an script intended to be imported as is, it is not a module nor a package.
- A small testing data set is provided in testing/test_snapshots/.
- The testing script testing/testing.jl shows one example use of every function, how to import the script, and provides a sanity check, as it should run without errors.
- A Manifest.toml and Project.toml files in testing/ provide the dependencies.

## Functions

There are four tiers of functions:

- AUXILIARY FUNCTIONS: These are only for internal use, and compensate some lack of functionality in Base and other used libraries, e.g. [Unitful.jl](https://github.com/PainterQubits/Unitful.jl).
- DATA ACQUISITION FUNCTIONS: These are only for internal use. They take the raw data, apply some transformation and return it as a familiar data estructures.
- PLOTTING FUNCTIONS: These functions take data in the format outputted by the data acquisition functions and return plot objects.
- PIPELINE FUNCTIONS: These are the ones intended to be externally used. They take the location of the snapshot files, and configuration parameters, and produce a series of figures/gif/videos automatically.

NOTE: Some of the pipeline functions may produce by default a large number of images (but it can be configure to do less), and they may take a long time to run, especially if the function uses the `pgfplotsx` backend of [Plots.jl](https://github.com/JuliaPlots/Plots.jl).

## Documentation

Each function is documented within the script, where a docstring explains each function, its arguments and its returns.
For examples on how to use the functions refer to testing/testing.jl, note that it expect a particular file structure, namely:

    .
    ├── gadget_plotting.jl
    └── testing  
        ├── results
        │   ├── ...
        │   └── ...
        ├── testing.jl
        ├── Manifest.toml 
        └── Project.toml

Which is the one present in this repo.

## References

[GADGET2](https://wwwmpa.mpa-garching.mpg.de/gadget/)

[GADGET4](https://wwwmpa.mpa-garching.mpg.de/gadget4/)

## No guarantee

This script may break at any moment, and some functions are intended for data produce by GADGET3 which is not a public code. So no guaranties are given.