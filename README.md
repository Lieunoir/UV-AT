# UV-AT

## Usage

To build : `mkdir build; cd build; cmake ..; make`
To run : `./build/bin/at-uv -i bunnyhead.obj`

You can use the button `Compute initial AT-UV` to compute the initial parameterization and cut and then run the method once. It results with a uv mapping and its corresponding v. It then enables cutting buttons.

`Cut` is used to compute the cuts from v.

When a cut is done, we can either use `Refine AT-UV` to run the method again from the result and obtain a new compute uv and v, or `Finalize AT-UV` to get a final result (force v at 1).

`Compute whole` is a shortcut for initialization, running our method once, cutting and finalizing.
