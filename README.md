# Joint optimization of distortion and cut location for mesh parameterization using an Ambrosio-Tortorelli functional

Source code of the paper: "[Joint optimization of distortion and cut location for mesh parameterization using an Ambrosio-Tortorelli functional](https://perso.liris.cnrs.fr/david.coeurjolly/publication/uv-at/uv-at.pdf)", Colin Weill–Duflos, David Coeurjolly, Fernando de Goes and Jacques-Olivier Lachaud, Computer Aided Geometric Design (proc. GMP 2023), 2023.

![IMG_0384](https://github.com/Lieunoir/UV-AT/assets/700165/3cbaa97c-ef53-4914-aee8-293fdd83ea95)


Authors of the code:
* Colin Weill-Duflos
* Jacques-Olivier Lachaud
* David Coeurjolly


## Building

Once cloned, fetch the dependencies (geometry-central, libigl, polyscope and DGtal) using git submodules:

```
git submodule update --init --recursive
```
or directly,
```
git clone --recursive https://github.com/Lieunoir/UV-AT.git
```


Then, on  unbutu/macos:

```
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j 8
```


## Usage

To run : `./build/bin/at-uv -i bunnyhead.obj`

You can use the button `Compute initial AT-UV` to compute the initial parameterization and cut and then run the method once. It results with a uv mapping and its corresponding v. It then enables cutting buttons.

`Cut` is used to compute the cuts from v.

When a cut is done, we can either use `Refine AT-UV` to run the method again from the result and obtain a new compute uv and v, or `Finalize AT-UV` to get a final result (force v at 1).

`Compute whole` is a shortcut for initialization, running our method once, cutting and finalizing.

The option `--withoutGUI` can also be used to run the method in headless mode, and put the results in an `output` directory. Use `./batch.sh ./path/to/obj/folder` can be used to run the method over all the obj files contained in a mesh. `./generate_tex.sh` can then be used to generate a tex file containing the results.
