# NICEOPS

NICEOPS (Noise Impact Computation Environment for Operations) is a tool for computing noise impact along aircraft trajectories over geographic areas. Written in the D programming language and compiled with the GCC D compiler.

NICEOPS was developed by Rob Rau and Dr. Eric Greenwood. This readme was written by Demi Zachos to document the steps necessary to compile and run NICEOPS. This is a living document for future use and reference by new and existing users.

## Compiling on Linux

### 1) Required packages

Certain packages are required to build NICEOPS. They are listed here and can be downloaded using package manager. All should be installed with -dev command.

- CGAL
- netCDF
- Boost
- libeigen3
- CMAKE
- DMD

For ubuntu as of 10/12/2021, CGAL, Boost, and libeigen3 are up to date and available. NetCDF will need to be installed from the website. Steps for these installations on ubuntu are written below:

- **I)** `sudo apt-get install libeigen3-dev`
- **II)** `sudo apt-get install libcgal-dev`
- **III)** `sudo apt-get install libboost-all-dev`
- **IV)** `sudo apt-get install libnetcdf-dev libnetcdff-dev`
- **V)** `sudo apt install cmake`
- **VI)** install DMD with the following:

  ```text
  sudo wget https://netcologne.dl.sourceforge.net/project/d-apt/files/d-apt.list -O /etc/apt/sources.list.d/d-apt.list
  sudo apt-get update --allow-insecure-repositories
  sudo apt-get -y --allow-unauthenticated install --reinstall d-apt-keyring
  sudo apt-get update && sudo apt-get install dmd-compiler dub
  ```

- **VII)** download ldc2 from github: https://github.com/ldc-developers/ldc/releases/tag/v1.27.1

  - a) unzip ldc2 using `tar xd /folder/to/ldc2-1.XX.X-linux-86_64.tar.xz`
  - b) move ldc2 to `./.local/bin` directory. You may need to make this. The command is `mkdir -p ./.local/bin` from within your home directory (navigated to using `cd ~/`)
  - c) Add symbolic link to ldc2 `ln -s /path/to/ldc2-1.XX.X-linux-x86_64/bin/ldc2 ./`
  - d) add ldc2 to your `.bashrc` (add `export PATH=$HOME/.local/bin:$PATH` to end of file)
  - e) update linux using `source .bashrc`

### 2) CMake

The command to compile NICEOPS is just cmake. However your first time might require directing the cmake files to your netcdf folder. The first time you run cmake use this command:

```text
cmake -DNETCDF=/usr/lib/x86_64-linux-gnu/libnetcdf.so ./
```

After the first time just the `cmake ./` command can be used.

### 3) Build

Build the program using `make`.

If the program doesn't build after these steps make sure to clear the cache (deleting the `CMakeFiles` folder and `CMakeCache` file) then try `cmake ./` and `make` again.

## NICEOPS basics

Run NICEOPS from within the NICEOPS folder using `./niceops`.

View the built-in manual by adding `--help` to the niceops run command.

## Command line arguments

- `-t` Tracks file location
- `-v` Database location
- `-d` .nod database (noise output database)
- `--qsam=true` Enable quasistatic aeroacoustic model
- `--loglevel=trace` Detailed logging output
- `-o` Output CSV file
- `-f` Write grid output using matplotlib (supports EPS)
- `-F` AC file location
- `-S` Save aircraft states to CSV
- `-m` Metrics file location
