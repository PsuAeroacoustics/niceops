NICEOPS was developed by Rob Rau and Dr. Eric Greenwood. This readme was written by Demi Zachos to document the steps necessary to compile and run NICEOPS. This is a living document for future use and reference by new and existing users.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Compiling on Linux
-------------------------------------------
1) Certain packages are required to build NICEOPS. They are listed here and can be downloaded using package manager. All should be installed with -dev command.
	a) CGAL
	b) netCDF
	c) Boost
	d) libeigen3

For ubuntu as of 10/12/2021, CGAL, Boost, and libeigen3 are up to date and available. NetCDF will need to be installed from the website. Steps for these installations on ubuntu are written below:
		I)   sudo apt-get install libeigen3-dev
		II)  sudo apt-get install libcgal-dev
		III) sudo apt-get install libboost-all-dev
		IV)  sudo apt-get install libnetcdf-dev libnetcdff-dev
		V)   download ldc2 from github: https://github.com/ldc-developers/ldc/releases/tag/v1.27.1
			a) unzip ldc2 using "tar xd /folder/to/ldc2-1.XX.X-linuxx86_64.tar.xz"
			b) move ldc2 to ./.local/bin directory
				You may need to make this. The command is "mkdir -p ./.local/bin" from within your home directory (navigated to using cd  ~/)
			c) add ldc2 to your .bashrc (add "export PATH=$HOME/.local/bin:$PATH")
			d) update linux using "source .bashrc"

2) The command to compile NICEOPS is just cmake. However your first time might require directing the cmake files to your netcdf folder. THe first time you run cmake use this command:
"cmake -DNETCDF=/usr/lib/x86_64-linux-gnu/libnetcdf.so ./"

After the first time just the "cmake ./" command can be used.

3) Build the program using "make"

If the program doesn't build after these steps make sure to clear the cache (deleting CMakeFiles folder and CMackeCache file) then try "cmake ./" and "make" again.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NICEOPS basics
-------------------------------------------
Run NICEOPS from within the NICEOPS folder using "./niceops" command.
View the built-in manual by adding "--help" to the niceops run command.
		
