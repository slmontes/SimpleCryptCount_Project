# SimpleCryptCount_Project

## Crypt count algorithm using segmented images

This algorithm was created to approximate the number of crypts on intestinal organoids by analysing boundaries extracted from segmented organoid images.
In our paper: "In-silico and in-vitro morphometric analysis of intestinal organoids", we applied this algorithmm to both _in-silico_ and _in-vitro_ images.

The main algorithm function is: CountingCrypts_wCircularityFun.m
The following parameters are inputs required:
- Type: Choose between "In vitro" or "In silico"
- MaskSet_name: Name of the organoid mask to read
Parameter values to identify crypts:
- Fourier_harmonic_term: See Table 2 in the paper
- Crypt_parameters: See Table 2 in the paper (Minimum crypt area, maximum
                    crypt area and minimum arc length)

## Model of crypt formation on 2D intestinal organoids with 4 different cell types 

Before looking at this, you may wish to look at some of the https://chaste.cs.ox.ac.uk/trac/wiki/UserTutorials for Chaste.

### Getting the code and installing dependencies 

Before running this model you will need to install Chaste following these instructionds: https://chaste.cs.ox.ac.uk/trac/wiki/InstallGuides/InstallGuide
The easiest way to do this is using an Ubuntu machine (or an Ubuntu virtual machine) as discussed on https://chaste.cs.ox.ac.uk/trac/wiki/InstallGuides/UbuntuPackage
Note that Chaste is only fully supported on Linux/Unix systems, so users of Windows or Mac OS X may need to follow the virtual machine route.
For manual installation of each dependency, on any version of Linux, see DeveloperInstallGuide. Alternatively, Chaste can also be run using Docker, and the dockerfile can be accessed at https://github.com/Chaste/chaste-docker, where you'll find instructions for installation and usage.

You will also need the source for the CryptFission4CellTypes project.  This can be done by checking out the version from the repository by using the command 
```sh
git clone https://github.com/slmontes/SimpleCryptCount_Project.git 
```
and add the CryptFission4CellTypes folder found in 2D Model to the projects folder of the Chaste directory:

```sh
Chaste/projects/CryptFission4CellTypes/TestCryptFission4CellTypes_EC10toTA70.hpp
```

### Documentation 
There are two folders - `src` and `test`.
 1. The `src` folder contains the classes necesary to run the simulation. These define the aditional forces and boundary conditions not in the core chaste code.
 1. The `test` folder contains:
  - TestCryptFission4CellTypes_EC10toTA70.hpp this file can be run to generate the in-silico results
  
### Running tests
You can then run tests and simulations from the terminal with:
```sh
cd <Chaste_build path>
make -j4 TestCryptFission4CellTypes_EC10toTA70
ctest -j4 -V -R TestCryptFission4CellTypes_EC10toTA70
```
----
> Note: the paper was developed with release version 2021, but will work on release version 2019 and 2018.1. It will not work with with release version 3.4 or under.

For further information on using Chaste, see https://chaste.cs.ox.ac.uk/trac/wiki/ChasteGuides.
You may also wish to look at some of the https://chaste.cs.ox.ac.uk/trac/wiki/UserTutorials.

### Dataset
In-vitro images used can be found at: https://doi.org/10.6084/m9.figshare.23301137.v1
