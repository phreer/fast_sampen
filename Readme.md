# fast_sampen: Fast Computation of Sample Entropy

This repository includes a library and a program for fast computation of Sample Entropy, based on kd tree, Range-kd tree and randomly sample (Monte Carlo and quasi-Monte Carlo) method.

## Requirements

- Linux, macOS or other UNIX-like OS
- C++ compiler supporting C++11
- CMake (version >= 3.5)
- GSL (for quasi-random number generation)
- Magick++7 (for image manipulation such as io and resize)

## Compile

### Step 1 Install the GNU Scientific Library (GSL)

We need the quasi-random number generator provided by GSL, so the first step is to install it. You have two options to
install GSL: install from source or install via package manager.

#### Install via Package Manager

On macOS, if HomeBrew is available, issue the following command in terminal to install GSL.

```bash
brew install gsl
```

On Ubuntu or other Linux distributions, replace `brew` to the package manager such as `apt` or `yum` shipped by the
distribution.

#### Install from Source

Download source code of GSL from [this cite](https://www.gnu.org/software/gsl/)
and issue the following commands to compile and install it (Suppose that version 2.4 is downloaded).

```bash
tar -zxvf gsl-2.4.tar.gz
cd gsl-2.4
./configure
# Or if administrative privilege is not available, add --prefix argument to 
# specify the location to install GSL, i.e., 
# ./configure --prefix=<INSTALL_LOCATION>
make -j4 # Issue 4 threads for compiling.
sudo make install
# If administrative privilege is unavailable, use the following instead:
# make install
```

### Step 2 Compile the Program

Once the GSL is installed, you can use the following command to compile the
`sampen_kdtree` program.

```bash
# Clone this repository and enter the root directory of it.
git clone https://github.com/phreer/sampen_kdtree.git
cd sampen_kdtree
# Create a building directory and enter it.
mkdir build
cd build
# Run cmake to generate Makefile.
cmake -DCMAKE_BUILD_TYPE=Release ..
# This command generates for release build. If you need a debug build, just
# change Release to Debug. If a customized C++ compiler is needed, add argument 
# -DCMAKE_CXX_COMPILER=<PATH_TO_COMPILER> to the command. If the GSL is not 
# installed system-widely, add -DGSL_ROOT_PATH=<GSL_INSTALL_LOCATION>.
# Finally, run make to execute compilation.
make -j4
```

The built binary will be located in `bin` directory in the building directory.

## Usage

Run programs with `--help` for usage.
