# sampen_kdtree: Fast Computation of Sample Entropy
This repository includes a program for fast computation of Sample Entropy, 
based on KD tree and randomly sampling method.

## Requirements
- Linux, macOS or other UNIX-like OS
- C++ compiler supporting C++11
- CMake (version >= 3.5)

## Compile
### Step 1 Install the GNU Scientific Library (GSL)
We need the quasi-random number generator provided by GSL, so the first step 
is to install it.
You have two options to install GSL: install from source or install via 
package manager.

#### Install via Package Manager
On macOS, if HomeBrew is available, issue the following command in terminal
to install GSL.
```bash
brew install gsl
```

On Ubuntu or other Linux distributions, replece `brew` to the package manager 
such as `apt` or `yum` shipped by the distribution.

#### Install from Source
Download source code of GSL from [this cite](https://www.gnu.org/software/gsl/)
and issue the following commands to compile and install it (Suppose that 
version 2.4 is downloaded).
```bash
tar -zxvf gsl-2.4.tar.gz
cd gsl-2.4
./configure
# Or if administrative privilege is not availible, add --prefix argument to 
# specify the location to install GSL, i.e., 
# ./configure --prefix=<INSTALL_LOCATION>
make -j4 # Issue 4 threads for compling.
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
```
Usage: build/kdtree_sample --input <INPUT> --input-type {simple, multirecord}\
                   -r <THRESHOLD> -m <TEMPLATE_LENGTH>\
                   -n <N> [-output-level {1,2,3}]\
                   --sample-size <SAMPLE_SIZE> --sample-num <SAMPLE_NUM>

Options and arguments:
--input <INPUT>         The file name of the input file.
--input-format <FORMAT> The format of the input file. Should be either simple
                        or multirecord. If set to simple, then each line of the
                        input file contains exactly one column; if set to
                        multirecord, then each line contains <NUM_RECORD> + 1
                        columns, of which the first indicates the line number
                        and the remaining columns are instances of the records.
                        The default value is simple.
--input-type <TYPE>     The data type of the input data, either int or float.
                        Default: double.
-r <R>                  The threshold argument in sample entropy.
-m <M>                  The template length argument of sample entropy. Note
                        that this program only supports 2 <= m <= 10.
-n <N>                  If the length of the signal specified by <FILENAME> is
                        greater than <N>, then it would be truncated to be of
                        length <N>. If <N> is 0, then the the original length
                        is employed. The default value is 0.
--sample-size <N0>      The number of points to sample.
--sample-num <N1>       The number of computations where the average is taken.
--random                If this option is enabled, the random seed will be set
                        randomly.
-q                      If this option is enabled, the quasi-Monte Carlo based
                        method is conducted.
--variance              If this option is enabled, then the variance of the
                        results of sampling methods will be computed.
--quasi-type <TYPE>     The type of the quasi-random sequence for sampling,
                        can be one of the following: sobol, halton,
                        reversehalton, niederreiter_2 or grid. Default: sobol.
-u                      If this option is enabled, the Monte Carlo based
                        using uniform distribution is conducted.
--output-level <LEVEL>  The amount of information printed. Should be one of
                        {0,1,2}. Level 0 is most silent while level 2 is for
                        debugging.
```