INSTALLATION
============

The `ldpsiz` package is written in C and Python. It runs under Linux
and OSX. I have not tried to compile it under Windows.

The package is available at [github](github.com/alanrogers/ldpsiz) in
various formats. Before compiling, you must install two libraries:
`pthreads` and [`gsl`](http://www.gnu.org/software/gsl). You will need
not only the libraries themselves but also several header files, such
as `pthread.h`. I didn't need to install `pthreads`, because it came
bundled with the Gnu C compiler. But the gsl was an extra. Under ubuntu
Linux, you can install it like this:

    sudo apt-get install libgsl0-dev

On the mac, using homebrew, the command is

    brew install gsl

By default, the executable files will be copied into a directory named
`bin` in your home directory. If you want them to go somewhere else,
edit the first non-comment line of src/Makefile.

Then 

1. Cd into the src directory.
2. Type "make".
3. Type "make install".

This will try to place the executables into directory "bin" in the
user's home directory. Make sure this directory appears in your
PATH, so that the shell can find it.

This installation will work under unix-like operating systems, such as
linux and Apple's osx. I haven't tried to port this software to
Windows. 

The directory `test` contains a unit test for many of the .c files in
directory `src`. Within this directory, type

1. make xhill
2. ./xhill

to test the source file `hill.c`.  To run all unit tests, type
"make". This will take awhile, as some of the unit tests are slow.

The unit tests will not compile if `NDEBUG` is defined. If
optimization is turned on during optimization, some of the unit tests
may be removed by the optimizer. To avoid these problems, comment out
the relevant line at the top of src/Makefile.



