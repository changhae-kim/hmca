# PyHMCA

## Overview

This is a Python wrapper of the HMCA code.
First, compile the C code as a shared library.

    gcc -fpic -shared -o libhmca.so hmca.c

Keep a copy of `libhmca.so` and `hmca.py` in your working directory.
Then, you can import the wrapper functions in `hmca.py` to your own script.

    import hmca

For examples, see `test_*.py`.

Like the C code, the Python wrapper only computes the right-hand side and the Jacobian.
To use this code in a chemical kinetic simulation, you need to pick a suitable numerical integrator.
Fortunately, `scipy` has a Python wrapper of CVODE.

More detailed explanations of the usage to be added . . .
