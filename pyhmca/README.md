# PyHMCA

This is a Python wrapper of the original HMCA code in C.
First, you need to compile the C code as a shared library.

  gcc -fpic -shared -o libhmca.so hmca.c
