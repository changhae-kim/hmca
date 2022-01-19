# PyHMCA

This code is a Python wrapper of the original HMCA code.
First, compile the C code as a shared library.

    gcc -fpic -shared -o libhmca.so hmca.c

Keep a copy of `libhmca.so` and `hmca.py` in your working directory.
Then, you can import the wrapper functions in `hmca.py` to your own script.

    import hmca

For examples, see `test_*.py`.
