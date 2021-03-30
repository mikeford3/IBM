Interior Ballistics Model

Currently version 0.1, a rough framework that sort-of works.

To use first compile the code using make on a system with a C++ compiler (tested on linux only with GCC), and then run the executable created, ibm.

Todo:
    * Read input file with parameters for Gun/Shell/Propellant/Gas.
    * Update physical laws, especially everything related to burn rate and dynamic friction.
    * Consider probabilistic model if the inputs can be specified to a reasonable level of precision,
      would need to take temperature, burn rate, energy release etc as distributions.
