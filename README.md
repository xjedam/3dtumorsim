3dtumorsim
==========

This is a 3d tumor simulator based on CPM model.

Version
----

1.0 beta

Installation
----

To compile this application you will need:
  - GCC compiler installed - https://gcc.gnu.org/
  - Make (if not present with GCC)
  - freeglut - http://freeglut.sourceforge.net/ (note that there are some libraries present in freeglut folder - replace them with freeglut libraries aprioprate for your system)

Once you have all the above installed you can compile the application by entering `make clean` and then `make`

Usage
----

Application will be compiled to an executable file located in the `bin` directory. You must provide simulation starting state file as a parameter when running the application. For example: `bin/3dtumorsim model.dat`.

The application starts in pause mode. In order to run the simulation you need to press the "p" key.

In order to change the simulation parameters you need to edit apprioprate globals or local variables and recompile.
