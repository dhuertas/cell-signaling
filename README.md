# Cell Signaling Project For Omnet++

The aim of this project is to model and simulate the communication channel of signaling cells over short distances (also knows as paracrine signaling). This code is the end result of my Final Project for the Telecommunication Engineering degree at Polytechnic University of Catalonia (BCN Telecom).

To use this code the user must have installed the Omnet++ Simulation Framework and import the code as a new project.

## Installation

Download the project on your Omnet++ Workplace directory and run opp_makemake togenerate the Makefile:

```
opp_makemake -f --deep -O out -I/home/dani/Applications/omnetpp-4.3.1/src/tkenv \
     -I/usr/include/tcl8.5 \
     -I/home/dani/Applications/omnetpp-4.3.1/include/platdep \
     -I/home/dani/Applications/omnetpp-4.3.1/src/common \
     -I/home/dani/Applications/omnetpp-4.3.1/src/envir \
     -- /usr/lib/libtcl8.5.so.0
```
Then call `make` to compile the project.

## Issues

* Change space cell identifiers to datastructures using unsigned integers instead

## Future work

## Release

* r1 - Initial release

## License

TBD

