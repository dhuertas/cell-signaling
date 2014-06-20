# Cell Signaling Project For Omnet++

The aim of this project is to model and simulate the communication channel of signaling cells over short distances (also knows as paracrine signaling). This code is the end result of my Final Project for the Telecommunication Engineering degree at Polytechnic University of Catalonia (BCN Telecom).

![Pulse example](https://dl.dropboxusercontent.com/u/1690746/images/cell-signaling-pulse-example.png "Pulse example")

To use this code the user must have installed the Omnet++ Simulation Framework and import the code as a new project.

## Publications

### IJCTE 2014

An article has been published in International Journal of Computer Theory and Engineering (IJCTE), Volume 6, Number 4 (Aug. 2014).

Daniel Huertas González and Alfonso Rojas Espinosa, "Simulation of Cell Signaling Communications Using Event-Driven Algorithms," International Journal of Computer Theory and Engineering vol. 6, no. 4, pp. 307-312, 2014.

* ![Article at IJCTE archive](http://www.ijcte.org/index.php?m=content&c=index&a=show&catid=57&id=1051)
* ![Download article from project](https://github.com/dhuertas/cell-signaling/blob/master/doc/880-F033.pdf?raw=true)


### UPC COMMONS



## Installation

Download the project on your Omnet++ Workplace directory and export the omnet installation path. Then run `opp_makemake` to generate the Makefile:

```
export OMNET_BASE=/path/to/omnet++/folder
opp_makemake -f --deep -O out \
	-I$OMNET_BASE/src/tkenv \
    -I$OMNET_BASE/include/platdep \
    -I$OMNET_BASE/src/common \
    -I$OMNET_BASE/src/envir \
    -I/usr/include/tcl8.5 \
    -- /usr/lib/libtcl8.5.so.0
```
**Note**: project includes some of the tcl/tk libraries but they are expected to be unnecessary in a near future.

Then call `make` to compile the project:

```
make MODE=release CONFIGNAME=gcc-release all
```

## Known issues

* Emitter may enter in an infinite loop if molecules are not small enough to be added in the domain
* Some Manager parameters are not asserted and the simulation might fail if care is not taken
* Web server does not terminate gracefully when calling modules' finish function
* Web client visualization is not well oriented

## Future work

Below follows a list of things that would be nice to have, in no specific order:

* Change space cell identifiers to use datastructures using unsigned integers instead of multiple integers (one integer per direction)
* Enable the cell list algorithm to use quadtrees
* Parallelization (particle/domain decomposition)
* Add ellipsoid molecules
* Periodic boundaries support
* Preload signaling molecules during emitter initialization instead.

## Release

* r1 - Initial release

## License

View LICENSE file.

