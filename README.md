# Cell Signaling Project For Omnet++

The aim of this project is to model and simulate the communication channel of signaling cells over short distances (also knows as paracrine signaling). This code is the end result of my Final Project for the Telecommunication Engineering degree at Polytechnic University of Catalonia (BCN Telecom).

It uses a hard-sphere model and brownian motion to reproduce the diffusion mechanism of signaling molecules in a fluid. It also includes a web based 3D visualization component using [three.js](http://threejs.org) library and JSON streaming.

Note that to use this code the user must have installed the Omnet++ Simulation Framework and import the code as a new project.

|Pulse received by a cell at 250nm from the emitter|
|:------------------------------------------------:|
|![Pulse-example](https://dl.dropboxusercontent.com/u/1690746/images/cell-signaling-pulse-example.png "Pulse example")|

|Snapshot of the 3D visualization component using Chromiun web browser|
|:-------------------------------------------------------------------:|
|![Snapshot](https://dl.dropboxusercontent.com/u/1690746/images/webclient-capture.png "Snapshot")|

## Publications

### IJCTE 2014

An article has been published in International Journal of Computer Theory and Engineering (IJCTE), Volume 6, Number 4 (Aug. 2014).

Cite: Daniel Huertas González and Alfonso Rojas Espinosa, "Simulation of Cell Signaling Communications Using Event-Driven Algorithms," International Journal of Computer Theory and Engineering vol. 6, no. 4, pp. 307-312, 2014.

* [Article at IJCTE archive](http://www.ijcte.org/index.php?m=content&c=index&a=show&catid=57&id=1051)
* [Download article from project](https://github.com/dhuertas/cell-signaling/blob/master/doc/880-F033.pdf?raw=true)


### UPC COMMONS

Final project (to be added).

## Installation

Download the project to your Omnet++ workplace directory and export the omnet installation path:

```
cd /path/to/workplace/folder
git clone https://github.com/dhuertas/cell-signaling ./cell-signaling
cd /path/to/omnet++/folder
. setenv
```

Start Omnet++ IDE (Eclipse) and create an empty Omnet++ project. Set the folder project to use the cell-signaling folder. To build the source code the project must be configured to link against libpthread:

1. Go to project properties and select Omnet++ section. 
2. Go to Makemake, select the source folder and choose "Makemake" option. 
3. Click the options button and navigate to the "Link" tab.
4. Add "pthread" to the link box.

Finally build and run the simulation environment.

## Known issues

* Emitter may enter in an infinite loop if molecules are not small enough to be added in the domain
* Some Manager parameters are not asserted and the simulation might fail if care is not taken
* Web server does not terminate gracefully when calling modules' finish function
* Web client visualization is not well oriented

## Future work

Below follows a list of things that would be nice to have, in no specific order:

* Change space cell identifiers to use datastructures using unsigned integers instead of multiple integers (one integer per direction)
* Enable the cell list algorithm to use octrees
* Parallelization (particle/domain/event decomposition)
* Add ellipsoid molecules
* Add periodic boundaries

## Release

* r1 - Initial release

## License

View LICENSE file.

