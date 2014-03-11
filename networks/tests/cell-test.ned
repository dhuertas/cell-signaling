
package cellsignaling.networks.tests;

import cellsignaling.src.Cell;
import cellsignaling.src.Manager;
import cellsignaling.src.Molecule;

network CellTest
{
	@display("bgb=400,300");

	double spaceSizeX = default(1000);
	double spaceSizeY = default(1000);
	double spaceSizeZ = default(1000);

	int numberOfInitialMolecules = default(0);

	submodules:
		manager: Manager {
			@display("p=30,30;i=block/cogwheel");
		}
		emitter: Cell;
		receiver: Cell;
		molecule[10]: Molecule;

	connections allowunconnected:

}
