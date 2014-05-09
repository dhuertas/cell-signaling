
package cellsignaling.networks.tests;

import cellsignaling.src.Cell;
import cellsignaling.src.Manager;
import cellsignaling.src.Molecule;
import cellsignaling.src.Probe;

network domain
{
	@display("bgb=400,300");

	double spaceSizeX = default(1000);
	double spaceSizeY = default(1000);
	double spaceSizeZ = default(1000);

	int numberOfInitialMolecules = default(0);

	submodules:
		manager: Manager {
			@display("p=-30,-30;i=block/cogwheel");
		}

		cell[1]: Cell;
		molecule[numberOfInitialMolecules]: Molecule;
		probe[10]: Probe;

	connections allowunconnected:

}
