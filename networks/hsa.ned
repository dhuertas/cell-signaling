
package cellsignaling.networks;

import cellsignaling.src.Cell;
import cellsignaling.src.Manager;
import cellsignaling.src.Molecule;

network hsa
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

		molecule[numberOfInitialMolecules]: Molecule;

	connections allowunconnected:

}
