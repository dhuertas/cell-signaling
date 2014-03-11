
package cellsignaling.networks.tests;

import cellsignaling.src.SimpleCell;
import cellsignaling.src.Manager;
import cellsignaling.src.Molecule;

network BasicTest
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
		emitter: SimpleCell;
		receiver: SimpleCell;
		molecule: Molecule;

	connections allowunconnected:

}
