//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with this program.  If not, see http://www.gnu.org/licenses/.
// 

package cellsignaling.networks;

import cellsignaling.src.Cell;
import cellsignaling.src.Manager;
import cellsignaling.src.Molecule;

network test5
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
