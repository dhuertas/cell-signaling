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

package cellsignaling.src;

simple Cell {

	parameters:
		double emitEvery;
		double xpos;			// x coordinate initial position
		double ypos;			// y coordinate initial position
		double zpos;			// z coordinate initial position

		double vx = default(0);	// x coordinate initial velocity
		double vy = default(0);	// y coordinate initial velocity
		double vz = default(0);	// z coordinate initial Velocity

		double radius;			// Radius of the cell

		double mass;

		double listRadius = default(0);

	gates:
		input in;
		output out;

}

simple Molecule {

	parameters:
		double xpos;			// x coordinate initial position
		double ypos;			// y coordinate initial position
		double zpos;			// z coordinate initial position

		double vx;				// x coordinate initial velocity
		double vy;				// y coordinate initial velocity
		double vz;				// z coordinate initial Velocity

		double radius;			// Radius of the molecule

		double mass;

		double listRadius = default(0);

	gates:
		input in;
		output out;
}

simple Manager {
	
	parameters:
		// The C++ class for Manager has spaceSizeX, spaceSizeY and spaceSizeZ
		// attributes. Nevertheless this attributes are retrieved from the 
		// network .ned file as network parameters.

		// Manager name
		string name = default("manager");

		// Molecule Dynamics Mode. Options are: 
		// 1 - Only use Cell Lists Method
		// 2 - Combination of Near-Neighbor Lists and Cell Lists methods
		int mode = default(2);
		
		bool performanceSimulation = default(false);

	gates:
	    input in;
		output out;
}