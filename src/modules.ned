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
		double xpos;		// x coordinate initial position
		double ypos;		// y coordinate initial position
		double zpos;		// z coordinate initial position

		double vx = default(0);
		double vy = default(0);
		double vz = default(0);

		double radius;		// Radius of the cell

		double mass;

	gates:
		input in;
		output out;

}

simple Molecule {

	parameters:
	    double xpos;		// x coordinate initial position
		double ypos;		// y coordinate initial position
		double zpos;		// z coordinate initial position

		double radius;		// Radius of the molecule

		double vx;			// x coordinate initial velocity
		double vy;			// y coordinate initial velocity
		double vz;			// z coordinate initial Velocity

		double mass;

	gates:
		input in;
		output out;
}

simple Manager {
	
	parameters:
		// The C++ class for Manager has spaceSizeX, spaceSizeY and spaceSizeZ
		// attributes. Nevertheless this attributes are retrieved from the 
		// network .ned file as network parameters.
	gates:
	    input in;
		output out;
}