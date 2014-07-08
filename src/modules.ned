//  This file is part of the Cell-Signaling project. Cell-Signaling is an
//  Omnet++ project to simulate cell signaling communications.
//  Copyright (C) 2014  Daniel Huertas
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

package cellsignaling.src;

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

		// Space Cell Size. Sets the size of the space cells in which the
		// simulation space will be divided. Note that the space cell size must
		// be at least as wide as the widest diameter of the particles it 
		// contains. If set to 0, the manager will adjust the space cell size
		// by searching the maximum radius of the particles (2*MaxRadius).
		double spaceCellSize = default(0);

		// Verlet list radius. When set, overwrites the list radius of the 
		// particles when using mode 2 (near neighbor lists, or verlet lists).
		double listRadius = default(0);

		// Place the particles following a specific pattern. When left unset 
		// particles get their position from the ini file. Options are:
		// - undefined
		// - uniform
		// - cube
		// - sphere
		// - twosided
		// - highdensity
		string particleDistribution = default("undefined");

		// Delta time for Brownian Motion
		double deltaTime = default(0);

		// Tk refresh rate. Tells the manager when to update the position
		// string of the particles for the tk environment.
		double tkRefreshRate = default(5000);

		// Statistics refresh rate.
		double statsRefreshRate = default(0);

		// Enable the Manager Web Server in order to visualize a 3D simulation 
		// space from a web browser (WARNING: enable this option only if your
		// browser supports WebGL).
		int enableWebServer = default(0);

	gates:
		input in;
		output out;
}

simple Molecule {

	parameters:

		// The position of the particle at the time of the last event
		double xpos;			// x position
		double ypos;			// y position
		double zpos;			// z position

		// The velocity vector of the particle
		double vx = default(0);	// x component velocity
		double vy = default(0);	// y component velocity
		double vz = default(0);	// z component velocity

		// Radius of the particle
		double radius = default(1);

		// Mass of the particle
		double mass = default(1);

		// Diffusion
		double diffusion = default(0);

		// Inertia
		double inertia = default(0);

		// Viscosity 
		double viscosity = default(0);

		// Radius for the near neighbor list when using the Verlet list
		// algorithm (mode = M_NNLIST).
		double listRadius = default(0);

		// Radius used to update the neighbor list of the molecule
		double refreshListRadius = default(0);

		// Amount of time that the particle exists in the simulation space. An
		// event (EV_TTLEXPIRED) is used to manage when a particle must leave.
		// Default 0 (does not expire).
		double timeToLive = default(0);

		// Particle behavior upon collision with a boundary element (cube walls
		// for now). Possible modes are:
		// 1 - Elastic wall collision (default). The particle bounces off the
		//     wall.
		// 2 - Particle expires. When the particle reaches a boundary it leaves
		//     the simulation space.
		// 3 - Periodic boundaries. The particle leaves the simulation space 
		//     and appears on the opposite boundary, with the same velocity
		//     vector (not implemented yet). Instead of the particle surface,
		//     the center of the particle is used to detect the boundary
		//     transfer.
		int boundariesMode = default(1);

		// Statistics refresh rate
		double statsRefreshRate = default(0);

		// Log collisions
		int logCollisions = default(0);

		// Log position
		int logPosition = default(0);

	gates:
		input in;
		output out;
}

simple MoleculeEmitter {
	
	parameters:

		// Whether the cell will emit molecules or not
		bool enabled = default(false);

		// Choose the way molecules must be emitted
		// TODO implement functions
		string emissionFunction = default("pulse");

		// The time to start emitting molecules
		double emissionStart = default(-1);

		// How long the emitter releases molecules
		double emissionDuration = default(0);

		// How many molecules are released per simulation time unit
		double emissionRate = default(0);

		// Radius of the emitted molecules
		double emissionParticleRadius = default(1);

		// Mass of the emitted molecules
		double emissionParticleMass = default(1);

		// How long the molecules must remain in the simulation before expiring
		double emissionTimeToLive = default(0);

		// Whether emitted molecules must leave the simulation space or not
		int emissionBoundariesMode = default(1);

		// Molecule initial velocity
		double emissionVelocity = default(1);

		// Verlet list radius parameter for the emitted molecules
		double emissionListRadius = default(0);

		// Verlet refresh list radius for the emitted molecules
		double emissionRefreshListRadius = default(0);

		// Diffusion coefficient of molecules when using brownian motion
		// Note: manager module must have a deltaTime value greater than 0
		// TODO check that manager.deltaTime is set during initialization process
		double emissionDiffusion = default(0);

		// Create molecule modules during the initialization process. This is 
		// useful when the emitter must release a lot of molecules.
		// Note: This uses more memory at the beginning of the simulation, but
		// speeds up the simulation execution.
		bool preloadMolecules = default(false);

	gates:
		input in;
		output out;

}

simple MoleculeReceiver {

	parameters:

		bool enabled = default(false);

		// Statistics refresh rate
		double statsRefreshRate = default(0);

	gates:
		input in;
		output out;

}

simple SimpleCell {

	parameters:
		//
		// Mobility parameters
		//
		double xpos;
		double ypos;
		double zpos;

		double vx = default(0);
		double vy = default(0);
		double vz = default(0);

		// Radius of the cell
		double radius = default(1);

		// Mass of the cell
		double mass = default(1);

		// Diffusion
		double diffusion = default(0);

		// Inertia
		double inertia = default(0);

		// Viscosity 
		double viscosity = default(0);

		// Radius for the near neighbor list when using the Verlet list
		// algorithm (mode = M_NNLIST).
		double listRadius = default(0);

		double refreshListRadius = default(0);

		// Amount of time that the cell lives. Default 0 (does not expire).
		double timeToLive = default(0);

		// Same as the Molecule boundariesMode parameter.
		int boundariesMode = default(1);

		// Statistics refresh rate in, 1/1000 of a time unit
		double statsRefreshRate = default(0);

		// Log collision times
		int logCollisions = default(0);

	gates:
		input in;
		output out;

}

module Cell {

	gates:
		input in;
		output out;

	submodules:

		mobility: SimpleCell;
		emitter: MoleculeEmitter;
		receiver: MoleculeReceiver;

	connections allowunconnected:

}

simple Probe {
	
	parameters:
		// Probe name
		string name = default("probe");

		int type = default(6); // T_SIGNALING

		double xpos;
		double ypos;
		double zpos;

		double radius = default(1);

		// Statistics refresh rate in, 1/1000 of a time unit
		double statsRefreshRate = default(0);

	gates:
		input in;
		output out;

}
