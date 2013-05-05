#ifndef MANAGER_H
#define MANAGER_H

#include <cmessage.h>
#include <csimplemodule.h>
#include <cqueue.h>
#include <omnetpp.h>

#include <vector>
#include <list>

#include "Particle.h"

class Manager : public cSimpleModule {

	private:

		// A list of particles contained in the simulation space. Every new
		// particle must be subscribed to (and unsubscribed from).
		std::list<Particle*> particles;

		// Space is divided into cells, each of which contains a list of
		// particles belonging to it.
		std::vector<std::list<Particle*> > spaceCells;

		// the simulation space size in each direction
		double spaceSizeX;
		double spaceSizeY;
		double spaceSizeZ;

		// The number of space cells in each direction. They are used to
		// access the spaceCells vector of vectors.
		int Nx; 
		int Ny;
		int Nz;

		// The space cell size
		double spaceCellSize;

		// TK environment refresh rate
		// It allows the manager module to send self-messages in order to
		// update the position of each particle.
		double tkEnvRefreshRate;

	protected:

	public:

		~Manager();

		// Every particle must be subscribed in order to access its attributes
		// during simulation time.
		void subscribe(Particle *);

		// Unsubcribe particles. Particles may expire, be received, leave the 
		// area, etc.
		void unsubscribe(Particle *);

		void attachParticleToSpaceCell(Particle *);
		void detachParticleFromSpaceCell(Particle *);

		// Move one particle from one space cell to another.
		void transferParticle(Particle *);

		// cSimpleModule inheritance
		virtual void initialize(int stage);
		virtual int numInitStages() const;
		virtual void handleMessage(cMessage *);
		virtual void finish();

		// Update the tk environment
		void tkEnvUpdateNetwork(void);

		// Gets and sets
		double getSpaceSizeX(void) { return spaceSizeX; };
		double getSpaceSizeY(void) { return spaceSizeY; };
		double getSpaceSizeZ(void) { return spaceSizeZ; };
		double getSpaceCellSize(void) { return spaceCellSize; };
		int getNumberOfSpaceCellsX(void) { return Nx; };
		int getNumberOfSpaceCellsY(void) { return Ny; };
		int getNumberOfSpaceCellsZ(void) { return Nz; };

		std::list<Particle *> getSpaceCellParticles(int);

		void setSpaceSizeX(double sx) { spaceSizeX = sx; };
		void setSpaceSizeY(double sy) { spaceSizeY = sy; };
		void setSpaceSizeZ(double sz) { spaceSizeZ = sz; };
		void setSpaceCellSize(double sc) { spaceCellSize = sc; };
		void setNumberOfSpaceCellsX(int n) { Nx = n; };
		void setNumberOfSpaceCellsY(int n) { Ny = n; };
		void setNumberOfSpaceCellsZ(int n) { Nz = n; };

};

#endif
