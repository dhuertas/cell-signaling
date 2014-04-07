#ifndef PARTICLE_H
#define PARTICLE_H

#include "Defines.h"
#include <vector>

class Particle {

	private:

	protected:

		bool active;

		point_t position;

		vect_t velocity;

		double mass;

		double lastCollisionTime;

		int particleId;

		int particleType;

		// Cell List attributes
		int spaceCell; // TODO remove this and start using spaceCellIdx instead

		index_t spaceCellIdx;

		int prevSpaceCell;

		index_t prevSpaceCellIdx;

		// Near Neighbor List of particles
		std::vector<Particle*> neighborParticles;

		// Near Neighbor List radius
		double listRadius;

		// Near Neighbor List refresh list radius
		double refreshListRadius;

		// Boundaries mode. Whether the particle must bounce, expire, etc.
		int boundariesMode;

		//
		// Brownian motion parameters
		// 
		double diffusion;

		double inertia;

		double viscosity;

		double BMStdDev;

	public:

		Particle();

		Particle(double, double, double, double, double);

		// Gets and sets
		inline int getParticleId(void) { return particleId; };

		inline int getParticleType(void) { return particleType; };

		inline double getX(void) { return position.x; };

		inline double getY(void) { return position.y; };

		inline double getZ(void) { return position.z; };

		inline point_t *getPosition(void) { return &position; };

		inline double getVx(void) { return velocity.x; };

		inline double getVy(void) { return velocity.y; };

		inline double getVz(void) { return velocity.z; };

		inline vect_t *getVelocity(void) { return &velocity; };

		inline double getMass(void) { return mass; };

		inline double *getDiffusion(void) { return &diffusion; };

		inline double *getInertia(void) { return &inertia; };

		inline double *getViscosity(void) { return &viscosity; };

		inline double *getBMStdDev(void) { return &BMStdDev; };

		// TODO remove this once the space cells indexes use data structures (spaceCellIdx)
		inline int getSpaceCell(void) { return spaceCell; };

		inline index_t getSpaceCellIdx(void) { return spaceCellIdx; };

		// TODO remove this once the space cells indexes use data structures (spaceCellIdx)
		inline int getPrevSpaceCell(void) { return prevSpaceCell; };

		inline index_t getPrevSpaceCellIdx(void) {return prevSpaceCellIdx; };

		inline double getListRadius(void) { return listRadius; };

		inline double getRefreshListRadius(void) { return refreshListRadius; };

		virtual double getRadius(void) = 0;

		inline double getLastCollisionTime(void) { return lastCollisionTime; };

		inline std::vector<Particle*> getNeighborParticles(void) { return neighborParticles; };

		inline int getBoundariesMode(void) const { return boundariesMode; };

		inline bool isActive (void) { return active; };
		
		inline void setParticleId(int id) { particleId = id; };

		inline void setParticleType(int t) { particleType = t; };

		inline void setX(double x) { position.x = x; };

		inline void setY(double y) { position.y = y; };

		inline void setZ(double z) { position.z = z; };

		inline void setPosition(point_t p) { position = p; };

		inline void setVx(double vx) { velocity.x = vx; };

		inline void setVy(double vy) { velocity.y = vy; };

		inline void setVz(double vz) { velocity.z = vz; };

		inline void setVelocity(vect_t v) { velocity = v; };

		inline void setMass(double m) { mass = m; };

		inline void setDiffusion(double d) { diffusion = d; };

		inline void setInertia(double i) { inertia = i; };

		inline void setViscosity(double v) { viscosity = v; };

		inline void setBMStdDev(double sd) { BMStdDev = sd; };

		inline void setSpaceCell(int c) { spaceCell = c; };

		inline void setSpaceCellIdx(index_t sci) { spaceCellIdx = sci; };

		inline void setPrevSpaceCell(int c) { prevSpaceCell = c; };

		inline void setListRadius(double r) { listRadius = r; };

		inline void setRefreshListRadius(double r) { refreshListRadius = r; };

		virtual void setRadius(double) = 0;

		inline void setLastCollisionTime(double tc) { lastCollisionTime = tc; };

		inline void setBoundariesMode(int bm) { boundariesMode = bm; };

		inline void setActive(bool a) { active = a; };

		//
		// tk Environment methods
		//
		virtual void tkEnvDrawShape(void) = 0;

		virtual void tkEnvUpdatePosition(void) = 0;

		virtual void tkEnvUpdatePosition(double) = 0;

		//
		// Event related methods
		//
		virtual void initMobilityMessages() = 0;

		virtual void initializeMobility() = 0;

		virtual void finishMobility() = 0;

		virtual void createNearNeighborList() = 0;

		virtual void updateNearNeighborList() = 0;

		virtual void expire() = 0;
};

#endif
