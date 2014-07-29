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

#ifndef PARTICLE_H
#define PARTICLE_H

#include "Defines.h"
#include <vector>

class Particle {

	private:

	protected:

		bool active_;

		point_t position_;

		vect_t velocity_;

		double mass_;

		double lastCollisionTime_;

		int particleId_;

		int particleType_;

		// Cell List attributes
		int spaceCell_; // TODO remove this and start using spaceCellIdx instead

		index_t spaceCellIdx_;

		int prevSpaceCell_;

		index_t prevSpaceCellIdx_;

		// Near Neighbor List of particles
		std::vector<Particle*> neighborParticles_;

		// Near Neighbor List radius
		double listRadius_;

		// Near Neighbor List refresh list radius
		double refreshListRadius_;

		// Boundaries mode. Whether the particle must bounce, expire, etc.
		int boundariesMode_;

		//
		// Brownian motion parameters
		// 
		double diffusion_;

		double inertia_;

		double viscosity_;

		double BMStdDev_;

	public:

		Particle();

		Particle(double, double, double, double, double);

		// Gets and sets
		inline int getParticleId(void) { return particleId_; };

		inline int getParticleType(void) { return particleType_; };

		inline double getX(void) { return position_.x; };

		inline double getY(void) { return position_.y; };

		inline double getZ(void) { return position_.z; };

		inline point_t *getPosition(void) { return &position_; };

		inline double getVx(void) { return velocity_.x; };

		inline double getVy(void) { return velocity_.y; };

		inline double getVz(void) { return velocity_.z; };

		inline vect_t *getVelocity(void) { return &velocity_; };

		inline double getMass(void) { return mass_; };

		inline double *getDiffusion(void) { return &diffusion_; };

		inline double *getInertia(void) { return &inertia_; };

		inline double *getViscosity(void) { return &viscosity_; };

		inline double *getBMStdDev(void) { return &BMStdDev_; };

		// TODO remove this once the space cells indexes use data structures (spaceCellIdx)
		inline int getSpaceCell(void) { return spaceCell_; };

		inline index_t getSpaceCellIdx(void) { return spaceCellIdx_; };

		// TODO remove this once the space cells indexes use data structures (spaceCellIdx)
		inline int getPrevSpaceCell(void) { return prevSpaceCell_; };

		inline index_t getPrevSpaceCellIdx(void) { return prevSpaceCellIdx_; };

		inline double getListRadius(void) { return listRadius_; };

		inline double getRefreshListRadius(void) { return refreshListRadius_; };

		virtual double getRadius(void) = 0;

		inline double getLastCollisionTime(void) { return lastCollisionTime_; };

		inline std::vector<Particle*> getNeighborParticles(void) { return neighborParticles_; };

		inline int getBoundariesMode(void) const { return boundariesMode_; };

		inline bool isActive (void) { return active_; };
		
		inline void setParticleId(int id) { particleId_ = id; };

		inline void setParticleType(int t) { particleType_ = t; };

		inline void setX(double x) { position_.x = x; };

		inline void setY(double y) { position_.y = y; };

		inline void setZ(double z) { position_.z = z; };

		inline void setPosition(point_t p) { position_ = p; };

		inline void setVx(double vx) { velocity_.x = vx; };

		inline void setVy(double vy) { velocity_.y = vy; };

		inline void setVz(double vz) { velocity_.z = vz; };

		inline void setVelocity(vect_t v) { velocity_ = v; };

		inline void setMass(double m) { mass_ = m; };

		inline void setDiffusion(double d) { diffusion_ = d; };

		inline void setInertia(double i) { inertia_ = i; };

		inline void setViscosity(double v) { viscosity_ = v; };

		inline void setBMStdDev(double sd) { BMStdDev_ = sd; };

		inline void setSpaceCell(int c) { spaceCell_ = c; };

		inline void setSpaceCellIdx(index_t sci) { spaceCellIdx_ = sci; };

		inline void setPrevSpaceCell(int c) { prevSpaceCell_ = c; };

		inline void setPrevSpaceCellIdx(index_t psci) { prevSpaceCellIdx_ = psci; };

		inline void setListRadius(double r) { listRadius_ = r; };

		inline void setRefreshListRadius(double r) { refreshListRadius_ = r; };

		virtual void setRadius(double) = 0;

		inline void setLastCollisionTime(double tc) { lastCollisionTime_ = tc; };

		inline void setBoundariesMode(int bm) { boundariesMode_ = bm; };

		inline void setActive(bool a) { active_ = a; };

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
