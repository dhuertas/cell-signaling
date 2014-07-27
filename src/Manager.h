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

#ifndef MANAGER_H
#define MANAGER_H

#include <cmessage.h>
#include <csimplemodule.h>
#include <cqueue.h>
#include <omnetpp.h>

#include <cmath>
#include <vector>
#include <list>

#include "Particle.h"
#include "ParticleDistributionHelper.h"
#include "WebServer.h"

class Manager : public cSimpleModule {

  private:

    // Molecule Dynamics Mode of operation
    int mode_;

    // The number of space cells in each direction. They are used to
    // access the spaceCells vector of vectors.
    int Nx_;
    int Ny_;
    int Nz_;

    // File Descriptors to communicate with the Web Server thread
    int quitFd_[2];

    // the simulation space size in each direction
    vect_t spaceSize_;

    // Delta time
    double deltaTime_;

    // The space cell size
    double spaceCellSize_;

    // Overwrite particles' list radius.
    double listRadius_;

    // Contains the particle Id from the last added particle to the domain.
    int lastParticleId_;

    // TK environment refresh rate
    // It allows the manager module to send self-messages in order to
    // update the position of each particle.
    double tkEnvRefreshRate_;

    // Update the cOutVectors periodically.
    double statsRefreshRate_;

    int enableWebServer_;

    pthread_t webServerThread_;

    struct arg_struct webServerArgs_;

    // Manager name
    std::string name_;

    // A list of particles contained in the simulation space. Every new
    // particle must be subscribed to (and unsubscribed from).
    std::list<Particle*> particles_;

    // Space is divided into cells, each of which contains a list of
    // particles belonging to it.
    std::vector<std::list<Particle*> > spaceCells_;

  protected:

    statistics_t stats_;

    cOutVector allCollisionsVector_;

    cOutVector particleCollisionsVector_;

    cOutVector wallCollisionsVector_;

    cOutVector transfersVector_;

    cOutVector expiresVector_;

  public:

    ~Manager();

    // Every particle must be subscribed in order to access its attributes
    // during simulation time.
    void subscribe(Particle *);

    // Unsubcribe particles. Particles may expire, be received, leave the 
    // area, etc.
    void unsubscribe(Particle *);

    void attachParticleToSpaceCell(Particle *, int);

    void detachParticleFromSpaceCell(Particle *, int);

    // Move one particle from one space cell to another.
    void transferParticle(Particle *, int, int);

    // Update the tk environment
    void tkEnvUpdateNetwork(void);

    //
    // Web Server related methods
    //
    void startWebServerThread(void);

    void stopWebServerThread(void);

    // Gets and sets
    vect_t *getSpaceSize(void) { return &spaceSize_; };

    double getSpaceSizeX(void) { return spaceSize_.x; };

    double getSpaceSizeY(void) { return spaceSize_.y; };

    double getSpaceSizeZ(void) { return spaceSize_.z; };

    double getSpaceCellSize(void) { return spaceCellSize_; };

    int *getNumberOfSpaceCellsX(void) { return &Nx_; };

    int *getNumberOfSpaceCellsY(void) { return &Ny_; };

    int *getNumberOfSpaceCellsZ(void) { return &Nz_; };

    int getMode(void) { return mode_; };

    double getDeltaTime(void) { return deltaTime_; };

    double getListRadius(void) { return listRadius_; };

    std::list<Particle *> *getSpaceCellParticles(int);

    int getNextParticleId(void);

    int getLastParticleId(void) { return lastParticleId_; };

    void setSpaceSize(vect_t vsz) { spaceSize_ = vsz; };

    void setSpaceSizeX(double sx) { spaceSize_.x = sx; };

    void setSpaceSizeY(double sy) { spaceSize_.y = sy; };

    void setSpaceSizeZ(double sz) { spaceSize_.z = sz; };

    void setSpaceCellSize(double sc) { spaceCellSize_ = sc; };

    void setNumberOfSpaceCellsX(int n) { Nx_ = n; };

    void setNumberOfSpaceCellsY(int n) { Ny_ = n; };

    void setNumberOfSpaceCellsZ(int n) { Nz_ = n; };

    void setMode(int m) { mode_ = m; };

    void setDeltaTime(double dt) { deltaTime_ = dt; };

    void setListRadius(double lr) { listRadius_ = lr; };

    void clearStatistics(void);

    void registerCollision(void);

    void registerWallCollision(void);

    void registerTransfer(void);

    void registerExpire(void);

    //
    // cSimpleModule inheritance
    //
    virtual void initialize(int stage);

    virtual int numInitStages() const;

    virtual void handleMessage(cMessage *);

    virtual void finish();
};

#endif
