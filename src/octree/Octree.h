#include "Defines.h"
#include "Particle.h"

#include <vector>
#include <list>

class Octree {

  typedef std::vector<std::list<Particle*> > layer_t;

 private:

 protected:

  unsigned int numLayers_;

  double maxSideLength_;

  vect_t spaceSize_;

  std::list<Particle*> particles_;

  std::vector<layer_t> spaceCellsTree_;

 public:

  Octree();

  Octree(unsigned int nl, double msl);

  ~Octree();

  unsigned int getNumLayers(void) { return numLayers_; };

  double getMaxSideLength(void) { return maxSideLength_; };

  vect_t getSpaceSize(void) { return spaceSize_; };

  double getSpaceCellSideLength(index_t);

  void setNumLayers(unsigned int nl) { numLayers_ = nl; };

  void setMaxSideLength(double sl) { maxSideLength_ = sl; };

  void setSpaceSize(vect_t ss) { spaceSize_ = ss; }; 

  std::list<Particle *> *getSpaceCellParticles(index_t idx);

  std::list<Particle *> getNeighborParticles(index_t idx);

  void transferParticle(Particle *p, index_t from, index_t to);

  void subscribe(Particle *p);

  void unsubscribe(Particle *p);

  void attachParticleToSpaceCell(Particle *, index_t idx);

  void detachParticleFromSpaceCell(Particle *, index_t idx);

};