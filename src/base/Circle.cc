#include "Circle.h"

using namespace std;

Circle::Circle(
    double x,
    double y,
    double vx,
    double vy,
    double r,
    double mass) : Particle(x, y, vx, vy, mass) {

    radius = r;

}
