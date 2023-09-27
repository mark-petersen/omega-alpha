#ifndef TEND_H
#define TEND_H

#include "Mesh.h"
#include "Meta.h"
#include "State.h"
#include <string>
#include <vector>

// Tend is derived from the State class, so it
// inherits all its variables.
class Tend : public State {
public:
  // tend variables: nothing extra!

  // constructor
  Tend(Mesh &m);
};
#endif
