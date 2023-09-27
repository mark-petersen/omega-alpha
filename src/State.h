#ifndef STATE_H
#define STATE_H

#include "Config.h"
#include "Mesh.h"
#include "Meta.h"
#include <string>
#include <vector>

class State {
public:
  // State variables
  // could have add and remove functions for each variable
  std::vector<double> normalVelocity;
  std::vector<double> layerThickness;

  // constructor
  State(Mesh &m);
  void init(Config &config, Mesh &m);
};
#endif
