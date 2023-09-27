#include "Meta.h"
#include "Config.h"
#include <iostream>
#include <string>
#include <vector>

Meta::Meta(Config &config) {
  LOG(4, "Meta constructor")

  if (config.timestep_method == "forward_Euler") {
    stateLevelsInMemory = 2;
    tendLevelsInMemory = 1;
  }

  timeArrayIndex.resize(stateLevelsInMemory); // index to large arrays in time
  for (int n = 0; n < stateLevelsInMemory; n++) {
    timeArrayIndex[n] = n;
  }
  LOG(4, "timeArrayIndex: " << timeArrayIndex[0] << timeArrayIndex[1])
}
