#ifndef META_H
#define META_H

#include "Config.h"
#include <iostream>
#include <string>
#include <vector>

class Meta {
public:
  size_t timeIndex = 0;               // index for time stepping
  std::vector<size_t> timeArrayIndex; // index to large arrays in time
  double time = 0.0;                  // model time in seconds
  size_t stateLevelsInMemory = 2;
  size_t tendLevelsInMemory = 1;

  // for later:
  // string xtime; // type="text" dimensions="Time" description="model time,
  // with format 'YYYY-MM-DD_HH:MM:SS'" string simulationStartTime; //
  // type="text" dimensions="" default_value="'no_date_available'"
  // description="start time of first simulation, with format
  // 'YYYY-MM-DD_HH:MM:SS'" double daysSinceStartOfSim; // type="real"
  // dimensions="Time" units="days" description="Time since simulationStartTime,
  // for plotting"

  // constructor
  Meta(Config &config);
};

#endif
