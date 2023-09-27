/*                                              __      __
  ____  ____ ___  ___  ____ _____ _      ____ _/ /___  / /_  ____ _
 / __ \/ __ `__ \/ _ \/ __ `/ __ `/_____/ __ `/ / __ \/ __ \/ __ `/
/ /_/ / / / / / /  __/ /_/ / /_/ /_____/ /_/ / / /_/ / / / / /_/ /
\____/_/ /_/ /_/\___/\__, /\__,_/      \__,_/_/ .___/_/ /_/\__,_/
                    /____/                   /_/                  */

// Omega-alpha is a preliminary version of the Omega ocean model. The modeling
// capability is the same as Omega-0: the shallow water equations with
// identical vertical layers and inactive tracers. Omega-alpha runs on a single
// cpu with no domain decomposition, no GPU capability, and with native c++
// vectors rather than YAKL array objects. The purpose is to work out the
// design of the classes, functions, files, naming conventions, and code style
// in a very simple setting. The design of all Omega versions and the model
// equations may be found at:
// https://github.com/E3SM-Project/Omega/blob/develop/components/omega/doc/design/ShallowWaterOmega0.md
// Mesh generation instructions in Mesh.cpp
//
// Copyright (c) 2023, Mark Petersen
// All rights reserved.
// This source code is licensed under the BSD-style license found in the
// LICENSE file in the root directory of this source tree.

#include "Config.h"
#include "Diag.h"
#include "Mesh.h"
#include "Meta.h"
#include "State.h"
#include "Tend.h"
#include "timestep.h"
#include <iostream>
#include <string>

int main() {
  LOG(2, "Omega-Alpha initialization")

  //*******************************************************
  //  Initialization
  //*******************************************************
  Config config;
  Meta meta(config);
  Mesh m(config);
  std::vector<State> s;
  {
    State temporaryState(m);
    for (size_t n = 0; n < meta.stateLevelsInMemory; n++) {
      s.push_back(temporaryState);
    }
  } // temporaryState is destroyed here.
  Diag d(m);
  std::vector<Tend> tend;
  {
    Tend temporaryTend(m);
    for (size_t n = 0; n < meta.tendLevelsInMemory; n++) {
      tend.push_back(temporaryTend);
    }
  } // temporaryTend is destroyed here.

  s[0].init(config, m);
  d.compute(config, m, s[0]);
  if (config.output_on_startup) {
    // write_output(config, meta, m,s,d, tend);
  }

  //*******************************************************
  //  time step loop
  //*******************************************************
  for (meta.timeIndex = 0; meta.timeIndex < config.n_timesteps;
       meta.timeIndex++) {
    timestep(config, meta, m, s, d, tend);

    if (meta.timeIndex % config.output_frequency == 0) {
      // write_output(config, meta, m,s,d, tend);
    }
  }

  //*******************************************************
  //  output
  //*******************************************************
}
