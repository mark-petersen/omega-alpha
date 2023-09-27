#include "timestep.h"
#include "Config.h"
#include "Diag.h"
#include "Mesh.h"
#include "Meta.h"
#include "State.h"
#include "Tend.h"
#include "io.h"
#include "tendencies.h"
#include <cmath> // for exp
#include <iostream>
#include <string>
#include <vector>

void forward_Euler_timestep(Config &config, Meta &meta, Mesh &m,
                            std::vector<State> &s, Diag &d,
                            std::vector<Tend> &tend) {
  LOG(4, "-> forward_Euler_timestep")
  // y(n+1) = y(n) + dt*RHS(n)

  size_t tCur = meta.timeArrayIndex[0];
  size_t tNew = meta.timeArrayIndex[1];
  compute_velocity_tendencies(config, meta, m, s[tCur], d, tend[0]);
  compute_thickness_tendencies(config, meta, m, s[tCur], d, tend[0]);

  size_t K = m.nVertLevels;
  for (size_t e = 0; e < m.nEdges; e++) {
    for (size_t k = 0; k < K; k++) {
      s[tNew].normalVelocity[e * K + k] =
          s[tCur].normalVelocity[e * K + k] +
          config.dt * tend[0].normalVelocity[e * K + k];
    }
  }
  for (size_t i = 0; i < m.nCells; i++) {
    for (size_t k = 0; k < K; k++) {
      s[tNew].layerThickness[i * K + k] =
          s[tCur].layerThickness[i * K + k] +
          config.dt * tend[0].layerThickness[i * K + k];
    }
  }
  d.compute(config, m, s[tNew]);

  // switch array indices for next step:
  meta.timeArrayIndex[0] = tNew;
  meta.timeArrayIndex[1] = tCur;
}

void timestep(Config &config, Meta &meta, Mesh &m, std::vector<State> &s,
              Diag &d, std::vector<Tend> &tend) {
  LOG(4, "-> timestep")

  if (config.timestep_method == "forward_Euler") {
    forward_Euler_timestep(config, meta, m, s, d, tend);
  }

  // check against exact solution
  size_t tOld = meta.timeArrayIndex[1];
  size_t tNew = meta.timeArrayIndex[0];
  double curTime = (meta.timeIndex) * config.dt;
  double sol = config.initial_condition_amplitude *
               std::exp(-config.hTend_del2_coef * curTime *
                        pow(2 * M_PI / config.Lx, 2));
  printXRow(s[tOld].layerThickness, 16, "comp  ");
  //    std::cout<<s[tNew].layerThickness[9]<<"
  //    "<<s[tNew].layerThickness[12]<<std::endl;
  //    for (size_t i=0; i<16; i++) {
  //        std::cout << i<<s[tNew].layerThickness[i]<< "u";
  //    }
  //    std::cout<<std::endl;
  printf("Sol:                                                 %.8f \n", sol);
}
