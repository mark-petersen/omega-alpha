#ifndef TIMESTEP_H
#define TIMESTEP_H

#include "Config.h"
#include "Diag.h"
#include "Mesh.h"
#include "Meta.h"
#include "State.h"
#include "Tend.h"
#include "vector"

void timestep(Config &config, Meta &meta, Mesh &m, std::vector<State> &s,
              Diag &diag, std::vector<Tend> &tend);

#endif
