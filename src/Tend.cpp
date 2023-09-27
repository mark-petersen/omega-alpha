// Using netcdf-c version of the netcdf libraries, and all
// arrays are c++ vector containers, 1D with computed index offset

#include "Tend.h"
#include "Mesh.h"
#include "State.h"

// constructor
Tend::Tend(Mesh &m) : State(m) { LOG(4, "-> Tend::Tend") }
