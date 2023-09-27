#ifndef CONFIG_H
#define CONFIG_H

/* VERBOSITY LEVELS
0. No printing
1. warnings only
2. High level code locations - init, timestep, finalize only (like mpas log).
3. Mid level information - summary diagnostics on simulation in text output per
time step.
4. Print upon entry to most subroutines
5. Max printing - individual variable reading
*/

#define VERBOSITY 3
#define LOG(verbose_level, message)                                            \
  {                                                                            \
    if (VERBOSITY >= verbose_level)                                            \
      std::cout << message << std::endl;                                       \
  }
#define ERRORMESSAGE(m)                                                        \
  {                                                                            \
    printf("Error: %s\n", m);                                                  \
    exit(0);                                                                   \
  }

#include <string>

inline constexpr double gravity =
    9.807; // m/s^2 -- gravity is light on this planet

class Config {
public:
  //*******************************************************
  //    Initialization
  //*******************************************************
  bool verbose{true}; // print verbose output during run
  std::string dirName{"links_su/"};
  // std::string fileName { "ocean.QU.240km.151209.nc" };
  // std::string fileName { "mpaso.EC30to60E2r3.230313.nc" };
  // std::string fileName { "oRRS18to6v3.171116.nc" };
  std::string fileName{"mpas_mesh_16x16.nc"};

  // std::string initial_condition = "init_file";
  size_t initialize_nVertLevels =
      1; // nVertLevels if initial_condition!="init_file"
  // std::string initial_condition = "constant";
  double initial_condition_amplitude = 1.0;
  std::string initial_condition = "sinx";
  double Lx = 64.0e3 * 16.0; // change later

  //*******************************************************
  //    time stepping
  //*******************************************************
  std::string timestep_method = "forward_Euler";
  double dt = 0.001; // time step [s]
  size_t n_timesteps = 10;

  //*******************************************************
  //    output
  //*******************************************************
  bool output_on_startup = true;
  size_t output_frequency = 1; // time step interval for output

  //*******************************************************
  //    velocity tendency terms
  //*******************************************************
  bool uTend_advection_enable = true;

  bool uTend_ssh_gradient_enable = true;
  bool uTend_KE_gradient_enable = true;

  double coriolis =
      1e-4; // Coriolis parameter f (to do: change to variable later)

  bool uTend_del2_enable = true;
  double uTend_del2_coef = 1e-4; // m^2/s
  bool uTend_del4_enable = true;
  bool uTend_bottom_drag_enable = true;

  bool uTend_Rayleigh_enable = true;
  double uTend_Rayleigh_drag = 0.1; // coefficient in -Ra*u term [1/s]

  bool uTend_wind_forcing_enable = true;

  //*******************************************************
  //    thickness tendency terms
  //*******************************************************
  bool hTend_advection_enable = true;

  bool hTend_decay_enable = false; // -c*h ( for testing only)
  double hTend_decay_coef = 0.1;   // coefficient Ra [1/s]

  bool hTend_del2_enable = true;  // kappa del2(h)
  double hTend_del2_coef = 1.0e7; // coefficient in -Ra*u term [1/s]
};
#endif
