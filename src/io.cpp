// Using netcdf-c version of the netcdf libraries, and all
// arrays are c++ vector containers, 1D with computed index offset

#include "io.h"
#include "Config.h"
#include <iostream>
#include <netcdf>
#include <string>
#include <vector>

void printXRow(std::vector<double> &var, size_t Nx, std::string text) {
  std::cout << text << ": ";
  for (size_t i = 0; i < 16; i++) {
    printf(" %.8f", var[i]);
  }
  printf("\n");
}

std::vector<int> readNCInt(int ncid, std::string varName, size_t dim,
                           bool minus1) {
  int varid, retval;
  std::vector<int> var;
  var.resize(dim);
  LOG(5, varName << ": ");
  if ((retval = nc_inq_varid(ncid, varName.c_str(), &varid)))
    ERR(retval);
  if ((retval = nc_get_var_int(ncid, varid, &var[0])))
    ERR(retval);
  if (minus1) {
    for (size_t i = 0; i < dim; i++) {
      var[i]--;
    }
  }
  LOG(5, var[0] << ", " << var[dim - 1])
  return var;
}

std::vector<double> readNCDouble(int ncid, std::string varName, size_t dim) {
  int varid, retval;
  std::vector<double> var;
  var.resize(dim);
  LOG(5, varName << ": ");
  if ((retval = nc_inq_varid(ncid, varName.c_str(), &varid)))
    ERR(retval);
  if ((retval = nc_get_var_double(ncid, varid, &var[0])))
    ERR(retval);
  LOG(5, var[0] << ", " << var[dim - 1])
  return var;
}

void fillNCDouble(int ncid, std::string varName, std::vector<double> &var) {
  int varid, retval;
  LOG(5, varName << ": ");
  if ((retval = nc_inq_varid(ncid, varName.c_str(), &varid)))
    ERR(retval);
  if ((retval = nc_get_var_double(ncid, varid, &var[0])))
    ERR(retval);
  LOG(5, var[0] << ", " << var[1])
}

size_t readNCDim(int ncid, std::string dimName) {
  int dimid, retval;
  size_t dim;
  LOG(5, dimName << ": ");
  if ((retval = nc_inq_dimid(ncid, dimName.c_str(), &dimid)))
    ERR(retval);
  if ((retval = nc_inq_dimlen(ncid, dimid, &dim)))
    ERR(retval);
  LOG(5, dim)
  return dim;
}
