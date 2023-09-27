// Using netcdf-c version of the netcdf libraries, and all
// arrays are c++ vector containers, 1D with computed index offset

#include "Mesh.h"
#include "Config.h"
#include "io.h"
#include <iostream>
#include <netcdf>
#include <string>

// constructor
Mesh::Mesh(Config &config) {
  LOG(4, "-> Mesh::Mesh")

  std::string meshFileName = config.dirName + config.fileName;

  int ncid, retval;
  if (config.verbose)
    std::cout << "** Opening file: " << meshFileName << " **" << std::endl;
  if ((retval = nc_open((meshFileName).c_str(), NC_NOWRITE, &ncid)))
    ERR(retval);
  if (config.verbose)
    std::cout << std::endl;

  LOG(4, "** Read in dimensions **")
  nCells = readNCDim(ncid, "nCells");
  nEdges = readNCDim(ncid, "nEdges");
  nVertices = readNCDim(ncid, "nVertices");
  maxEdges = readNCDim(ncid, "maxEdges");
  maxEdges2 = readNCDim(ncid, "maxEdges2");
  vertexDegree = readNCDim(ncid, "vertexDegree");
  // set nVertLevels, either from init file or from config.
  if (config.initial_condition == "init_file") {
    nVertLevels = readNCDim(ncid, "nVertLevels");
  } else {
    nVertLevels = config.initialize_nVertLevels;
  }
  // inner dimension shortcuts
  K = nVertLevels;
  ME = maxEdges;
  ME2 = maxEdges2;
  VD = vertexDegree;

  LOG(4, "** Read in mesh variables **")
  latCell = readNCDouble(ncid, "latCell", nCells);
  lonCell = readNCDouble(ncid, "lonCell", nCells);
  xCell = readNCDouble(ncid, "xCell", nCells);
  yCell = readNCDouble(ncid, "yCell", nCells);
  zCell = readNCDouble(ncid, "zCell", nCells);
  bottomDepth = readNCDouble(ncid, "bottomDepth", nCells);
  latEdge = readNCDouble(ncid, "latEdge", nEdges);
  lonEdge = readNCDouble(ncid, "lonEdge", nEdges);
  xEdge = readNCDouble(ncid, "xEdge", nEdges);
  yEdge = readNCDouble(ncid, "yEdge", nEdges);
  zEdge = readNCDouble(ncid, "zEdge", nEdges);
  latVertex = readNCDouble(ncid, "latVertex", nVertices);
  lonVertex = readNCDouble(ncid, "lonVertex", nVertices);
  xVertex = readNCDouble(ncid, "xVertex", nVertices);
  yVertex = readNCDouble(ncid, "yVertex", nVertices);
  zVertex = readNCDouble(ncid, "zVertex", nVertices);
  areaCell = readNCDouble(ncid, "areaCell", nCells);
  angleEdge = readNCDouble(ncid, "angleEdge", nEdges);
  dcEdge = readNCDouble(ncid, "dcEdge", nEdges);
  dvEdge = readNCDouble(ncid, "dvEdge", nEdges);
  areaTriangle = readNCDouble(ncid, "areaTriangle", nVertices);
  meshDensity = readNCDouble(ncid, "meshDensity", nCells);
  weightsOnEdge = readNCDouble(ncid, "weightsOnEdge", maxEdges2 * nEdges);
  // description="Reconstruction weights associated with each of the
  // edgesOnEdge,
  //   used to reconstruct the tangentialVelocity from normalVelocities on
  //   neighboring edges."

  // Cell pointers. These need to be reduced by 1 for Fortran->C
  indexToCellID = readNCInt(ncid, "indexToCellID", nCells, true);
  indexToEdgeID = readNCInt(ncid, "indexToEdgeID", nEdges, true);
  indexToVertexID = readNCInt(ncid, "indexToVertexID", nVertices, true);
  cellsOnCell = readNCInt(ncid, "cellsOnCell", nCells * maxEdges, true);
  edgesOnCell = readNCInt(ncid, "edgesOnCell", nCells * maxEdges, true);
  verticesOnCell = readNCInt(ncid, "verticesOnCell", nCells * maxEdges, true);
  edgesOnEdge = readNCInt(ncid, "edgesOnEdge", nEdges * maxEdges2, true);
  cellsOnEdge = readNCInt(ncid, "cellsOnEdge", nEdges * 2.0, true);
  verticesOnEdge = readNCInt(ncid, "verticesOnEdge", nEdges * 2.0, true);
  cellsOnVertex =
      readNCInt(ncid, "cellsOnVertex", nVertices * vertexDegree, true);
  edgesOnVertex =
      readNCInt(ncid, "edgesOnVertex", nVertices * vertexDegree, true);

  // These do not need to be reduced by 1 for Fortran->C
  nEdgesOnCell = readNCInt(ncid, "nEdgesOnCell", nCells, false);
  nEdgesOnEdge = readNCInt(ncid, "nEdgesOnEdge", nEdges, false);

  // mesh quality variables. Not in init file.
  // cellQuality = readNCDouble(ncid, "cellQuality", nCells);
  // gridSpacing = readNCDouble(ncid, "gridSpacing", nCells);
  // triangleQuality = readNCDouble(ncid, "triangleQuality", nCells);
  // triangleAngleQuality = readNCDouble(ncid, "triangleAngleQuality", nCells);
  // boundaryVertex = readNCInt(ncid, "boundaryVertex", nVertices);
  // obtuseTriangle = readNCInt(ncid, "obtuseTriangle", nCells);

  LOG(4, "** Closing file: " << meshFileName << " **")
  if ((retval = nc_close(ncid)))
    ERR(retval);
  if (config.verbose)
    std::cout << std::endl;

  {
    edgeSignOnCell.resize(nCells * maxEdges);
    size_t iCell, i, cell1, cell2, iEdge;
    for (iCell = 0; iCell < nCells; iCell++) {
      for (i = 0; i < nEdgesOnCell[iCell]; i++) {
        iEdge = edgesOnCell[iCell * ME + i];
        cell1 = cellsOnEdge[iEdge * 2];
        cell2 = cellsOnEdge[iEdge * 2 + 1];
        // Vectors point from lower to higher cell number.
        // If my cell number is higher than my neighbor, then
        // vector points towards me and the edge sign is positive.
        if (iCell == std::max(cell1, cell2)) {
          edgeSignOnCell[iCell * ME + i] = 1;
        } else {
          edgeSignOnCell[iCell * ME + i] = -1;
        }
      }
    }
  }

  {
    edgeSignOnVertex.resize(nVertices * VD);
    size_t iVertex, i, iEdge;
    for (iVertex = 0; iVertex < nVertices; iVertex++) {
      for (i = 0; i < VD; i++) {
        iEdge = edgesOnVertex[iVertex * VD + i];
        // Vectors point from lower to higher vertex number.
        // If my vertex number is higher than my neighbor, then
        // vector points towards me and the edge sign is positive.
        if (iVertex == verticesOnEdge[iEdge * 2]) {
          edgeSignOnVertex[iVertex * VD + i] = -1;
        } else {
          edgeSignOnVertex[iVertex * VD + i] = 1;
        }
      }
    }
  }
}
// mrp add printing routing later
//  std::cout << "IC mesh xCell: ";
//  for (size_t i=0; i<16; i++) {
//    std::cout << xCell[i]<< " ";
//  }
//  std::cout << std::endl;
