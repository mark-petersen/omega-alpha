#include "Diag.h"
#include "Config.h"
#include "Mesh.h"
#include "State.h"
#include <cmath> // for pow
#include <vector>

// constructor
Diag::Diag(Mesh &m) {
  LOG(4, "-> Diag::Diag")

  tangentialVelocity.resize(m.nEdges * m.nVertLevels, -1.0e32);
  kineticEnergyCell.resize(m.nCells * m.nVertLevels, -1.0e32);
  divergence.resize(m.nCells * m.nVertLevels, -1.0e32);
  relativeVorticity.resize(m.nVertices * m.nVertLevels, -1.0e32);
}

void Diag::compute(Config &config, Mesh &m, State &s) {
  LOG(4, "-> Diag::compute")

  //*******************************************************
  // Loop over edges: tangential velocity
  //*******************************************************
  {
    size_t i, iEdge, k, eoe;
    for (iEdge = 0; iEdge < m.nEdges; iEdge++) {
      for (k = 0; k < m.K; k++) {
        tangentialVelocity[iEdge * m.K + k] = 0.0;
      }
      for (i = 0; i < m.nEdgesOnEdge[iEdge]; i++) {
        eoe = m.edgesOnEdge[iEdge * m.ME2 + i];
        for (k = 0; k < m.K; k++) {
          tangentialVelocity[iEdge * m.K + k] +=
              m.weightsOnEdge[iEdge * m.ME2 + i] *
              s.normalVelocity[eoe * m.K + k];
        }
      }
    }
  }

  //*******************************************************
  // Loop over vertices: relativeVorticity
  //*******************************************************
  {
    size_t iVertex, iEdge, i, k;
    double invAreaTri1;
    for (iVertex = 0; iVertex < m.nVertices; iVertex++) {
      invAreaTri1 = 1.0 / m.areaTriangle[iVertex];
      for (k = 0; k < m.K; k++) {
        relativeVorticity[iVertex * m.K + k] = 0.0;
      }
      for (i = 0; i < m.VD; i++) {
        iEdge = m.edgesOnVertex[iVertex * m.VD + i];
        for (k = 0; k < m.K; k++) {
          relativeVorticity[iVertex * m.K + k] +=
              m.edgeSignOnVertex[iVertex * m.VD + i] * m.dcEdge[iEdge] *
              s.normalVelocity[iEdge * m.K + k] * invAreaTri1;
          // maybe later:
          // kineticEnergyVertex[k,iVertex] +=
          // m.dcEdge[iEdge]*pow(m.dvEdge[iEdge],0.25)/areaTriangle[iVertex]*std::pow(s.normalVelocity[k,iEdge],2);
        }
      }
    }
  }

  //*******************************************************
  // Loop over cells: divergence, kineticEnergyCell
  //*******************************************************
  {
    size_t iEdge, i, k;
    double invAreaCell, r_tmp;
    signed char edgeSignOnCell_temp;
    for (size_t iCell = 0; iCell < m.nCells; iCell++) {
      invAreaCell = 1.0 / m.areaCell[iCell];
      for (k = 0; k < m.K; k++) {
        divergence[iCell * m.K + k] = 0.0;
        kineticEnergyCell[iCell * m.K + k] = 0.0;
      }
      for (i = 0; i < m.nEdgesOnCell[iCell]; i++) {
        iEdge = m.edgesOnCell[iCell * m.ME + i];
        edgeSignOnCell_temp = m.edgeSignOnCell[iCell * m.ME + i];
        for (k = 0; k < m.K; k++) {
          r_tmp =
              m.dvEdge[iEdge] * s.normalVelocity[iEdge * m.K + k] * invAreaCell;

          divergence[iCell * m.K + k] -= edgeSignOnCell_temp * r_tmp;
          kineticEnergyCell[iCell * m.K + k] +=
              0.25 * r_tmp * m.dcEdge[iEdge] *
              s.normalVelocity[iEdge * m.K + k];
          // maybe later:
          // div_hu[k] -= layerThicknessEdgeFlux[k,iEdge]*
          // edgeSignOnCell_temp*r_tmp;
        }
      }
    }
  }
}

//*******************************************************
// Loop over edges
//*******************************************************
//      do iVertex = 1, nVerticesAll
//         invAreaTri1 = 1.0 / areaTriangle(iVertex)
//         do i = 1, vertexDegree
//            iEdge = edgesOnVertex(i, iVertex)
//            do k = minLevelVertexTop(iVertex), maxLevelVertexBot(iVertex)
//              r_tmp = dcEdge(iEdge) * normalVelocity(k, iEdge)
//              relativeVorticity(k, iVertex) = relativeVorticity(k, iVertex) +
//              edgeSignOnVertex(i, iVertex) * r_tmp * invAreaTri1
//            end do
//         end do
//      end do
//    //kineticEnergyVertex, for KE mix, if I want it later? (optional, consider
//    later)
//      do iVertex = 1, nVerticesAll
//         kineticEnergyVertex(:, iVertex) = 0.0
//         do i = 1, vertexDegree
//            iEdge = edgesOnVertex(i, iVertex)
//            rTmp  = dcEdge(iEdge)*dvEdge(iEdge)*0.25/ &
//                    areaTriangle(iVertex)
//            do k = 1, nVertLevels
//               kineticEnergyVertex(k,iVertex) = &
//               kineticEnergyVertex(k,iVertex) + rTmp* &
//                            normalVelocity(k,iEdge)**2
//            end do
//         end do
//      end do
// relativeVorticity
// kineticEnergyVertex, for KE mix, if I want it later? (optional, consider
// later)
//      do iVertex = 1, nVerticesAll
//         kineticEnergyVertex(:, iVertex) = 0.0
//         do i = 1, vertexDegree
//            iEdge = edgesOnVertex(i, iVertex)
//            rTmp  = dcEdge(iEdge)*dvEdge(iEdge)*0.25/ &
//                    areaTriangle(iVertex)
//            do k = 1, nVertLevels
//               kineticEnergyVertex(k,iVertex) = &
//               kineticEnergyVertex(k,iVertex) + rTmp* &
//                            normalVelocity(k,iEdge)**2
//            end do
//         end do
//      end do
//*******************************************************
// Loop over edges
//*******************************************************
// tangentialVelocity
//   do iEdge = 1, nEdges
//      tangentialVelocity(:,iEdge) = 0.0
//      do i = 1, nEdgesOnEdge(iEdge)
//         eoe = edgesOnEdge(i,iEdge)
//         weightsOnEdge_temp = weightsOnEdge(i, iEdge)
//         do k = kmin,kmax
//            tangentialVelocity(k,iEdge) = &
//            tangentialVelocity(k,iEdge) + weightsOnEdge_temp* &
//                                          normalVelocity(k, eoe)
//         end do
//      end do
//   end do
//*******************************************************
//  Loop over cells
//*******************************************************
// to do: kineticEnergyCell, from KE vertex mix (optional, consider later)

//      do iCell = 1, nCells
//         divergence(:,iCell) = 0.0
//         kineticEnergyCell(:,iCell) = 0.0
//         div_hu(:) = 0.0
//         invAreaCell1 = invAreaCell(iCell)
//         do i = 1, nEdgesOnCell(iCell)
//            iEdge = edgesOnCell(i, iCell)
//            edgeSignOnCell_temp = edgeSignOnCell(i, iCell)
//            dcEdge_temp = dcEdge(iEdge)
//            dvEdge_temp = dvEdge(iEdge)
//            do k = kmin,kmax
//               r_tmp = dvEdge_temp*normalVelocity(k,iEdge)*invAreaCell1
//
//               divergence(k,iCell) = divergence(k,iCell) -
//               edgeSignOnCell_temp*r_tmp div_hu(k) = div_hu(k) -
//               layerThicknessEdgeFlux(k,iEdge)* edgeSignOnCell_temp*r_tmp
//               kineticEnergyCell(k,iCell) = kineticEnergyCell(k,iCell) +
//               0.25*r_tmp*dcEdge_temp* normalVelocity(k,iEdge)
//            end do
//         end do
//    //kineticEnergyCell, from KE vertex mix (optional, consider later)
//      do iCell = 1, nCells
//         kineticEnergyVertexOnCells(:, iCell) = 0.0
//         invAreaCell1 = invAreaCell(iCell)
//         do i = 1, nEdgesOnCell(iCell)
//            j = kiteIndexOnCell(i, iCell)
//            iVertex = verticesOnCell(i, iCell)
//            do k = 1, nVertLevels
//               kineticEnergyVertexOnCells(k,iCell) = &
//               kineticEnergyVertexOnCells(k,iCell) + &
//                        kiteAreasOnVertex(j,iVertex)* &
//                        kineticEnergyVertex(k,iVertex)*invAreaCell1
//            end do
//         end do
//      end do
//      do iCell = 1, nCells
//      do k = 1, nVertLevels
//         kineticEnergyCell(k,iCell) = 5.0/8.0* &
//                                      kineticEnergyCell(k,iCell) &
//                                    + 3.0/8.0* &
//                                      kineticEnergyVertexOnCells(k,iCell)
//      end do
//      end do

//    tangentialVelocity.resize(m.nEdges * m.nVertLevels, -1.0e32); //
//    type="real" dimensions="nVertLevels nEdges Time" units="m s^-1"
//    description="horizontal velocity, tangential to an edge"
//    kineticEnergyCell.resize(m.nCells * m.nVertLevels, -1.0e32); //
//    type="real" dimensions="nVertLevels nCells Time" units="m^2 s^-2"
//    description="kinetic energy of horizontal velocity on cells"
//    divergence.resize(m.nCells * m.nVertLevels, -1.0e32); // type="real"
//    dimensions="nVertLevels nCells Time" units="s^-1" description="divergence
//    of horizontal velocity" relativeVorticity.resize(m.nVertices *
//    m.nVertLevels, -1.0e32); // type="real" dimensions="nVertLevels nVertices
//    Time" units="s^-1" description="curl of horizontal velocity, defined at
//    vertices"

// unused for now:
//    relativeVorticityCell.resize(m.nCells * m.nVertLevels, -1.0e32); //
//    type="real" dimensions="nVertLevels nCells Time" units="s^-1"
//    description="curl of horizontal velocity, averaged from vertices to cell
//    centers" normalizedRelativeVorticityEdge.resize(m.nEdges * m.nVertLevels,
//    -1.0e32); // type="real" dimensions="nVertLevels nEdges Time" units="s^-1"
//    description="curl of horizontal velocity divided by layer thickness,
//    averaged from vertices to edges"
//    normalizedPlanetaryVorticityEdge.resize(m.nEdges * m.nVertLevels,
//    -1.0e32); // type="real" dimensions="nVertLevels nEdges Time" units="s^-1"
//    description="earth's rotational rate (Coriolis parameter, f) divided by
//    layer thickness, averaged from vertices to edges"
//    normalizedRelativeVorticityCell.resize(m.nCells * m.nVertLevels, -1.0e32);
//    // type="real" dimensions="nVertLevels nCells Time" units="s^-1"
//    description="curl of horizontal velocity divided by layer thickness,
//    averaged from vertices to cell centers"
//    layerThicknessEdgeMean.resize(m.nEdges * m.nVertLevels, -1.0e32); //
//    type="real" dimensions="nVertLevels nEdges Time" units="m"
//    description="layer thickness averaged from cell center to edges"
//    layerThicknessEdgeFlux.resize(m.nEdges * m.nVertLevels, -1.0e32); //
//    type="real" dimensions="nVertLevels nEdges Time" units="m"
//    description="layer thickness used for fluxes through edges. May be
//    centered, upwinded, or a combination of the two."
//    layerThicknessVertex.resize(m.nVertices * m.nVertLevels, -1.0e32); //
//    type="real" dimensions="nVertLevels nVertices Time" units="m"
//    description="layer thickness averaged from cell center to vertices"
//    circulation.resize(m.nVertices * m.nVertLevels, -1.0e32); // type="real"
//    dimensions="nVertLevels nVertices Time" units="m^2 s^-1"
//    description="area-integrated vorticity"
