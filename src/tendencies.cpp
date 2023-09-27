#include "tendencies.h"
#include "Config.h"
#include "Diag.h"
#include "Mesh.h"
#include "Meta.h"
#include "State.h"
#include "Tend.h"
#include <iostream>

void uTend_KE_gradient(Config &config, Meta &meta, Mesh &m, State &s, Diag &d,
                       Tend &tend) {
  LOG(4, "-> uTend_KE_gradient")
  if (!config.uTend_KE_gradient_enable)
    return;

  size_t k, cell1, cell2;
  for (size_t iEdge = 0; iEdge < m.nEdges; iEdge++) {
    cell1 = m.cellsOnEdge[iEdge * 2];
    cell2 = m.cellsOnEdge[iEdge * 2 + 1];
    for (k = 0; k < m.K; k++) {
      tend.normalVelocity[iEdge * m.K + k] +=
          d.kineticEnergyCell[cell2 * m.K + k] -
          d.kineticEnergyCell[cell1 * m.K + k];
      // note, MPAS-O uses edgeMask(k,iEdge) as coefficient.
    }
  }
}

void uTend_ssh_gradient(Config &config, Meta &meta, Mesh &m, State &s, Diag &d,
                        Tend &tend) {
  LOG(4, "-> uTend_ssh_gradient")
  if (!config.uTend_ssh_gradient_enable)
    return;

  size_t k, cell1, cell2;
  for (size_t iEdge = 0; iEdge < m.nEdges; iEdge++) {
    cell1 = m.cellsOnEdge[iEdge * 2];
    cell2 = m.cellsOnEdge[iEdge * 2 + 1];
    for (k = 0; k < m.K; k++) {
      tend.normalVelocity[iEdge * m.K + k] +=
          gravity *
          ((-m.bottomDepth[cell2] + s.layerThickness[cell2 * m.K + k]) -
           (-m.bottomDepth[cell1] + s.layerThickness[cell1 * m.K + k]));
      // note, MPAS-O uses edgeMask(k,iEdge) as coefficient.
    }
  }
}

void uTend_advection(Config &config, Meta &meta, Mesh &m, State &s, Diag &d,
                     Tend &tend) {
  LOG(4, "-> uTend_advection")
  if (!config.uTend_advection_enable)
    return;

  // Option 1: \left[ \frac{\omega_v +f_v}{[h_i]_v}\right]_e\left([h_i]_e
  // u_e^{\perp}\right) mrp add option 1 later.
  //         do iEdge = 1, nEdgesAll
  //         do k=1,nVertLevels
  //            tmpVorticity(k,iEdge) = normRelVortEdge(k,iEdge) + &
  //                                 normPlanetVortEdge(k,iEdge)
  //         end do
  //         end do
  //         do iEdge = 1, nEdgesOwned
  //         cell1 = cellsOnEdge(1,iEdge)
  //         cell2 = cellsOnEdge(2,iEdge)
  //         invLength = 1.0_RKIND / dcEdge(iEdge)
  //
  //         do k = kmin, kmax
  //            tend(k,iEdge) = tend(k,iEdge) + &
  //                        edgeMask(k,iEdge)* (qArr(k,iEdge) - &
  //                           (kineticEnergyCell(k,cell2) &
  //                          - kineticEnergyCell(k,cell1))*invLength)
  //         end do
  //      end do
  //      do iEdge = 1, nEdgesOwned
  //         cell1 = cellsOnEdge(1,iEdge)
  //         cell2 = cellsOnEdge(2,iEdge)
  //         invLength = 1.0_RKIND / dcEdge(iEdge)
  //         do k = kmin, kmax
  //            tend(k,iEdge) = tend(k,iEdge) + &
  //                        edgeMask(k,iEdge)* (qArr(k,iEdge) - &
  //                           (kineticEnergyCell(k,cell2) &
  //                          - kineticEnergyCell(k,cell1))*invLength)
  //         end do
  //      end do

  // Option 2: ([\omega_v]_e +f_e) u_e^\perp
  {
    size_t k, vertex1, vertex2;
    double relativeVorticityEdge;
    for (size_t iEdge = 0; iEdge < m.nEdges; iEdge++) {
      vertex1 = m.verticesOnEdge[iEdge * 2];
      vertex2 = m.verticesOnEdge[iEdge * 2 + 1];
      for (k = 0; k < m.K; k++) {
        relativeVorticityEdge = 0.5 * (d.relativeVorticity[vertex1 * m.K + k] +
                                       d.relativeVorticity[vertex2 * m.K + k]);
        tend.normalVelocity[iEdge * m.K + k] +=
            (relativeVorticityEdge + config.coriolis) *
            d.tangentialVelocity[iEdge * m.K + k];
      }
    }
  }
}

void uTend_del2(Config &config, Meta &meta, Mesh &m, State &s, Diag &d,
                Tend &tend) {
  LOG(4, "-> uTend_del2")
  if (!config.uTend_del2_enable)
    return;

  //   do iEdge = 1, nEdgesOwned
  //      cell1 = cellsOnEdge(1,iEdge)
  //      cell2 = cellsOnEdge(2,iEdge)
  //      vertex1 = verticesOnEdge(1,iEdge)
  //      vertex2 = verticesOnEdge(2,iEdge)
  //
  //      dcEdgeInv = 1.0_RKIND / dcEdge(iEdge)
  //      dvEdgeInv = 1.0_RKIND / dvEdge(iEdge)
  //
  //      visc2 =  viscDel2*meshScalingDel2(iEdge)
  //
  //      do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
  //
  //         ! Here -( relativeVorticity(k,vertex2) -
  //         !         relativeVorticity(k,vertex1) ) / dvEdge(iEdge)
  //         ! is - \nabla relativeVorticity pointing from vertex 2
  //         ! to vertex 1, or equivalently
  //         ! + k \times \nabla relativeVorticity pointing from cell1
  //         ! to cell2.
  //
  //         uDiff = (div(k,cell2) - div(k,cell1))*dcEdgeInv &
  //                -(relVort(k,vertex2)-relVort(k,vertex1))*dvEdgeInv
  //
  //         tend(k,iEdge) = tend(k,iEdge) + &
  //                         edgeMask(k,iEdge)*visc2*uDiff
  //
  //      end do
  //   end do

  size_t k, cell1, cell2;
  size_t vertex1, vertex2;
  double dcEdgeInv, dvEdgeInv;
  for (size_t iEdge = 0; iEdge < m.nEdges; iEdge++) {
    cell1 = m.cellsOnEdge[iEdge * 2];
    cell2 = m.cellsOnEdge[iEdge * 2 + 1];
    vertex1 = m.verticesOnEdge[iEdge * 2];
    vertex2 = m.verticesOnEdge[iEdge * 2 + 1];
    dcEdgeInv = 1.0 / m.dcEdge[iEdge];
    dvEdgeInv = 1.0 / m.dvEdge[iEdge];
    //      visc2 =  viscDel2*meshScalingDel2(iEdge)
    for (k = 0; k < m.K; k++) {
      // Here -( relativeVorticity(k,vertex2) -
      //         relativeVorticity(k,vertex1) ) / dvEdge(iEdge)
      // is - \nabla relativeVorticity pointing from vertex 2
      // to vertex 1, or equivalently
      // + k \times \nabla relativeVorticity pointing from cell1
      // to cell2.

      // MPAS-O uses edgeMask(k,iEdge)
      tend.normalVelocity[iEdge * m.K + k] +=
          config.uTend_del2_coef *
              (d.divergence[cell2 * m.K + k] - d.divergence[cell1 * m.K + k]) *
              dcEdgeInv -
          (d.relativeVorticity[vertex2 * m.K + k] -
           d.relativeVorticity[vertex1 * m.K + k]) *
              dvEdgeInv;
    }
  }
}

void uTend_del4(Config &config, Meta &meta, Mesh &m, State &s, Diag &d,
                Tend &tend) {
  LOG(4, "-> uTend_del4")
  if (!config.uTend_del4_enable)
    return;

  size_t k, cell1, cell2;
  for (size_t iEdge = 0; iEdge < m.nEdges; iEdge++) {
    for (k = 0; k < m.K; k++) {
      // tend.normalVelocity[iEdge*m.K+k] += config.uTend_del4 *
      // s.normalVelocity[iEdge*m.K+k];
    }
  }
}

void uTend_bottom_drag(Config &config, Meta &meta, Mesh &m, State &s, Diag &d,
                       Tend &tend) {
  LOG(4, "-> uTend_bottom_drag")
  if (!config.uTend_bottom_drag_enable)
    return;

  size_t k, cell1, cell2;
  for (size_t iEdge = 0; iEdge < m.nEdges; iEdge++) {
    for (k = 0; k < m.K; k++) {
      // tend.normalVelocity[iEdge*m.K+k] += config.uTend_bottom_drag *
      // s.normalVelocity[iEdge*m.K+k];
    }
  }
}

void uTend_wind_forcing(Config &config, Meta &meta, Mesh &m, State &s, Diag &d,
                        Tend &tend) {
  LOG(4, "-> uTend_wind_forcing")
  if (!config.uTend_wind_forcing_enable)
    return;

  size_t k, cell1, cell2;
  for (size_t iEdge = 0; iEdge < m.nEdges; iEdge++) {
    for (k = 0; k < m.K; k++) {
      // tend.normalVelocity[iEdge*m.K+k] += config.uTend_wind_forcing *
      // s.normalVelocity[iEdge*m.K+k];
    }
  }
}

void uTend_Rayleigh(Config &config, Meta &meta, Mesh &m, State &s, Diag &d,
                    Tend &tend) {
  LOG(4, "-> uTend_Rayleigh")
  if (!config.uTend_Rayleigh_enable)
    return;

  size_t k, cell1, cell2;
  for (size_t iEdge = 0; iEdge < m.nEdges; iEdge++) {
    for (k = 0; k < m.K; k++) {
      tend.normalVelocity[iEdge * m.K + k] -=
          config.uTend_Rayleigh_drag * s.normalVelocity[iEdge * m.K + k];
    }
  }
}

void hTend_advection(Config &config, Meta &meta, Mesh &m, State &s, Diag &d,
                     Tend &tend) {
  LOG(4, "-> hTend_advection")
  if (!config.hTend_advection_enable)
    return;

  size_t iEdge, i, k;
  double invAreaCell, r_tmp;
  signed char edgeSignOnCell_temp;
  for (size_t iCell = 0; iCell < m.nCells; iCell++) {
    invAreaCell = 1.0 / m.areaCell[iCell];
    for (i = 0; i < m.nEdgesOnCell[iCell]; i++) {
      iEdge = m.edgesOnCell[iCell * m.ME + i];
      edgeSignOnCell_temp = m.edgeSignOnCell[iCell * m.ME + i];
      for (k = 0; k < m.K; k++) {
        r_tmp =
            m.dvEdge[iEdge] * s.normalVelocity[iEdge * m.K + k] * invAreaCell;

        tend.layerThickness[iCell * m.K + k] -= edgeSignOnCell_temp * r_tmp;
      }
    }
  }
}

void hTend_decay(Config &config, Meta &meta, Mesh &m, State &s, Diag &d,
                 Tend &tend) {
  LOG(4, "-> hTend_decay")
  if (!config.hTend_decay_enable)
    return;

  size_t k, cell1, cell2;
  for (size_t iCell = 0; iCell < m.nCells; iCell++) {
    for (k = 0; k < m.K; k++) {
      tend.layerThickness[iCell * m.K + k] -=
          config.hTend_decay_coef * s.layerThickness[iCell * m.K + k];
    }
  }
}

void hTend_del2(Config &config, Meta &meta, Mesh &m, State &s, Diag &d,
                Tend &tend) {
  LOG(4, "-> hTend_del2")
  if (!config.hTend_del2_enable)
    return;

  size_t iCell, i, k, cell1, cell2, iEdge;
  double invAreaCell, r_tmp, tracer_turb_flux, flux;
  for (iCell = 0; iCell < m.nCells; iCell++) {
    invAreaCell = 1.0 / m.areaCell[iCell];
    for (i = 0; i < m.nEdgesOnCell[iCell]; i++) {
      iEdge = m.edgesOnCell.at(iCell * m.ME + i);
      cell1 = m.cellsOnEdge.at(iEdge * 2);
      cell2 = m.cellsOnEdge.at(iEdge * 2 + 1);
      r_tmp = config.hTend_del2_coef * m.dvEdge[iEdge] / m.dcEdge[iEdge];
      for (k = 0; k < m.K; k++) {
        //        ! \kappa_2 \nabla \phi on edge
        tracer_turb_flux = s.layerThickness.at(cell2 * m.K + k) -
                           s.layerThickness.at(cell1 * m.K + k);
        //        ! div(h \kappa_2 \nabla \phi) at cell center, but no h
        //        coefficient here
        flux = tracer_turb_flux * r_tmp;
        tend.layerThickness[iCell * m.K + k] -=
            m.edgeSignOnCell[iCell * m.ME + i] * flux * invAreaCell;
      }
    }
  }
  //    do iCell = 1, nCells
  //    invAreaCell = 1.0_RKIND / areaCell(iCell)
  //    do i = 1, nEdgesOnCell(iCell)
  //      iEdge = edgesOnCell(i, iCell)
  //      cell1 = cellsOnEdge(1,iEdge)
  //      cell2 = cellsOnEdge(2,iEdge)
  //
  //      r_tmp = meshScalingDel2(iEdge) * eddyDiff2 * dvEdge(iEdge) /
  //      dcEdge(iEdge)
  //
  //      do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
  //      do iTracer = 1, num_tracers
  //        ! \kappa_2 \nabla \phi on edge
  //        tracer_turb_flux = tracers(iTracer, k, cell2) - tracers(iTracer, k,
  //        cell1)
  //
  //        ! div(h \kappa_2 \nabla \phi) at cell center
  //        flux = layerThickEdgeMean(k, iEdge) * tracer_turb_flux * r_tmp
  //
  //        tend(iTracer, k, iCell) = tend(iTracer, k, iCell) -
  //        edgeSignOnCell(i, iCell) * flux * invAreaCell
  //      end do
  //      end do
  //
  //    end do
  //    end do
}

void compute_velocity_tendencies(Config &config, Meta &meta, Mesh &m, State &s,
                                 Diag &d, Tend &tend) {
  LOG(4, "-> compute_velocity_tendencies")

  fill(tend.normalVelocity.begin(), tend.normalVelocity.end(), 0.0);
  uTend_advection(config, meta, m, s, d, tend);
  uTend_ssh_gradient(config, meta, m, s, d, tend);
  uTend_del2(config, meta, m, s, d, tend);
  uTend_del4(config, meta, m, s, d, tend);
  uTend_bottom_drag(config, meta, m, s, d, tend);
  uTend_wind_forcing(config, meta, m, s, d, tend);
  uTend_Rayleigh(config, meta, m, s, d, tend);
}

void compute_thickness_tendencies(Config &config, Meta &meta, Mesh &m, State &s,
                                  Diag &d, Tend &tend) {
  LOG(4, "-> compute_thickness_tendencies")

  fill(tend.layerThickness.begin(), tend.layerThickness.end(), 0.0);
  hTend_advection(config, meta, m, s, d, tend);
  hTend_del2(config, meta, m, s, d, tend);
  hTend_decay(config, meta, m, s, d, tend);
}
