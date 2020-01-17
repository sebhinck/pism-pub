/* Copyright (C) 2016, 2017, 2018, 2019 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "remove_narrow_tongues.hh"

#include "pism/util/IceGrid.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/util/Mask.hh"
#include "pism/util/ConfigInterface.hh"

namespace pism {

/** Remove tips of one-cell-wide ice tongues ("noses")..
 *
 * The center icy cell in ice tongues like this one (and equivalent)
 *
 * @code
   O O ?
   X X O
   O O ?
   @endcode
 *
 * where "O" is ice-free and "?" is any mask value, are removed.
 * Ice tongues like this one
 *
 * @code
   # O ?
   X X O
   # O ?
   @endcode
 * where one or two of the "#" cells are ice-filled, are not removed.
 *
 * See the code for the precise rule, which uses `ice_free_ocean()` for the "O"
 * cells if the center cell has grounded ice, and uses `ice_free()` if the
 * center cell has floating ice.
 *
 * @note We use `pism_mask` (and not ice_thickness) to make decisions.
 * This means that we can update `ice_thickness` in place without
 * introducing a dependence on the grid traversal order.
 *
 * @param[in,out] mask cell type mask
 * @param[in,out] ice_thickness modeled ice thickness
 *
 * @return 0 on success
 */
void remove_narrow_tongues(const Geometry &geometry,
                           IceModelVec2S &ice_thickness) {

  auto &mask       = geometry.cell_type;
  auto &bed        = geometry.bed_elevation;
  auto &sea_level  = geometry.sea_level_elevation;
  auto &lake_level = geometry.lake_level_elevation;

  IceGrid::ConstPtr grid = mask.grid();
  Config::ConstPtr config = grid->ctx()->config();
  GeometryCalculator gc(*config);

  IceModelVec::AccessList list{&mask, &bed, &sea_level, &lake_level, &ice_thickness};

  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    const double water_level = gc.water_level(sea_level(i, j), bed(i, j), lake_level(i, j));
    if (mask.ice_free(i,j) or
        (mask.grounded_ice(i,j) and bed(i,j) >= water_level)) {
      continue;
    }

    bool
      ice_free_N  = false, ice_free_E  = false,
      ice_free_S  = false, ice_free_W  = false,
      ice_free_NE = false, ice_free_NW = false,
      ice_free_SE = false, ice_free_SW = false;

    if (mask.grounded_ice(i,j)) {
      // if (i,j) is grounded ice then we will remove it if it has
      // exclusively ice-free ocean neighbors
      ice_free_N  = mask.ice_free_ocean(i, j + 1);
      ice_free_E  = mask.ice_free_ocean(i + 1, j);
      ice_free_S  = mask.ice_free_ocean(i, j - 1);
      ice_free_W  = mask.ice_free_ocean(i - 1, j);
      ice_free_NE = mask.ice_free_ocean(i + 1, j + 1);
      ice_free_NW = mask.ice_free_ocean(i - 1, j + 1);
      ice_free_SE = mask.ice_free_ocean(i + 1, j - 1);
      ice_free_SW = mask.ice_free_ocean(i - 1, j - 1);
    } else if (mask.floating_ice(i,j)) {
      // if (i,j) is floating then we will remove it if its neighbors are
      // ice-free, whether ice-free ocean or ice-free ground
      ice_free_N  = mask.ice_free(i, j + 1);
      ice_free_E  = mask.ice_free(i + 1, j);
      ice_free_S  = mask.ice_free(i, j - 1);
      ice_free_W  = mask.ice_free(i - 1, j);
      ice_free_NE = mask.ice_free(i + 1, j + 1);
      ice_free_NW = mask.ice_free(i - 1, j + 1);
      ice_free_SE = mask.ice_free(i + 1, j - 1);
      ice_free_SW = mask.ice_free(i - 1, j - 1);
    }

    if ((not ice_free_W and
         ice_free_NW    and
         ice_free_SW    and
         ice_free_N     and
         ice_free_S     and
         ice_free_E)    or
        (not ice_free_N and
         ice_free_NW    and
         ice_free_NE    and
         ice_free_W     and
         ice_free_E     and
         ice_free_S)    or
        (not ice_free_E and
         ice_free_NE    and
         ice_free_SE    and
         ice_free_W     and
         ice_free_S     and
         ice_free_N)    or
        (not ice_free_S and
         ice_free_SW    and
         ice_free_SE    and
         ice_free_W     and
         ice_free_E     and
         ice_free_N)) {
      ice_thickness(i, j) = 0.0;
    }
  }
}

} // end of namespace pism
