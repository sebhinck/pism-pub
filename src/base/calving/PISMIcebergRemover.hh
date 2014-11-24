/* Copyright (C) 2013, 2014 PISM Authors
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

#ifndef _PISMICEBERGREMOVER_H_
#define _PISMICEBERGREMOVER_H_

#include "PISMComponent.hh"
#include "iceModelVec.hh"

namespace pism {

/*! \brief PISM iceberg remover */
/*!
 * Identifies and removes free-floating icebergs, which cause
 * well-posedness problems for stress solvers.
 *
 * Icebergs are, in this context, floating regions that are _not_
 * attached, through a chain of positive thickness ice-filled cells,
 * to at least one grounded cell.
 *
 * They cause the SSA operator to have a nontrivial null space.
 *
 * They are observed to cause unrealistically large velocities that
 * may affect ice velocities elsewhere.
 *
 * This class uses a serial connected component labeling algorithm to
 * remove "icebergs".
 */
class IcebergRemover : public Component
{
public:
  IcebergRemover(IceGrid &g, const Config &conf);
  virtual ~IcebergRemover();

  virtual void init(Vars &vars);
  void update(IceModelVec2Int &pism_mask, IceModelVec2S &ice_thickness);
protected:

  virtual void add_vars_to_output(const std::string &keyword, std::set<std::string> &result);
  virtual void define_variables(const std::set<std::string> &vars, const PIO &nc,
                                          IO_Type nctype);
  virtual void write_variables(const std::set<std::string> &vars, const PIO& nc);
  
  IceModelVec2S m_iceberg_mask;
  Vec m_mask_p0;

  IceModelVec2Int *m_bcflag;
};

} // end of namespace pism

#endif /* _PISMICEBERGREMOVER_H_ */
