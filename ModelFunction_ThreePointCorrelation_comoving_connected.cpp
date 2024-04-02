/*******************************************************************
 *  Copyright (C) 2017 by Federico Marulli                         *
 *  federico.marulli3@unibo.it                                     *
 *                                                                 *
 *  This program is free software; you can redistribute it and/or  *
 *  modify it under the terms of the GNU General Public License as *
 *  published by the Free Software Foundation; either version 2 of *
 *  the License, or (at your option) any later version.            *
 *                                                                 *
 *  This program is distributed in the hope that it will be useful,*
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the  *
 *  GNU General Public License for more details.                   *
 *                                                                 *
 *  You should have received a copy of the GNU General Public      *
 *  License along with this program; if not, write to the Free     *
 *  Software Foundation, Inc.,                                     *
 *  59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.      *
 *******************************************************************/

/** 
 *  @file
 *  Modelling/ThreePointCorrelation/ModelFunction_ThreePointCorrelation_comoving_connected.cpp
 *
 *  @brief Functions to model the connected three-point correlation
 *  function in comoving coordinates
 *
 *  This file contains the implementation of the functions used to
 *  model the connected three-point correlation function in comoving
 *  coordinates
 *
 *  @authors Federico Marulli
 *
 *  @authors federico.marulli3@unibo.it
 */

#include "ModelFunction_ThreePointCorrelation.h"

#include "ModelFunction_ThreePointCorrelation_comoving_connected.h"

using namespace std;

using namespace cbl;


// =====================================================================================


std::vector<double> cbl::modelling::threept::zeta_RSD (const std::vector<double> theta, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  
  shared_ptr<STR_data_model_threept> pp = static_pointer_cast<STR_data_model_threept>(inputs);

  // input likelihood parameters

  // bias parameters
  
  double b1 = parameter[0];
  double b2 = parameter[1];
  double bt = parameter[2];

  // beta = f/b1
  
  double beta = parameter[3];

  cosmology::ThreePointCorrelation three_point_correlation(*pp->cosmology);
  
  std::vector<double> kk = logarithmic_bin_vector(200, 1.e-4, 100.);
  std::vector<double> Pk_matter = pp->cosmology->Pk_matter(kk, "CAMB", false, 1.);

  vector<double> model = three_point_correlation.zeta_RSD(pp->r1, pp->r2, theta, b1, b2, bt, beta, kk, Pk_matter);

  return model;

}


// =====================================================================================


std::vector<double> cbl::modelling::threept::zeta_RSD_all_scales (const std::vector<double> index, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  
  shared_ptr<STR_data_model_threept> pp = static_pointer_cast<STR_data_model_threept>(inputs);

  // input likelihood parameters

  // bias parameters
  
  double b1 = parameter[0];
  double b2 = parameter[1];
  double bt = parameter[2];

  // beta = f/b1
  
  double beta = parameter[3];

  cosmology::ThreePointCorrelation three_point_correlation(*pp->cosmology);
  
  std::vector<double> kk = logarithmic_bin_vector(200, 1.e-4, 100.);
  std::vector<double> Pk_matter = pp->cosmology->Pk_matter(kk, "CAMB", false, 1.);

  vector<double> model = three_point_correlation.zeta_RSD_all_scales(index, pp->r12, pp->r13, pp->r23, b1, b2, bt, beta, kk, Pk_matter);

  return model;

}