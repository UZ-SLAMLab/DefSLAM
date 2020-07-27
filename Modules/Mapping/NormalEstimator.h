/**
* This file is part of DefSLAM.
*
* Copyright (C) 2017-2020 Jose Lamarca Peiro <jlamarca at unizar dot es>, J.M.M. Montiel (University
*of Zaragoza) && Shaifali Parashar, Adrien Bartoli (Université Clermont Auvergne)
* For more information see <https://github.com/unizar/DefSLAM>
*
* DefSLAM is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* DefSLAM is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with DefSLAM. If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef NORMALESTIMATOR_H
#define NORMALESTIMATOR_H

#include "PolySolver.h"
#include "Surface.h"
#include "SurfacePoint.h"
#include "WarpDatabase.h"
#include "ceres/ceres.h"

namespace defSLAM
{

  class PolySolver;
  class SchwarpDatabase;
  class Surface;
  class SurfacePoint;
  using ORB_SLAM2::KeyFrame;

  class NormalEstimator
  {
  public:
    // Construtor by default deleted.
    NormalEstimator() = delete;
    // The program is initialized with the database of warps that contain the information needed
    // to estimate the normals.
    NormalEstimator(WarpDatabase *WarpDatabase);

    ~NormalEstimator();
    /***********
     * Core function to estimate K1 and K2 that are the two first components of the normal.
     * n = [k1, k2 , 1- k1*u-k2*v]
     * *********/
    void ObtainK1K2();

  public:
    WarpDatabase *WarpDB;
  };

} // namespace defSLAM

#endif // DEFORMATION MAPPOINT_H
