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

#include "SurfacePoint.h"

namespace defSLAM
{
  // Constructor.
  SurfacePoint::SurfacePoint()
      : NormalOn(false) {}

  SurfacePoint::~SurfacePoint() {}

  void SurfacePoint::setNormal(cv::Vec3f &N)
  {
    Normal = N;
    NormalOn = true;
  }

  bool SurfacePoint::getNormal(cv::Vec3f &N)
  {
    if (NormalOn)
    {
      N = Normal;
      return true;
    }
    else
    {
      return false;
    }
  }
  void SurfacePoint::setDepth(cv::Vec3f &x3d) { x3D = x3d; }

  void SurfacePoint::getDepth(cv::Vec3f &x3d) { x3d = x3D; }

  bool SurfacePoint::thereisNormal() { return NormalOn; }
} // namespace defSLAM
