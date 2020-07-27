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

#ifndef SURFACEPOINT_H
#define SURFACEPOINT_H

#include <opencv2/core/core.hpp>

namespace defSLAM
{

  class SurfacePoint
  {
  public:
    // Constructor.
    SurfacePoint();

    // Destructor.
    ~SurfacePoint();

    // Set normal
    void setNormal(cv::Vec3f &N);

    // Get normal
    bool getNormal(cv::Vec3f &N);

    // Check if the normal was set
    bool thereisNormal();

    // Set depth
    void setDepth(cv::Vec3f &N);

    // Get depth
    void getDepth(cv::Vec3f &x3D);

  private:
    bool NormalOn;    // Normal is set.
    cv::Vec3f Normal; // Normal
    cv::Vec3f x3D;    // Position 3D
  };

} // namespace defSLAM

#endif // SURFACEPOINT_H
