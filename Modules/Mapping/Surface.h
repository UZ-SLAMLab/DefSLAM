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

#ifndef SURFACE_H
#define SURFACE_H

#include "iostream"

#include "SurfacePoint.h"
#include "Thirdparty/BBS/bbs.h"
#include <Eigen/Dense>
#include <opencv2/core/core.hpp>

namespace defSLAM
{

  class SurfacePoint;

  class Surface
  {
  public:
    // No constructor by default
    Surface() = delete;

    // Number of keypoints in the image
    Surface(uint numberOfPoints);

    // Destructor
    ~Surface();

    /***********
     * Save the array with the depth estimated by 
     ************/
    void saveArray(double *Array, BBS::bbs_t &bbss);

    /***********
     * Check if enough normals have been estimated
     * to define the surface with Shape-from-normals.
     ************/
    bool enoughNormals();

    // Return normals saved
    uint getnumberofNormals();

    /***********
    * Set a normal saved. In DefSLAM ind makes reference to
    *the keypoint index
    ************/
    void setNormalSurfacePoint(size_t ind, cv::Vec3f &N);

    // Get the normal of the keypoint with index ind.
    bool getNormalSurfacePoint(size_t ind, cv::Vec3f &N);

    // Set the 3D position of the keypoint with index ind.
    void set3DSurfacePoint(size_t ind, cv::Vec3f &x3D);

    // Get the 3D position of the keypoint with index ind.
    void get3DSurfacePoint(size_t ind, cv::Vec3f &x3D);

    // Apply scale to the entire surface. Called after surface
    // registration
    void applyScale(double s22);

    // Discretize the surface in xs*ys vertex. It is used to
    // create a mesh of xs columns and ys rows.
    void getVertex(std::vector<cv::Mat> &NodesSurface, uint, uint);

  private:
    double *nodesDepth_; // Depth of the nodes register.

    std::vector<SurfacePoint> surfacePoints_; // Points correspoding to the keypoints of the keyframe

    uint numberofNormals_; // Number of normals registered.

    const uint normalsLimit_; // Limit of normals to estimate a surface.

    BBS::bbs_t bbs; // B-BSpline that represent the surface.
  };
} // namespace defSLAM

#endif
