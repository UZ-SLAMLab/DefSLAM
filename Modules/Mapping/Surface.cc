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

#include "Surface.h"

namespace BBS
{
 
  void EvalEigen(bbs_t* bbs, const std::vector<double>& ctrlpts, std::vector<double>& u, std::vector<double>& v,
      int nval, Eigen::MatrixXd& val, int du, int dv);
}

namespace defSLAM
{
  // Constructor
  Surface::Surface(unsigned int numberOfPoints)
      : nodesDepth_(), surfacePoints_(numberOfPoints, SurfacePoint()),
        numberofNormals_(0), normalsLimit_(10)
  {
  }

  // Destructor
  Surface::~Surface()
  {
    surfacePoints_.clear();
   
  }

  /***********
  * Save the array with the depth estimated by 
  ************/
  void Surface::saveArray(const std::vector<double>& Array, BBS::bbs_t& bbss)
  {
      bbs = bbss;
      nodesDepth_.resize(bbs.nptsu * bbs.nptsv);
      copy(Array.begin(), Array.begin()+ bbs.nptsu * bbs.nptsv, nodesDepth_.begin()); // NOTE should till the END, but I don\t want to change logic

  }

  /***********
  * Check if enough normals have been estimated
  * to define the surface with Shape-from-normals.
  ************/
  bool Surface::enoughNormals()
  {
    std::cout << "Number Of normals " << numberofNormals_ << " " << this
              << std::endl;
    return (numberofNormals_ >= normalsLimit_);
  }

  // Return normals saved
  unsigned int Surface::getnumberofNormals() { return numberofNormals_; }

  /***********
    * Set a normal saved. In DefSLAM ind makes reference to
    * the keypoint index
    ************/
  void Surface::setNormalSurfacePoint(size_t ind, cv::Vec3f &N)
  {
    // Count normals just once
    if (!surfacePoints_[ind].thereisNormal())
    {
      numberofNormals_++;
    }

    surfacePoints_[ind].setNormal(N);
  }

  // Get the normal of the keypoint with index ind.
  bool Surface::getNormalSurfacePoint(size_t ind, cv::Vec3f &N)
  {
    return surfacePoints_[ind].getNormal(N);
  }

  // Set the 3D position of the keypoint with index ind.
  void Surface::set3DSurfacePoint(size_t ind, cv::Vec3f &x3D)
  {
    surfacePoints_[ind].setDepth(x3D);
  }

  // Get the 3D position of the keypoint with index ind.
  void Surface::get3DSurfacePoint(size_t ind, cv::Vec3f &x3D)
  {
    surfacePoints_[ind].getDepth(x3D);
  }

  // Apply scale to the entire surface. Called after surface
  // registration
  void Surface::applyScale(double s22)
  {
    for (unsigned int i(0); i < surfacePoints_.size(); i++)
    {
      cv::Vec3f x3D;
      surfacePoints_[i].getDepth(x3D);
      cv::Vec3f x3c = s22 * x3D;
      surfacePoints_[i].setDepth(x3c);
    }

    for (int i(0); i < bbs.nptsu * bbs.nptsv; i++)
    {
      nodesDepth_[i] = s22 * nodesDepth_[i];
    }
  }

  // Discretize the surface in xs*ys vertex. It is used to
  // create a mesh of xs columns and ys rows.
  void Surface::getVertex(std::vector<cv::Mat> &NodesSurface, unsigned int xs, unsigned int ys)
  {
    std::vector<double> arrayCU(xs * ys);
    std::vector<double> arrayCV(xs * ys);

    unsigned int us(0);
    double t(0.03);
    double umaxtemp = bbs.umax; //-0.20;
    double umintemp = bbs.umin; //+0.15;
    double vmaxtemp = bbs.vmax; //-0.20;
    double vmintemp = bbs.vmin; //+0.25;
    for (unsigned int x(0); x < xs; x++)
    {
      for (unsigned int j(0); j < ys; j++)
      {
        arrayCU[us] =
            double((umaxtemp - umintemp - 2 * t) * x) / (xs - 1) + (umintemp + t);
        arrayCV[us] =
            double((vmaxtemp - vmintemp - 2 * t) * j) / (ys - 1) + (vmintemp + t);
        us++;
      }
    }

    Eigen::MatrixXd Val2(xs * ys, 1);
    BBS::EvalEigen(&bbs, nodesDepth_, arrayCU,
                   arrayCV, xs * ys, Val2, 0, 0);
    NodesSurface.reserve(xs * ys);
    for (unsigned int x(0); x < xs * ys; x++)
    {
      cv::Mat x3D(4, 1, CV_32F);
      x3D.at<float>(0, 0) = arrayCU[x] * Val2(x, 0);
      x3D.at<float>(1, 0) = arrayCV[x] * Val2(x, 0);
      x3D.at<float>(2, 0) = Val2(x, 0);
      x3D.at<float>(3, 0) = 1;
      NodesSurface.push_back(x3D);
    }
  }
} // namespace defSLAM
