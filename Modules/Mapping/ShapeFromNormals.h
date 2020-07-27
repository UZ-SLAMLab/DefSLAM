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

#ifndef SHAPEFROMNORMALS_H
#define SHAPEFROMNORMALS_H

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "KeyFrame.h"
#include "Surface.h"
#include "WarpDatabase.h"
#include <ceres/ceres.h>

namespace defSLAM
{
  class Surface;
  using ORB_SLAM2::KeyFrame;
  using ORB_SLAM2::Map;

  // Structure for Eigen. It was created to avoid some issues
  // with the initialization
  struct SfNEigen
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Eigen::MatrixXd LinearSystem;
    Eigen::SparseMatrix<double> LinearSystemSparse;
    Eigen::MatrixXd M; // Matrix to impose the normals
    Eigen::MatrixXd B; // Bending regularizer.
    Eigen::MatrixXd X;
    Eigen::MatrixXd Solaux;
  };

  class ShapeFromNormals
  {
  public:
    // Constructor. It initializes the linear system.
    ShapeFromNormals(KeyFrame *refKf_,
                     double bendingWeight_);
    // destructor
    ~ShapeFromNormals() { delete sfnEigen_; }
    /*************
     * estimate(). Estimate a surface from the normals registered in
     * refKf_
     * ********/
    virtual bool estimate();

  private:
    /*************
     * obtainM(). Estimate M that is a matrix that makes the cross
     * product with the current control points and penalizes the 
     * surface estimated if does not accomplish with the normals.
     * ********/
    void obtainM(BBS::bbs_t &bbs, KeyFrame *Refdefkf, Eigen::MatrixXd &M);

  private:
    // Keyframe whose surface is estimated.
    KeyFrame *refKf_;

    double bendingWeight_; // Weight for the bending penalization

    std::vector<double> u_vector_;
    std::vector<double> v_vector_;

    int Mrows_;
    int Brows;

  private:
    void doInit() { sfnEigen_ = new SfNEigen; }
    SfNEigen *sfnEigen_;
  };

} // namespace defSLAM

#endif // LOCALMAPPING_H
