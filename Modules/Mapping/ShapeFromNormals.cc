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

#include "ShapeFromNormals.h"
#include "DefKeyFrame.h"
#include "Thirdparty/BBS/bbs.h"
#include "Thirdparty/BBS/bbs_coloc.h"
#include <Thirdparty/BBS/bbs_MAC.h>

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <chrono>
#include <ctime>
#include <Eigen/QR>
#include <Eigen/SVD>

namespace defSLAM
{
  // Constructor. It initializes the linear system.
  ShapeFromNormals::ShapeFromNormals(KeyFrame *refKf_,
                                     double bendingWeight_)
      : refKf_(refKf_), bendingWeight_(bendingWeight_)
  {
    DefKeyFrame *refdefkf = static_cast<DefKeyFrame *>(refKf_);
    this->doInit();

    BBS::bbs_t bbsM;
    bbsM.umin = refdefkf->umin;
    bbsM.umax = refdefkf->umax;
    bbsM.nptsu = _NumberOfControlPointsU;
    bbsM.nptsv = _NumberOfControlPointsV;
    bbsM.vmax = refdefkf->vmax;
    bbsM.vmin = refdefkf->vmin;
    bbsM.valdim = 1;

    std::vector<Eigen::SparseMatrix<double>> Bs;
    Eigen::MatrixXd Mr;
    this->obtainM(bbsM, refKf_, Mr);

    Eigen::SparseMatrix<double> smBendinga(
        (_NumberOfControlPointsU * _NumberOfControlPointsV),
        (_NumberOfControlPointsU * _NumberOfControlPointsV));
    BBS::BendingEigen(&bbsM, bendingWeight_, smBendinga);

    this->Mrows_ = Mr.rows();

    this->Brows = _NumberOfControlPointsU *
                  _NumberOfControlPointsV;

    sfnEigen_->B.resize(Mr.rows() + smBendinga.rows(), 1);
    sfnEigen_->B.setZero(Mr.rows() + smBendinga.rows(), 1);
    sfnEigen_->X.resize(_NumberOfControlPointsU * _NumberOfControlPointsV, 1);
    sfnEigen_->X.setZero(_NumberOfControlPointsU * _NumberOfControlPointsV, 1);
    sfnEigen_->LinearSystem.resize(Mr.rows() + smBendinga.rows(), Mr.cols());
    sfnEigen_->LinearSystem.block(0, 0, Mr.rows(), Mr.cols()) = Mr;
    sfnEigen_->LinearSystem.block(Mr.rows(), 0, smBendinga.rows(), smBendinga.cols()) =
        Eigen::MatrixXd(smBendinga);
  }
  /*************
     * estimate(). Estimate a surface from the normals registered in
     * refKf_
     * ********/
  bool ShapeFromNormals::estimate()
  {
    int Mrows_aux = Mrows_;
    int Browsaux = _NumberOfControlPointsU * _NumberOfControlPointsV;
    DefKeyFrame *refdefkf = static_cast<DefKeyFrame *>(refKf_);

    Eigen::MatrixXd Bauxaux(Mrows_aux + Browsaux + 1, 1);
    double MeanDepth(static_cast<DefKeyFrame *>(this->refKf_)
                         ->accMean);

    Bauxaux << sfnEigen_->B.block(0, 0, Mrows_aux + Browsaux, 1),
        Browsaux * MeanDepth;
    Eigen::MatrixXd LynearAuxaux(Mrows_aux + Browsaux + 1, Browsaux);
    LynearAuxaux << sfnEigen_->LinearSystem.block(0, 0, Mrows_aux + Browsaux,
                                                  Browsaux),
        Eigen::MatrixXd::Ones(1, Browsaux);

    Eigen::MatrixXd Solauxs = LynearAuxaux.householderQr().solve(Bauxaux);

    DefKeyFrame *defrefKf_ =
        static_cast<DefKeyFrame *>(this->refKf_);

    for (uint var = 0; var < this->refKf_->mvKeysUn.size(); ++var)
    {
      u_vector_.push_back(defrefKf_->mpKeypointNorm[var].pt.x);
      v_vector_.push_back(defrefKf_->mpKeypointNorm[var].pt.y);
    }

    if (u_vector_.size() == 0)
      return false;

    for (uint i(0); i < Solauxs.rows(); i++)
      if (std::isnan(Solauxs(i)))
      {
        std::cout << "nan fail" << std::endl;
        return false;
      }
    for (uint i(0); i < Solauxs.rows(); i++)
      if (std::isinf(Solauxs(i)))
      {
        std::cout << "inf fail" << std::endl;
        return false;
      }

    double *Array;
    Array = new double[_NumberOfControlPointsU * _NumberOfControlPointsV];

    std::vector<float> depthVec;
    depthVec.reserve(Solauxs.rows());
    for (uint i(0); i < _NumberOfControlPointsU * _NumberOfControlPointsV; i++)
    {
      depthVec.push_back(Solauxs(i));
    }
    std::sort(depthVec.begin(), depthVec.end());
    float corr = 1 / (depthVec[depthVec.size() / 2]);
    sfnEigen_->Solaux.resize(_NumberOfControlPointsU * _NumberOfControlPointsV,
                             1);
    for (uint i(0); i < _NumberOfControlPointsU * _NumberOfControlPointsV; i++)
    {
      Array[i] = corr * Solauxs(i);
      sfnEigen_->Solaux(i) = corr * Solauxs(i);
    }

    BBS::bbs_t bbs;
    bbs.umin = refdefkf->umin;
    bbs.umax = refdefkf->umax;
    bbs.nptsu = _NumberOfControlPointsU;
    bbs.nptsv = _NumberOfControlPointsV;
    bbs.vmax = refdefkf->vmax;
    bbs.vmin = refdefkf->vmin;
    bbs.valdim = 1;

    Eigen::MatrixXd Val(u_vector_.size(), 1);

    BBS::EvalEigen(&bbs, static_cast<double *>(Array), &u_vector_[0], &v_vector_[0],
                   u_vector_.size(), Val, 0, 0);

    for (uint i(0); i < Val.rows(); i++)
    {
      cv::Vec3f X3d;
      X3d(0) = u_vector_[i] * Val(i, 0);
      X3d(1) = v_vector_[i] * Val(i, 0);
      X3d(2) = Val(i, 0);
      defrefKf_->surface->set3DSurfacePoint(i, X3d);
    }

    defrefKf_->surface->saveArray(Array, bbs);

    delete[] Array;
    return true;
  }

  /*************
  * obtainM(). Estimate M that is a matrix that makes the cross
  * product with the current control points and penalizes the 
  * surface estimated if does not accomplish with the normals.
  * ********/
  void ShapeFromNormals::obtainM(BBS::bbs_t &bbs, KeyFrame *Refkf,
                                 Eigen::MatrixXd &M)
  {
    DefKeyFrame *refdefkf = static_cast<DefKeyFrame *>(Refkf);
    Surface *sref = (refdefkf->surface);
    std::vector<double> u_vector_;
    std::vector<double> v_vector_;

    std::vector<cv::Vec3f> Normals;
    u_vector_.reserve(Refkf->mvKeysUn.size());
    v_vector_.reserve(Refkf->mvKeysUn.size());
    Normals.reserve(Refkf->mvKeysUn.size());
    for (uint var = 0; var < Refkf->mvKeysUn.size(); ++var)
    {
      cv::Vec3f Normal;
      auto mp = Refkf->GetMapPoint(var);
      if (!mp)
        continue;
      if (mp->isBad())
        continue;

      if (sref->getNormalSurfacePoint(var, Normal))
      {
        if (Normal == Normal)
        {
          // if (!mp->covNorm)
          //  continue;
          Normals.push_back(Normal);
          u_vector_.push_back(refdefkf->mpKeypointNorm[var].pt.x);
          v_vector_.push_back(refdefkf->mpKeypointNorm[var].pt.y);
        }
      }
    }

    double u[u_vector_.size()];
    std::copy(u_vector_.begin(), u_vector_.end(), u);
    double v[v_vector_.size()];
    std::copy(v_vector_.begin(), v_vector_.end(), v);

    Eigen::MatrixXd coloc = Eigen::MatrixXd::Zero(
        u_vector_.size(), _NumberOfControlPointsU * _NumberOfControlPointsV);
    BBS::colocEigen(&bbs, u, v, u_vector_.size(), coloc);

    Eigen::MatrixXd coloc_du = Eigen::MatrixXd::Zero(
        u_vector_.size(), _NumberOfControlPointsU * _NumberOfControlPointsV);
    Eigen::MatrixXd coloc_dv = Eigen::MatrixXd::Zero(
        u_vector_.size(), _NumberOfControlPointsU * _NumberOfControlPointsV);

    BBS::coloc_derivEigen(&bbs, u, v, u_vector_.size(), 1, 0, coloc_du);
    BBS::coloc_derivEigen(&bbs, u, v, u_vector_.size(), 0, 1, coloc_dv);

    uint npts(u_vector_.size());
    uint nparams(_NumberOfControlPointsU * _NumberOfControlPointsV);
    Eigen::MatrixXd Mdense = Eigen::MatrixXd::Zero(2 * npts, nparams);
    M.resize(2 * npts, nparams);

    for (uint i(0); i < npts; i++)
    {
      Eigen::Vector3d n;
      n << Normals[i](0), Normals[i](1), Normals[i](2);
      n.normalize();
      Eigen::Vector3d etat, etatu, etatv;
      etat << u[i], v[i], 1;
      etatu << 1, 0, 0;
      etatv << 0, 1, 0;
      //  std::cout << "Normal : " <<n << std::endl;
      Eigen::MatrixXd sau(3, coloc_du.cols());
      sau << etat * coloc_du.row(i);

      Eigen::MatrixXd sa2(3, coloc_du.cols());
      sa2 << etatu * coloc.row(i);
      // std::cout << "sa2" << sa2 << std::endl;

      Eigen::MatrixXd Mi(2, coloc_du.cols());
      Mi << n.transpose() * (sau + sa2),
          (n.transpose() * (etat * coloc_dv.row(i) + etatv * coloc.row(i)));
      Eigen::MatrixXd Mi2(2, coloc_du.cols());
      Mi2 << Mi;
      Mdense.row(i) = Mi.row(0);
      Mdense.row(i + npts) = Mi.row(1);
    }
    M = Mdense;
  }
} // namespace defSLAM
