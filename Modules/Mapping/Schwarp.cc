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

#include "Schwarp.h"
#include <Eigen/SparseCore>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "Eigen/SparseCholesky"
#include "thread"
#include <chrono>
#include <ctime>
#include <iostream>
#include <Thirdparty/BBS/bbs_coloc.h>

namespace Warps
{

  Warp::Warp(const std::vector<cv::KeyPoint> &KP1,
             const std::vector<cv::KeyPoint> &KP2, std::vector<float> &invSigmas,
             double umin, double umax, double vmin, double vmax, int NCu, int NCv,
             int valdim, double fx, double fy)
      : CostFunction(), KP1(KP1), KP2(KP2), invSigmas(invSigmas), fx(fx), fy(fy)
  {
    bbs.umin = umin;
    bbs.umax = umax;
    bbs.nptsu = NCu;
    bbs.vmin = vmin;
    bbs.vmax = vmax;
    bbs.nptsv = NCv;
    bbs.valdim = valdim;
    double u[KP1.size()];
    double v[KP1.size()];

    Eigen::MatrixXd q2;
    q2.resize(KP1.size(), 2);
    uint i(0);
    Eigen::MatrixXd q1;
    q1.resize(KP2.size(), 2);
    for (std::vector<cv::KeyPoint>::const_iterator it = KP1.begin();
         it != KP1.end(); it++)
    {
      u[i] = it->pt.x;
      v[i] = it->pt.y;
      q1(i, 0) = it->pt.x;
      q1(i, 0) = it->pt.y;
      q2(i, 0) = KP2[i].pt.x;
      q2(i, 1) = KP2[i].pt.y;
      i++;
    }

    Eigen::MatrixXd ColocGen = Eigen::MatrixXd::Zero(
        KP1.size(), _NumberOfControlPointsU * _NumberOfControlPointsV);
    BBS::colocEigen(&bbs, u, v, KP1.size(), ColocGen);

    Jdata.resize(KP1.size() * 2 * _NumberOfControlPointsU *
                     _NumberOfControlPointsV * 2,
                 0.0);
    Eigen::MatrixXd J(ColocGen.rows() * 2, ColocGen.cols() * 2);

    J << -ColocGen * double(fx),
        Eigen::MatrixXd::Zero(ColocGen.rows(), ColocGen.cols()),
        Eigen::MatrixXd::Zero(ColocGen.rows(), ColocGen.cols()),
        -ColocGen * double(fy);

    for (uint i(0); i < KP1.size() * 2; i++)
    {
      for (uint j(0); j < _NumberOfControlPointsU * _NumberOfControlPointsV * 2;
           j++)
      {
        Jdata[i * _NumberOfControlPointsU * _NumberOfControlPointsV * 2 + j] =
            J(i, j);
      }
    }
    this->set_num_residuals(KP1.size() * 2);
    mutable_parameter_block_sizes()->push_back(_NumberOfControlPointsU *
                                               _NumberOfControlPointsV * 2);
  }

  void Warp::initialize(
      std::vector<cv::KeyPoint> &KP1, std::vector<cv::KeyPoint> &KP2, double lambda,
      double umin, double umax, double vmin, double vmax, int NCu, int NCv,
      int valdim,
      double (&x)[_NumberOfControlPointsU * _NumberOfControlPointsV * 2])
  {
    BBS::bbs_t bbs;

    bbs.umin = umin;
    bbs.umax = umax;
    bbs.nptsu = NCu;
    bbs.vmin = vmin;
    bbs.vmax = vmax;
    bbs.nptsv = NCv;
    bbs.valdim = valdim;
    double u[KP1.size()];
    double v[KP1.size()];

    Eigen::MatrixXd q2;
    q2.resize(KP1.size(), 2);
    uint i(0);
    Eigen::MatrixXd q1;
    q1.resize(KP2.size(), 2);
    for (std::vector<cv::KeyPoint>::iterator it = KP1.begin(); it != KP1.end();
         it++)
    {
      u[i] = it->pt.x;
      v[i] = it->pt.y;
      q1(i, 0) = it->pt.x;
      q1(i, 0) = it->pt.y;
      q2(i, 0) = KP2[i].pt.x;
      q2(i, 1) = KP2[i].pt.y;
      i++;
    }

    Eigen::SparseMatrix<double> smcoloc(KP1.size(), (NCu * NCv));
    Eigen::SparseMatrix<double> smcoloctsmcoloc((NCu * NCv), KP1.size());

    Eigen::MatrixXd ColocGen = Eigen::MatrixXd::Zero(
        KP1.size(), _NumberOfControlPointsU * _NumberOfControlPointsV);
    BBS::colocEigen(&bbs, u, v, KP1.size(), ColocGen);

    smcoloc = ColocGen.sparseView();
    smcoloctsmcoloc = smcoloc.transpose() * smcoloc;

    Eigen::SparseMatrix<double> smBending((NCu * NCv), (NCu * NCv));
    BBS::BendingEigen(&bbs, lambda, smBending);

    Eigen::SparseMatrix<double> cptsaux = smcoloctsmcoloc + smBending;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(cptsaux);
    Eigen::SparseMatrix<double> I((NCu * NCv), (NCu * NCv));
    I.setIdentity();
    Eigen::SparseMatrix<double> cptsaux_inv = solver.solve(I);
    Eigen::SparseMatrix<double> trans = smcoloc.transpose();
    Eigen::SparseMatrix<double> cpts2 = cptsaux_inv * trans;
    Eigen::Map<
        Eigen::Matrix<double, _NumberOfControlPointsU * _NumberOfControlPointsV,
                      2, Eigen::ColMajor>>
        a(x);
    a = Eigen::MatrixXd(cpts2) * q2;
  }

  std::vector<cv::KeyPoint> Warp::getEstimates(
      std::vector<cv::KeyPoint> &KP1, double umin, double umax, double vmin,
      double vmax, int NCu, int NCv, int valdim,
      double (&x)[_NumberOfControlPointsU * _NumberOfControlPointsV * 2],
      Eigen::MatrixXd &Val, uint a1, uint b)
  {
    BBS::bbs_t bbs;
    bbs.umin = umin;
    bbs.umax = umax;
    bbs.nptsu = _NumberOfControlPointsU;
    bbs.vmin = vmin;
    bbs.vmax = vmax;
    bbs.nptsv = _NumberOfControlPointsV;
    bbs.valdim = valdim;

    Eigen::Map<const Eigen::Matrix<
                   double, _NumberOfControlPointsU * _NumberOfControlPointsV, 2>,
               Eigen::ColMajor>
        ControlPoints(x);

    double Array[ControlPoints.rows() * 2];
    double u[KP1.size()];
    double v[KP1.size()];
    for (int n = 0; n < 2; ++n)
    {
      for (int l = 0; l < ControlPoints.rows(); ++l)
      {
        Array[bbs.valdim * l + n] =
            ControlPoints(l, n); // double(ControlPoints(l,n)); //A matrix is
                                 // filled/accessed as a linear array.
      }
    }
    ///// REPROJECTION ERROR

    for (uint i(0); i < KP1.size(); i++)
    {
      u[i] = KP1[i].pt.x;
      v[i] = KP1[i].pt.y;
    }

    Eigen::MatrixXd val;
    val.resize(KP1.size(), 2);

    BBS::EvalEigen(&bbs, static_cast<double *>(Array), u, v, KP1.size(), val, a1,
                   b);

    std::vector<cv::KeyPoint> a;
    for (uint i(0); i < KP1.size(); i++)
    {
      cv::KeyPoint kpoint;
      kpoint.pt.x = (val(i, 0));
      kpoint.pt.y = (val(i, 1));
      a.push_back(kpoint);
    }

    double X[_NumberOfControlPointsU * _NumberOfControlPointsV];
    double Y[_NumberOfControlPointsU * _NumberOfControlPointsV];
    uint us(0);
    for (int i(0); i < bbs.nptsu; i++)
    {
      for (int j(0); j < bbs.nptsv; j++)
      {
        X[us] = double((bbs.umax - bbs.umin) * i) / (bbs.nptsu - 1) + bbs.umin;
        Y[us] = double((bbs.vmax - bbs.vmin) * j) / (bbs.nptsv - 1) + bbs.vmin;
        us++;
      }
    }

    BBS::EvalEigen(&bbs, static_cast<double *>(Array), X, Y,
                   _NumberOfControlPointsU * _NumberOfControlPointsV, Val, 0, 0);

    return a;
  }

  bool Warp::Evaluate(double const *const *parameters, double *residuals,
                      double **jacobians) const
  {
    {
      const uint rows(_NumberOfControlPointsU * _NumberOfControlPointsV);
      Eigen::Map<
          const Eigen::Matrix<
              double, _NumberOfControlPointsU * _NumberOfControlPointsV, 2>,
          Eigen::RowMajor>
          ControlPoints(parameters[0]);

      BBS::bbs_t b(bbs);
      double Array[ControlPoints.rows() * 2];
      double u[KP1.size()];
      double v[KP1.size()];

      // Eigen::Matrix<double,_NumberOfControlPointsU*_NumberOfControlPointsV,2>
      // CP = ControlPoints.cast<double>();
      for (int n = 0; n < 2; ++n)
      {
        for (int l = 0; l < ControlPoints.rows(); ++l)
        {
          Array[bbs.valdim * l + n] =
              ControlPoints(l, n); // double(ControlPoints(l,n)); //A matrix is
                                   // filled/accessed as a linear array.
        }
      }

      ///// REPROJECTION ERROR
      for (uint i(0); i < KP1.size(); i++)
      {
        u[i] = KP1[i].pt.x;
        v[i] = KP1[i].pt.y;
      }

      Eigen::MatrixXd val;
      val.resize(KP1.size(), 2);

      BBS::EvalEigen(&b, (Array), u, v, KP1.size(), val, 0, 0);

      for (uint i(0); i < KP1.size(); i++)
      {
        residuals[i] = invSigmas[i] * (double(KP2[i].pt.x) - val(i, 0)) * fx;
        residuals[i + KP1.size()] =
            invSigmas[i] * (double(KP2[i].pt.y) - val(i, 1)) * fy;
      }

      if (jacobians != NULL && jacobians[0] != NULL)
      {
        for (uint i(0); i < KP2.size() * 2; i++)
        {
          for (uint j(0); j < rows * 2; j++)
          {
            jacobians[0][i * rows * 2 + j] = Jdata[i * rows * 2 + j];
          }
        }
        for (uint i(0); i < KP2.size(); i++)
        {
          for (uint j(0); j < rows * 2; j++)
          {
            jacobians[0][i * rows * 2 + j] = Jdata[i * rows * 2 + j];
            jacobians[0][(i + KP2.size()) * rows * 2 + j] =
                Jdata[i * rows * 2 + j];
          }
        }
      }
      return true;
    }
  }

  Schwarzian::Schwarzian(double lambda, double umin, double umax, double vmin,
                         double vmax, int valdim)
      : CostFunction(), lambda(lambda)
  {
    /******** 
     * Regularization for the nodes of the mesh. 
     * **********/
    this->bbs.umin = umin;
    this->bbs.umax = umax;
    this->bbs.nptsu = _NumberOfControlPointsU;
    this->bbs.vmin = vmin;
    this->bbs.vmax = vmax;
    this->bbs.nptsv = _NumberOfControlPointsV;
    this->bbs.valdim = valdim;
    uint rows(_NumberOfControlPointsU * _NumberOfControlPointsV);
    double X[rows];
    double Y[rows];
    uint us(0);
    for (int i(0); i < bbs.nptsu; i++)
    {
      for (int j(0); j < bbs.nptsv; j++)
      {
        X[us] = double((bbs.umax - bbs.umin) * i) / (bbs.nptsu - 1) + bbs.umin;
        Y[us] = double((bbs.vmax - bbs.vmin) * j) / (bbs.nptsv - 1) + bbs.vmin;
        us++;
      }
    }

    Eigen::MatrixXd coloc_duu1 = Eigen::MatrixXd::Zero(rows, rows);
    Eigen::MatrixXd coloc_duv1 = Eigen::MatrixXd::Zero(rows, rows);
    Eigen::MatrixXd coloc_dvv1 = Eigen::MatrixXd::Zero(rows, rows);

    Eigen::MatrixXd coloc_du1 = Eigen::MatrixXd::Zero(rows, rows);
    Eigen::MatrixXd coloc_dv1 = Eigen::MatrixXd::Zero(rows, rows);
    BBS::coloc_derivEigen(&this->bbs, X, Y, rows, 1, 1, coloc_duv1);
    BBS::coloc_derivEigen(&this->bbs, X, Y, rows, 0, 2, coloc_dvv1);
    BBS::coloc_derivEigen(&this->bbs, X, Y, rows, 2, 0, coloc_duu1);
    Eigen::Map<Eigen::Matrix<
        double, 4 * _NumberOfControlPointsU * _NumberOfControlPointsV,
        _NumberOfControlPointsU * _NumberOfControlPointsV * 3>>
        A21(A21v);
    Eigen::Map<Eigen::Matrix<
        double, 6 * _NumberOfControlPointsU * _NumberOfControlPointsV,
        _NumberOfControlPointsU * _NumberOfControlPointsV * 2>>
        B21(B21v);
    A21 << coloc_duu1, coloc_duv1, coloc_dvv1, coloc_duu1, coloc_duv1, coloc_dvv1,
        coloc_duu1, coloc_duv1, coloc_dvv1, coloc_duu1, coloc_duv1, coloc_dvv1;

    BBS::coloc_derivEigen(&this->bbs, X, Y, rows, 0, 1, coloc_dv1);
    BBS::coloc_derivEigen(&this->bbs, X, Y, rows, 1, 0, coloc_du1);

    B21 << coloc_du1, Eigen::MatrixXd(coloc_dv1), Eigen::MatrixXd(coloc_du1),
        Eigen::MatrixXd(coloc_dv1), Eigen::MatrixXd(coloc_du1),
        Eigen::MatrixXd(coloc_dv1), Eigen::MatrixXd(coloc_du1),
        Eigen::MatrixXd(coloc_dv1), Eigen::MatrixXd(coloc_du1),
        Eigen::MatrixXd(coloc_dv1), Eigen::MatrixXd(coloc_du1),
        Eigen::MatrixXd(coloc_dv1);
    this->set_num_residuals(_NumberOfControlPointsU * _NumberOfControlPointsV *
                            4);
    mutable_parameter_block_sizes()->push_back(_NumberOfControlPointsU *
                                               _NumberOfControlPointsV * 2);
  }

  bool Schwarzian::Evaluate(double const *const *parameters, double *residuals,
                            double **jacobians) const
  {
    {
      const uint rows(_NumberOfControlPointsU * _NumberOfControlPointsV);

      Eigen::Map<
          const Eigen::Matrix<
              double, _NumberOfControlPointsU * _NumberOfControlPointsV, 2>,
          Eigen::RowMajor>
          ControlPoints(parameters[0]);

      BBS::bbs_t b(bbs);
      double Array[ControlPoints.rows() * 2];

      for (int n = 0; n < 2; ++n)
      {
        for (int l = 0; l < ControlPoints.rows(); ++l)
        {
          Array[bbs.valdim * l + n] =
              ControlPoints(l, n); // double(ControlPoints(l,n)); //A matrix is
                                   // filled/accessed as a linear array.
        }
      }

      Eigen::Map<Eigen::Matrix<
          double, _NumberOfControlPointsU * _NumberOfControlPointsV * 4, 1>>
          ResidualSchwarzian(&residuals[0]);

      double X[rows];
      double Y[rows];
      uint us(0);
      for (int i(0); i < bbs.nptsu; i++)
      {
        for (int j(0); j < bbs.nptsv; j++)
        {
          X[us] = double((bbs.umax - bbs.umin) * i) / (bbs.nptsu - 1) + bbs.umin;
          Y[us] = double((bbs.vmax - bbs.vmin) * j) / (bbs.nptsv - 1) + bbs.vmin;
          us++;
        }
      }

      Eigen::MatrixXd xydu1(rows, 2), xydv1(rows, 2),
          xydu2(rows, 2), xydv2(rows, 2), xydudv(rows, 2); // First column contain dx
                                                           // second column dy

      BBS::EvalEigen(&b, Array, X, Y, rows, xydu1, 1, 0);
      BBS::EvalEigen(&b, Array, X, Y, rows, xydv1, 0, 1);

      BBS::EvalEigen(&b, Array, X, Y, rows, xydu2, 2, 0);
      BBS::EvalEigen(&b, Array, X, Y, rows, xydv2, 0, 2);
      BBS::EvalEigen(&b, Array, X, Y, rows, xydudv, 1, 1);

      // Schwarzian equations:
      // (d2x/du2*dy/du-d2y/dudu*dx/du)
      ResidualSchwarzian.block(0, 0, rows, 1) =
          ((xydu2.block(0, 0, rows, 1)
                .cwiseProduct(xydu1.block(0, 1, rows, 1)) -
            xydu2.block(0, 1, rows, 1)
                .cwiseProduct(xydu1.block(0, 0, rows, 1)))) *
          lambda;
      // d2y/dvdv*dx/dv-d2x/dv2*dy/dv
      ResidualSchwarzian.block(rows, 0, rows, 1) =
          ((xydv2.block(0, 1, rows, 1)
                .cwiseProduct(xydv1.block(0, 0, rows, 1)) -
            xydv2.block(0, 0, rows, 1)
                .cwiseProduct(xydv1.block(0, 1, rows, 1)))) *
          lambda;

      // (d2x/dudu*dy/dv-d2y/du2*dx/dv)+2*(dx/dudv*dy/du-dy/dudv*dx/du)
      ResidualSchwarzian.block(rows * 2, 0, rows, 1) =
          ((xydu2.block(0, 0, rows, 1)
                .cwiseProduct(xydv1.block(0, 1, rows, 1)) -
            xydu2.block(0, 1, rows, 1)
                .cwiseProduct(xydv1.block(0, 0, rows, 1)) +
            2 * (xydudv.block(0, 0, rows, 1)
                     .cwiseProduct(xydu1.block(0, 1, rows, 1)) -
                 xydudv.block(0, 1, rows, 1)
                     .cwiseProduct(xydu1.block(0, 0, rows, 1))))) *
          lambda;

      // (d2x/dvdv*dy/du-d2y/dv2*dx/du)+2*(dx/dudv*dy/dv-dy/dudv*dx/dv)
      ResidualSchwarzian.block(rows * 3, 0, rows, 1) =
          ((xydv2.block(0, 1, rows, 1)
                .cwiseProduct(xydu1.block(0, 0, rows, 1)) -
            xydv2.block(0, 0, rows, 1)
                .cwiseProduct(xydu1.block(0, 1, rows, 1)) +
            2 * (xydudv.block(0, 1, rows, 1)
                     .cwiseProduct(xydv1.block(0, 0, rows, 1)) -
                 xydudv.block(0, 0, rows, 1)
                     .cwiseProduct(xydv1.block(0, 1, rows, 1))))) *
          lambda;

      if (jacobians != NULL && jacobians[0] != NULL)
      {
        Eigen::MatrixXd B(6 * rows, rows * 2), A(4 * rows, rows * 3);

        Eigen::Map<const Eigen::Matrix<
            double, 4 * _NumberOfControlPointsU * _NumberOfControlPointsV,
            _NumberOfControlPointsU * _NumberOfControlPointsV * 3>>
            A21(A21v);
        Eigen::Map<const Eigen::Matrix<
            double, 6 * _NumberOfControlPointsU * _NumberOfControlPointsV,
            _NumberOfControlPointsU * _NumberOfControlPointsV * 2>>
            B21(B21v);
        std::thread j1([&xydu1, &xydv1, &A, &A21] {
          Eigen::Matrix<double,
                        _NumberOfControlPointsU * _NumberOfControlPointsV * 4, 1>
              A1;
          A1 << xydu1.col(0), xydu1.col(1), xydv1.col(0),
              xydv1.col(1);
          A = A1.asDiagonal() * A21;
        });
        std::thread j2([&xydu2, &xydudv, &xydv2, &B, &B21] {
          Eigen::MatrixXd B1(6 * rows, 1);
          B1 << xydu2.block(0, 0, rows, 1),
              xydu2.block(0, 1, rows, 1), xydudv.block(0, 0, rows, 1),
              xydudv.block(0, 1, rows, 1), xydv2.block(0, 0, rows, 1),
              xydv2.block(0, 1, rows, 1);
          B = B1.asDiagonal() * B21;
        });

        int r = rows;
        int c = rows;
        Eigen::MatrixXd jI, jJ, jM, jN;
        jI.resize(r, c * 2);
        jJ.resize(r, c * 2);
        jM.resize(r, c * 2);
        jN.resize(r, c * 2);
        j1.join();
        j2.join();

        j1 = std::thread([&jI, &jJ, &A, &B, &r, &c] {
          jI << (A.block(r, 0, c, c)) - B.block(r, 0, c, c),
              B.block(0, 0, c, c) - A.block(0, 0, c, c);
          jJ << B.block(5 * r, c, c, c) - A.block(3 * r, 2 * c, c, c),
              A.block(2 * r, 2 * c, c, c) - (B.block(4 * r, c, c, c));
        });
        j2 = std::thread([&jM, &jN, &A, &B, &r, &c] {
          jM << A.block(3 * r, 0, c, c) - (B.block(r, c, c, c)) +
                    2 * A.block(r, c, c, c) - 2 * (B.block(3 * r, 0, c, c)),
              B.block(0, c, c, c) - (A.block(2 * r, 0, c, c)) +
                  2 * (B.block(2 * r, 0, c, c)) - 2 * (A.block(0, c, c, c));
          jN << (B.block(5 * r, 0, c, c)) - (A.block(r, 2 * c, c, c)) -
                    2 * (A.block(3 * r, c, c, c)) + 2 * (B.block(3 * r, c, c, c)),
              (A.block(0, 2 * c, c, c)) - (B.block(4 * r, 0, c, c)) -
                  2 * (B.block(2 * r, c, c, c)) + 2 * A.block(2 * r, c, c, c);
        });

        j1.join();
        j2.join();

        Eigen::MatrixXd Jacobians;
        Jacobians.resize(jI.rows() * 4, jI.cols());
        Jacobians << lambda * jI, lambda * jJ, lambda * jM, lambda * jN;

        Eigen::SparseMatrix<double> S = Jacobians.sparseView();
        for (uint i(0); i < 4 * rows; i++)
        {
          for (uint j(0); j < rows * 2; j++)
          {
            jacobians[0][i * rows * 2 + j] = 0;
          }
        }
        for (int k = 0; k < S.outerSize(); ++k)
        {
          for (Eigen::SparseMatrix<double>::InnerIterator it(S, k); it; ++it)
          {
            jacobians[0][it.row() * (Jacobians.cols()) + it.col()] = it.value();
          }
        }
      }

      return true;
    }
  }
} // namespace Warps
