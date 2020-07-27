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

#include "NormalEstimator.h"
#include "DefKeyFrame.h"
#define MAX_LOG_LEVEL 1

namespace defSLAM
{
  // Construtor
  NormalEstimator::NormalEstimator(WarpDatabase *WarpDB) : WarpDB(WarpDB){};

  // Destructor
  NormalEstimator::~NormalEstimator() { WarpDB = nullptr; }

  /***********
  * Core function to estimate K1 and K2 that are the two first components of the normal.
  * n = [k1, k2 , 1- k1*u-k2*v]
  * *********/
  void NormalEstimator::ObtainK1K2()
  {
    // This function iterates over all the map points checking if they
    // have new observation to reestimate their normals.

    uint nummin(1);
    auto &diffDB = WarpDB->getDiffDatabase();
    auto &toProccess = WarpDB->getToProccess();
    uint counter(0);
    uint counterPoints(0);

    for (auto &process : toProccess) // Iteration over map points in WarpDB.
    {
      // Only process map point if there is a new observation
      if (!process.second)
        continue;
      process.second = false;

      MapPoint *mapPoint = process.first;
      if (!mapPoint)
        continue;
      if (mapPoint->isBad())
        continue;

      auto &kf2kfDiff = diffDB[process.first];
      if ((kf2kfDiff.size()) < nummin)
        continue;
      counterPoints++;
      ceres::Problem problem;

      double x[2];

      KeyFrame *refKF = mapPoint->GetReferenceKeyFrame();
      size_t idx = mapPoint->GetIndexInKeyFrame(refKF);
      if (idx < 0)
        continue;

      int count(0);
      // Set coefficient of the polynomials
      for (auto &diffkf2kf : kf2kfDiff)
      {
        if (refKF != diffkf2kf->KFToKF.first)
          continue;
        idx = diffkf2kf->idx1;
        float a, b, c, d, t1, t2, I1u, I1v, I2u, I2v, e1, e2;
        /***********
   *  Jacobian Image Points
   *  J12 = [[J12a,J12c],
   *         [J12b,J12d]]
   *  J21 = [[J21a,J21c],
   *         [J21b,J21d]]
   ***********/
        a = diffkf2kf->J12a;
        b = diffkf2kf->J12b;
        c = diffkf2kf->J12c;
        d = diffkf2kf->J12d;
        t1 = -diffkf2kf->J12b * diffkf2kf->H12vvx / 2 +
             diffkf2kf->J12a * diffkf2kf->H12vvy / 2;
        t2 = -(diffkf2kf->J12d * diffkf2kf->H12vvx) / 2 +
             (diffkf2kf->J12c * diffkf2kf->H12vvy) / 2;
        I1u = diffkf2kf->I1u;
        I1v = diffkf2kf->I1v;
        I2u = diffkf2kf->I2u;
        I2v = diffkf2kf->I2v;
        e2 = 1 + I2u * I2u + I2v * I2v;
        e1 = 1 + I1u * I1u + I1v * I1v;
        double eq1[10];
        double eq2[10];

        PolySolver::getCoefficients(a, b, c, d, t1, t2, e1, e2, I1u, I1v, I2u,
                                    I2v, 0, eq1);
        PolySolver::getCoefficients(a, b, c, d, t1, t2, e1, e2, I1u, I1v, I2u,
                                    I2v, 1, eq2);

        Eigen::Map<Eigen::Matrix<double, 1, 10>> eqM1(eq1);
        Eigen::Map<Eigen::Matrix<double, 1, 10>> eqM2(eq2);

        ceres::CostFunction *cost_function = new PolySolver(eqM1, eqM2);
        problem.AddResidualBlock(cost_function, nullptr, x);
        count++;
      }
      // Solve polynomials system with an optimization.
      if (count > 0) // Check that the point has equations to solve
      {
        cv::Vec3f Ni;
        // The initial solution is the last solution estimated if it exists,
        // and zero if it does not.
        if (static_cast<DefKeyFrame *>(refKF)
                ->surface->getNormalSurfacePoint(idx, Ni))
        {
          x[0] = Ni(0);
          x[1] = Ni(1);
        }
        else
        {
          x[0] = 0.0;
          x[1] = -0.0;
        }

        problem.AddParameterBlock(x, 2);
        // Solve it with ceres.
        ceres::Solver::Options options;
        options.check_gradients = false;
        options.gradient_tolerance = 1E-8;
        options.logging_type = ceres::PER_MINIMIZER_ITERATION;
        options.linear_solver_type = ceres::DENSE_NORMAL_CHOLESKY;
        options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
        options.num_threads = 1;
        options.function_tolerance = 1e-10;
        options.minimizer_progress_to_stdout = false;
        options.max_num_iterations = 200;
        ceres::Solver::Summary summary;

        Solve(options, &problem, &summary);

        ceres::Covariance::Options optionscov;
        ceres::Covariance covariance(optionscov);
        vector<pair<const double *, const double *>> covariance_blocks;
        covariance_blocks.push_back(make_pair(x, x));
        if (!covariance.Compute(covariance_blocks, &problem))
          continue;
        covariance.GetCovarianceBlock(x, x, mapPoint->covNorm);
        // Recover the normal in the rest of keyframes
        DefKeyFrame *defKf = static_cast<DefKeyFrame *>(refKF);
        float I1u = defKf->mpKeypointNorm[idx].pt.x;
        float I1v = defKf->mpKeypointNorm[idx].pt.y;
        cv::Vec3f normal;
        // Normal equation cross multiplying the vectors in the tangent space.
        normal(0) = x[0];
        normal(1) = x[1];
        normal(2) = 1 - x[0] * I1u - x[1] * I1v;
        defKf->surface->setNormalSurfacePoint(idx, normal);
      }

      // Once estimated we can calculate the solution in all the keyframes.
      for (auto &diffkf2kf : kf2kfDiff)
      {
        size_t idx1 = diffkf2kf->idx1;
        size_t idx2 = diffkf2kf->idx2;

        cv::Vec2d norm;
        if (refKF == diffkf2kf->KFToKF.first)
        {
          norm(0) = x[0];
          norm(1) = x[1];
        }
        else
        {
          cv::Vec3f Ni;
          if (static_cast<DefKeyFrame *>(diffkf2kf->KFToKF.first)
                  ->surface->getNormalSurfacePoint(idx1, Ni))
          {
            norm(0) = Ni(0);
            norm(1) = Ni(1);
          }
          else
          {
            continue;
          }
        }

        float j21_11, j21_21, j21_12, j21_22, t1, t2, a, b, c, d;
        j21_11 = diffkf2kf->J21a;
        j21_12 = diffkf2kf->J21c;
        j21_21 = diffkf2kf->J21b;
        j21_22 = diffkf2kf->J21d;
        a = diffkf2kf->J12a;
        b = diffkf2kf->J12b;
        c = diffkf2kf->J12c;
        d = diffkf2kf->J12d;
        float detJ12 = a * d - c * b;
        t1 = -b * diffkf2kf->H12vvx / 2 + a * diffkf2kf->H12vvy / 2;
        t2 = (d * diffkf2kf->H12uux) / 2 - (c * diffkf2kf->H12uuy) / 2;
        auto k1_kfi = j21_11 * norm(0) + j21_12 * norm(1) +
                      (d * t2 - b * t1) / (detJ12 * detJ12);
        auto k2_kfi = j21_21 * norm(0) + j21_22 * norm(1) +
                      (a * t1 - c * t2) / (detJ12 * detJ12);
        auto I2u = diffkf2kf->I2u;
        auto I2v = diffkf2kf->I2v;
        cv::Vec3f normal_i;
        normal_i(0) = k1_kfi;
        normal_i(1) = k2_kfi;
        normal_i(2) = 1 - k1_kfi * I2u - k2_kfi * I2v;
        DefKeyFrame *defKf_i =
            static_cast<DefKeyFrame *>(diffkf2kf->KFToKF.second);
        defKf_i->surface->setNormalSurfacePoint(idx2, normal_i);
      }
      counter++;
    }
    std::cout << "NORMALS REESTIMATED : " << counter << " - " << counterPoints
              << std::endl;
  } // namespace defSLAM
} // namespace defSLAM
