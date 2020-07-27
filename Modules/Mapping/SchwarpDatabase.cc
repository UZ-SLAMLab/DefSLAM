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

#include <list>
#include <mutex>
#include <set>
#include <thread>
#include <vector>

#include "Frame.h"
#include "KeyFrame.h"

#include "DefKeyFrame.h"
#include "DefORBmatcher.h"

#include "Schwarp.h"
#include "SchwarpDatabase.h"

namespace defSLAM
{
  /*********** 
     * It insert the keyframe to the database. There we 
     * search for its covisible anchor keyframes. Once, 
     * we have them, if they have a certain number of matches
     * it estimates a initial warp between the anchor keyframes 
     * and the new keyframe. With that warp, it search for more 
     * matches. Once it is done it recalculate the warp and then
     * it save the correspondences and its derivatives
     * in mapPointsDB_ (WarpDatabase).
     * ******/
  void SchwarpDatabase::add(KeyFrame *mpCurrentKeyFrame)
  {
    // Check a keyframe is not inserted twice
    CHECK(mpkeyframes.count(mpCurrentKeyFrame) == 0);
    // If it's the first keyframe don't do anything
    if (mpkeyframes.size() == 0)
    {
      mpkeyframes.insert(mpCurrentKeyFrame);
      return;
    }

    // Get matched points
    const vector<MapPoint *> vpMapPointMatches =
        mpCurrentKeyFrame->GetMapPointMatches();

    // Get keyframes of reference for the map points matched
    std::unordered_map<KeyFrame *, int> countKFMatches;

    for (size_t i = 0; i < vpMapPointMatches.size(); i++)
    {
      MapPoint *mapPoint = vpMapPointMatches[i];
      if (!mapPoint)
        continue;
      if (mapPoint->isBad())
        continue;
      // Keyframe in which the point was created.
      KeyFrame *refkf = mapPoint->GetReferenceKeyFrame();
      if (countKFMatches.count(refkf) == 0)
        countKFMatches[refkf] = 0;
      countKFMatches[refkf]++;
    }
    /// We search all the possible matches with the reference keyframes.
    /// there will be more points
    for (const auto &kv : countKFMatches)
    {
      // Only take into account those with more than 30 matches
      KeyFrame *refkf = kv.first;
      // Get matched points between keyframes
      vector<pair<size_t, size_t>> vMatchedIndices;
      for (size_t i = 0; i < vpMapPointMatches.size(); i++)
      {
        MapPoint *mapPoint = vpMapPointMatches[i];
        if (!mapPoint)
          continue;
        if (mapPoint->isBad())
          continue;
        /// Check that the point is in both keyframes
        if (mapPoint->IsInKeyFrame(mpCurrentKeyFrame) &&
            (mapPoint->IsInKeyFrame(refkf)))
        {
          const int idx1 = mapPoint->GetIndexInKeyFrame(refkf);
          const int idx2 = mapPoint->GetIndexInKeyFrame(mpCurrentKeyFrame);
          vMatchedIndices.push_back({idx1, idx2});
        }
      }
      if (vMatchedIndices.size() < 20)
        continue;
      // Create Schwarp
      double x[_NumberOfControlPointsU * _NumberOfControlPointsV * 2];
      DefORBmatcher Matcher;
      std::cout << "Finding by Schwarp" << std::endl;
      Matcher.findbyWarp(refkf, mpCurrentKeyFrame, vMatchedIndices, x, lambda_);
      if (vMatchedIndices.size() < 20)
        continue;
      std::cout << "Calculating final Schwarp" << std::endl;

      this->calculateSchwarps(refkf, mpCurrentKeyFrame, vMatchedIndices, x, lambda_);
      static_cast<DefKeyFrame *>(mpCurrentKeyFrame)->KeyframesRelated++;
    }
    int count(0);
    for (auto mp : newInformation_)
    {
      if (mp.second)
      {
        count++;
      }
    }
    mpkeyframes.insert(mpCurrentKeyFrame);
  }

  // Delete kf from database
  void SchwarpDatabase::erase(KeyFrame *mpCurrentKeyFrame)
  {
    mpkeyframes.erase(mpCurrentKeyFrame);
  }

  // clear the database
  void SchwarpDatabase::clear() { mpkeyframes.clear(); }

  /************************
   *  Function that estimates the warps between two keyframes 
   * given the correspondences vMatchesIndices and 
  * saves their correspondances and its derivatives. It uses 
  * the Schwarzian penalization with a lambda of lambda. 
  * ******************/
  void SchwarpDatabase::calculateSchwarps(
      KeyFrame *KFi, KeyFrame *KF2i,
      vector<pair<size_t, size_t>> &vMatchedIndices,
      double (&x)[_NumberOfControlPointsU * _NumberOfControlPointsV * 2],
      double lambda)
  {
    // Calculate From One To Two
    DefKeyFrame *KF = static_cast<DefKeyFrame *>(KFi);
    DefKeyFrame *KF2 = static_cast<DefKeyFrame *>(KF2i);
    std::vector<cv::KeyPoint> KP1;
    std::vector<cv::KeyPoint> KP2, kp2un;
    std::vector<float> invSigmas;

    KP1.reserve(vMatchedIndices.size());
    KP2.reserve(vMatchedIndices.size());
    kp2un.reserve(vMatchedIndices.size());
    double umin = KF->umin;
    double umax = KF->umax;
    double vmin = KF->vmin;
    double vmax = KF->vmax;

    for (uint ikp = 0; ikp < vMatchedIndices.size(); ikp++)
    {
      const auto &idx1 = vMatchedIndices[ikp].first;
      const auto &idx2 = vMatchedIndices[ikp].second;

      cv::KeyPoint kp1 = KF->mpKeypointNorm[idx1];
      cv::KeyPoint kp2 = KF2->mpKeypointNorm[idx2];
      cv::KeyPoint kp2un_i = KF2->mvKeysUn[idx2];

      KP1.push_back(std::move(kp1));
      KP2.push_back(std::move(kp2));
      kp2un.push_back(std::move(kp2un_i));
      const cv::KeyPoint &kpUn = KF->mvKeysUn[idx1];

      const float &invSigma2 = KF->mvInvLevelSigma2[kpUn.octave];
      invSigmas.push_back(sqrt(invSigma2));
    }

    Eigen::MatrixXd ControlPointsInitialKF(
        _NumberOfControlPointsU * _NumberOfControlPointsV, 2);
    uint us(0);
    for (int i(0); i < KF->NCu; i++)
    {
      for (int j(0); j < KF->NCv; j++)
      {
        ControlPointsInitialKF(us, 0) =
            double((KF->umax - KF->umin) * i) / (KF->NCu - 1) + KF->umin;
        ControlPointsInitialKF(us, 1) =
            double((KF->vmax - KF->vmin) * j) / (KF->NCv - 1) + KF->vmin;
        us++;
      }
    }

    ceres::CostFunction *Rep =
        new Warps::Warp(KP1, KP2, invSigmas, umin, umax, vmin, vmax, KF->NCu,
                        KF->NCv, KF->valdim, double(KF->fy), double(KF->fx));
    /// Analitical implementation
    ceres::CostFunction *Schwarzian =
        new Warps::Schwarzian(lambda, umin, umax, vmin, vmax, KF->valdim);

    ceres::Problem problem;
    problem.AddParameterBlock(x, _NumberOfControlPointsU *
                                     _NumberOfControlPointsV * 2);
    problem.AddResidualBlock(Rep, new ceres::HuberLoss(5.77), x);
    problem.AddResidualBlock(Schwarzian, nullptr, x);

    // Run the solver!
    ceres::Solver::Options options;
    options.dynamic_sparsity = true;
    options.num_threads = 1;
    options.max_num_iterations = 3;
    options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
    options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
    options.minimizer_progress_to_stdout = false;
    ceres::Solver::Summary summary;

    Solve(options, &problem, &summary);

    Eigen::Matrix<double, _NumberOfControlPointsU * _NumberOfControlPointsV, 2>
        ControlPointsinitial;

    us = 0;
    for (int i(0); i < KF->NCu; i++)
    {
      for (int j(0); j < KF->NCv; j++)
      {
        ControlPointsinitial(us, 0) =
            double((umax - umin) * i) / (KF->NCu - 1) + umin;
        ControlPointsinitial(us, 1) =
            double((vmax - vmin) * j) / (KF->NCv - 1) + vmin;
        us++;
      }
    }

    Eigen::MatrixXd ControlPoints(
        _NumberOfControlPointsU * _NumberOfControlPointsV, 2);

    std::vector<cv::KeyPoint> qe =
        Warps::Warp::getEstimates(KP1, umin, umax, vmin, vmax, KF->NCu, KF->NCv,
                                  KF->valdim, x, ControlPoints, 0, 0);
    // first order derivatives
    // dqu vector with keypoints <du/du dv/du >
    std::vector<cv::KeyPoint> dqu =
        Warps::Warp::getEstimates(KP1, umin, umax, vmin, vmax, KF->NCu, KF->NCv,
                                  KF->valdim, x, ControlPoints, 1, 0);
    // dqu vector with keypoints <dv/du dv/dv >
    std::vector<cv::KeyPoint> dqv =
        Warps::Warp::getEstimates(KP1, umin, umax, vmin, vmax, KF->NCu, KF->NCv,
                                  KF->valdim, x, ControlPoints, 0, 1);
    // second order derivatives
    std::vector<cv::KeyPoint> dquv =
        Warps::Warp::getEstimates(KP1, umin, umax, vmin, vmax, KF->NCu, KF->NCv,
                                  KF->valdim, x, ControlPoints, 1, 1);
    std::vector<cv::KeyPoint> dquu =
        Warps::Warp::getEstimates(KP1, umin, umax, vmin, vmax, KF->NCu, KF->NCv,
                                  KF->valdim, x, ControlPoints, 2, 0);
    std::vector<cv::KeyPoint> dqvv =
        Warps::Warp::getEstimates(KP1, umin, umax, vmin, vmax, KF->NCu, KF->NCv,
                                  KF->valdim, x, ControlPoints, 0, 2);
    uint counter(0);
    uint counter2(0);

    for (size_t ikp = 0; ikp < vMatchedIndices.size(); ikp++)
    {
      const size_t &idx1 = vMatchedIndices[ikp].first;
      const size_t &idx2 = vMatchedIndices[ikp].second;
      MapPoint *mapPoint = KF->GetMapPoint(idx1);
      MapPoint *mapPoint2 = KF2->GetMapPoint(idx2);
      if ((!mapPoint) or (!mapPoint2))
      {
        continue;
      }
      if ((mapPoint->isBad()) or (mapPoint2->isBad()))
      {
        continue;
      }

      auto error_i = qe[ikp].pt - KP2[ikp].pt;
      error_i.x *= KF->fx;
      error_i.y *= KF->fy;
      counter2++;

      if (cv::norm(error_i) > 10)
      {
        mapPoint2->EraseObservation(KF2);
        KF2->EraseMapPointMatch(idx2);
        continue;
      }
      /// All map points are used for the warp estimation, but only those
      /// whose reference keyframe is the estimated one are saved.
      auto pKfref = mapPoint->GetReferenceKeyFrame();
      if (pKfref != KFi)
        continue;
      mapPointsDB_[mapPoint].push_back(std::shared_ptr<DiffProp>(new DiffProp()));
      std::shared_ptr<DiffProp> diffVect = mapPointsDB_[mapPoint].back();
      diffVect->KFToKF = std::pair<KeyFrame *, KeyFrame *>(KFi, KF2i);
      diffVect->idx1 = size_t(idx1);
      diffVect->idx2 = size_t(idx2);
      // Set points
      diffVect->I1u = KP1[ikp].pt.x;
      diffVect->I1v = KP1[ikp].pt.y;
      diffVect->I2u = KP2[ikp].pt.x;
      diffVect->I2v = KP2[ikp].pt.y;
      /***********
   *  Jacobian Image Points
   *  J12 = [J12a,J12c]
   *        [J12b,J12d]
   *  J21 = [J21a,J21c]
   *        [J21b,J21d]
   ***********/
      diffVect->J12a = dqu[ikp].pt.x; // du/du
      diffVect->J12b = dqu[ikp].pt.y; // dv/du
      diffVect->J12c = dqv[ikp].pt.x; // du/dv
      diffVect->J12d = dqv[ikp].pt.y; // dv/dv

      // Inverse of the Jacobian
      diffVect->J21a = diffVect->J12d / (dqu[ikp].pt.x * dqv[ikp].pt.y -
                                         dqv[ikp].pt.x * dqu[ikp].pt.y);
      diffVect->J21b = -diffVect->J12c / (dqu[ikp].pt.x * dqv[ikp].pt.y -
                                          dqv[ikp].pt.x * dqu[ikp].pt.y);
      diffVect->J21c = -diffVect->J12b / (dqu[ikp].pt.x * dqv[ikp].pt.y -
                                          dqv[ikp].pt.x * dqu[ikp].pt.y);
      diffVect->J21d = diffVect->J12a / (dqu[ikp].pt.x * dqv[ikp].pt.y -
                                         dqv[ikp].pt.x * dqu[ikp].pt.y);

      /*************
       *  Hessians
       * H12u = [du^2/dudu du^2/dudv] = [H12uux H12uvx]
      *         [du^2/dvdu du^2/dvdv]   [H12uvx H12vvx]
      *  H12v = [dv^2/dudu dv^2/dudv] = [H12uuy H12uvy]
      *         [dv^2/dvdu dv^2/dvdv]   [H12uvy H12vvy]
      **************/
      diffVect->H12uux = dquu[ikp].pt.x;
      diffVect->H12uuy = dquu[ikp].pt.y;
      diffVect->H12uvx = dquv[ikp].pt.x;
      diffVect->H12uvy = dquv[ikp].pt.y;
      diffVect->H12vvx = dqvv[ikp].pt.x;
      diffVect->H12vvy = dqvv[ikp].pt.y;
      newInformation_[mapPoint] = true;
      counter++;
    }
    std::cout << "Points to reestimate : " << counter << " " << counter2
              << std::endl;
  }
} // namespace defSLAM
