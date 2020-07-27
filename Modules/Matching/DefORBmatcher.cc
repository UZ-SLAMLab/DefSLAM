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

#include "DefORBmatcher.h"
#include <limits.h>
#include <opencv2/core/core.hpp>
#include <opencv2/features2d/features2d.hpp>

#include "DefKeyFrame.h"
#include "DefMapPoint.h"
#include "Frame.h"
#include "KeyFrame.h"
#include "MapPoint.h"
#include "Schwarp.h"
#include "Thirdparty/BBS/bbs_coloc.h"
#include <Eigen/Core>
#include <opencv2/core/core.hpp>
#include <opencv2/features2d/features2d.hpp>

using namespace std;

namespace defSLAM
{
  // Constructor.
  DefORBmatcher::DefORBmatcher(float nnratio, bool checkOri)
      : ORBmatcher(nnratio, checkOri) {}

  // Find more matches using a Schwarp. x is the array with the nodes of the warp.
  void DefORBmatcher::findbyWarp(
      KeyFrame *Kf1, KeyFrame *Kf2, vector<pair<size_t, size_t>> &vMatchedIndices,
      double (&x)[_NumberOfControlPointsU * _NumberOfControlPointsV * 2],
      double lambda)
  {

    this->CalculateInitialSchwarp(Kf1, Kf2, vMatchedIndices, x, lambda);

    std::vector<pair<size_t, size_t>> vMatchedIndices2;

    this->searchBySchwarp(Kf1, Kf2, x, vMatchedIndices2);
    std::cout << "New points : " << vMatchedIndices2.size() << std::endl;
    for (uint i(0); i < vMatchedIndices2.size(); i++)
    {
      MapPoint *pMP = Kf1->GetMapPoint(vMatchedIndices2[i].first);
      if (pMP)
      {
        pMP->AddObservation(Kf2, vMatchedIndices2[i].second);
        Kf2->addMapPoint(pMP, vMatchedIndices2[i].second);
      }
    }

    vMatchedIndices.insert(vMatchedIndices.end(), vMatchedIndices2.begin(),
                           vMatchedIndices2.end());
  }

  /*********************
  * Initialize a warp and make a refinement with the Schwarzian equation
  * to search for map points. HERE IS INITIALIZED THE WARP. IF YOU DO NOT
  * CALL findbyWarp IN SchwarpDatabase.cc, YOU MUST INITIALIZE THE WARP
  * BEFORE OPTIMIZING THE SCHWARP. (To initialize a warp, use the function
  * Warps::Warp::initialize() from Mapping/Schwarp.h)
  **********************/
  void DefORBmatcher::CalculateInitialSchwarp(
      KeyFrame *Kf1i, KeyFrame *Kf2i,
      vector<pair<size_t, size_t>> &vMatchedIndices,
      double (&x)[_NumberOfControlPointsU * _NumberOfControlPointsV * 2],
      double lambda)
  {
    // Calculate From One To Two
    DefKeyFrame *KF = static_cast<DefKeyFrame *>(Kf1i);
    DefKeyFrame *KF2 = static_cast<DefKeyFrame *>(Kf2i);

    std::vector<cv::KeyPoint> KP1;
    std::vector<cv::KeyPoint> KP2;
    std::vector<float> invSigmas;

    for (size_t ikp = 0; ikp < vMatchedIndices.size(); ikp++)
    {
      const int idx1 = vMatchedIndices[ikp].first;
      const int idx2 = vMatchedIndices[ikp].second;
      const cv::KeyPoint kp1 = KF->mpKeypointNorm[idx1];
      const cv::KeyPoint kp2 = KF2->mpKeypointNorm[idx2];
      KP1.push_back(kp1);
      KP2.push_back(kp2);
      const cv::KeyPoint &kpUn = KF->mvKeysUn[idx1];
      const float &invSigma2 = KF->mvInvLevelSigma2[kpUn.octave];
      invSigmas.push_back(sqrt(invSigma2));
    }

    for (uint i(0); i < _NumberOfControlPointsU * _NumberOfControlPointsU * 2;
         i++)
    {
      x[i] = 0.0;
    }

    // This must be called before optimizing.
    Warps::Warp::initialize(KP1, KP2, lambda, KF->umin, KF->umax, KF->vmin, KF->vmax,
                            KF->NCu, KF->NCv, KF->valdim, x);

    for (uint i(0); i < _NumberOfControlPointsU * _NumberOfControlPointsU * 2;
         i++)
    {
      if (std::isnan(x[i]))
        x[i] = 0.0;
    }

    ceres::CostFunction *Rep =
        new Warps::Warp(KP1, KP2, invSigmas, KF->umin, KF->umax, KF->vmin,
                        KF->vmax, KF->NCu, KF->NCv, KF->valdim, KF->fx, KF->fy);

    ceres::Problem problem;
    problem.AddParameterBlock(x, KF->NCu * KF->NCv * 2);
    problem.AddResidualBlock(Rep, new ceres::HuberLoss(5.77), x);

    ceres::Problem::EvaluateOptions optionsEv;

    double total_cost = 0.0;
    vector<double> residuals;
    problem.Evaluate(optionsEv, &total_cost, &residuals, nullptr, nullptr);
    std::vector<size_t> outliers;
    // Remove points that bad predicted by the Warp.
    for (size_t i(0); i < vMatchedIndices.size(); i++)
    {
      double error = residuals[2 * i] * residuals[2 * i] +
                     residuals[2 * i + 1] * residuals[2 * i + 1];
      if (error > 20)
      {
        outliers.push_back(i);
      }
    }

    for (auto i : outliers)
    {
      const auto idx2 = vMatchedIndices[i].second;
      KF2->EraseMapPointMatch(idx2);
      vMatchedIndices.erase(vMatchedIndices.begin() + i);
    }
  }

  // It use the result of the schwarzian warp to search for map points.
  int DefORBmatcher::searchBySchwarp(
      KeyFrame *pKF1, KeyFrame *pKF2,
      double (&x)[_NumberOfControlPointsU * _NumberOfControlPointsV * 2],
      std::vector<pair<size_t, size_t>> &vMatchedPairs)
  {
    DefKeyFrame *dKF = static_cast<DefKeyFrame *>(pKF1);
    DefKeyFrame *dKF2 = static_cast<DefKeyFrame *>(pKF2);

    std::vector<int> listMapPoints;
    std::vector<cv::KeyPoint> lskeypoints;

    for (uint i(0); i < dKF->mpKeypointNorm.size(); i++)
    {
      auto pMP = dKF->GetMapPoint(i);
      if (!pMP)
        continue;
      if (pMP->isBad())
        continue;
      if (pMP->IsInKeyFrame(dKF2))
        continue;
      listMapPoints.push_back(i);
      lskeypoints.push_back(dKF->mpKeypointNorm[i]);
    }

    if (lskeypoints.size() < 1)
      return 0;

    // Find matches between not tracked keypoints
    int nmatches = 0;
    vector<int> vMatches12(listMapPoints.size(), -1);

    Eigen::Matrix<double, _NumberOfControlPointsU * _NumberOfControlPointsV, 2>
        ControlPointsinitial;

    uint us(0);
    for (int i(0); i < dKF->NCu; i++)
    {
      for (int j(0); j < dKF->NCv; j++)
      {
        ControlPointsinitial(us, 0) =
            double((dKF->umax - dKF->umin) * i) / (dKF->NCu - 1) + dKF->umin;
        ControlPointsinitial(us, 1) =
            double((dKF->vmax - dKF->vmin) * j) / (dKF->NCv - 1) + dKF->vmin;
        us++;
      }
    }

    Eigen::MatrixXd ControlPoints(
        _NumberOfControlPointsU * _NumberOfControlPointsV, 2);

    std::vector<cv::KeyPoint> KP2_e = Warps::Warp::getEstimates(
        lskeypoints, dKF->umin, dKF->umax, dKF->vmin, dKF->vmax, dKF->NCu,
        dKF->NCv, dKF->valdim, x, ControlPoints, 0, 0);

    // Iterations over the keypoints with map point not assigned in kf1
    for (uint i(0); i < lskeypoints.size(); i++)
    {
      float x = KP2_e[i].pt.x * dKF2->fx + dKF2->cx;
      float y = KP2_e[i].pt.y * dKF2->fy + dKF2->cy;

      if (!dKF2->IsInImage(x, y))
        continue;

      const cv::Mat &d1 = dKF->mDescriptors.row(listMapPoints[i]);
      int bestDist = TH_LOW;
      int bestIdx2 = -1;
      float th = 2;
      const auto &features = dKF2->GetFeaturesInArea(x, y, th);
      // Iterations over the keypoints without map point assigned in kf2
      for (uint j(0); j < features.size(); j++)
      {
        auto pMP = dKF2->GetMapPoint(features[j]);
        if (pMP)
          continue;

        const cv::Mat &d2 = dKF2->mDescriptors.row(features[j]);
        const int dist = DescriptorDistance(d1, d2);

        if (dist < TH_LOW && dist < bestDist)
        {
          bestIdx2 = features[j];
          bestDist = dist;
        }
      }

      if (bestIdx2 >= 0)
      {
        vMatches12[i] = bestIdx2;
        nmatches++;
      }
    }
    vMatchedPairs.clear();
    vMatchedPairs.reserve(nmatches);

    for (size_t i = 0, iend = vMatches12.size(); i < iend; i++)
    {
      if (vMatches12[i] < 0)
        continue;
      vMatchedPairs.push_back(make_pair(listMapPoints[i], vMatches12[i]));
    }

    return nmatches;
  }

  // Copy of function SearchByProjection of ORBmatcher.h. It is used in
  // in the tracking. It only get points embedded in the 3D template
  int DefORBmatcher::SearchByProjection(Frame &CurrentFrame,
                                        const Frame &LastFrame,
                                        const float th,
                                        const bool bMono)
  {
    int nmatches = 0;

    // Rotation Histogram (to check rotation consistency)
    vector<int> rotHist[HISTO_LENGTH];
    for (int i = 0; i < HISTO_LENGTH; i++)
      rotHist[i].reserve(500);
    const float factor = 1.0f / HISTO_LENGTH;

    const cv::Mat Rcw = CurrentFrame.mTcw.rowRange(0, 3).colRange(0, 3);
    const cv::Mat tcw = CurrentFrame.mTcw.rowRange(0, 3).col(3);

    const cv::Mat twc = -Rcw.t() * tcw;

    const cv::Mat Rlw = LastFrame.mTcw.rowRange(0, 3).colRange(0, 3);
    const cv::Mat tlw = LastFrame.mTcw.rowRange(0, 3).col(3);

    const cv::Mat tlc = Rlw * twc + tlw;

    const bool bForward = tlc.at<float>(2) > CurrentFrame.mb && !bMono;
    const bool bBackward = -tlc.at<float>(2) > CurrentFrame.mb && !bMono;
    for (int i = 0; i < LastFrame.N; i++)
    {
      MapPoint *pMP = LastFrame.mvpMapPoints[i];

      if (pMP)
      {
        if (!LastFrame.mvbOutlier[i])
        {
          if (pMP->isBad())
            continue;
          if (!static_cast<DefMapPoint *>(pMP)->getFacet())
            continue;
          // Project
          cv::Mat x3Dw = pMP->GetWorldPos();
          cv::Mat x3Dc = Rcw * x3Dw + tcw;

          const float xc = x3Dc.at<float>(0);
          const float yc = x3Dc.at<float>(1);
          const float invzc = 1.0 / x3Dc.at<float>(2);

          if (invzc < 0)
            continue;

          float u = CurrentFrame.fx * xc * invzc + CurrentFrame.cx;
          float v = CurrentFrame.fy * yc * invzc + CurrentFrame.cy;

          if (u < CurrentFrame.mnMinX || u > CurrentFrame.mnMaxX)
            continue;
          if (v < CurrentFrame.mnMinY || v > CurrentFrame.mnMaxY)
            continue;

          int nLastOctave = LastFrame.mvKeys[i].octave;

          // Search in a window. Size depends on scale
          float radius = th * CurrentFrame.mvScaleFactors[nLastOctave];

          vector<size_t> vIndices2;

          if (bForward)
            vIndices2 = CurrentFrame.GetFeaturesInArea(u, v, radius, nLastOctave);
          else if (bBackward)
            vIndices2 =
                CurrentFrame.GetFeaturesInArea(u, v, radius, 0, nLastOctave);
          else
            vIndices2 = CurrentFrame.GetFeaturesInArea(
                u, v, radius, nLastOctave - 1, nLastOctave + 1);

          if (vIndices2.empty())
            continue;

          const cv::Mat dMP = pMP->GetDescriptor();

          int bestDist = 256;
          int bestIdx2 = -1;

          for (vector<size_t>::const_iterator vit = vIndices2.begin(),
                                              vend = vIndices2.end();
               vit != vend; vit++)
          {
            const size_t i2 = *vit;
            if (CurrentFrame.mvpMapPoints[i2])
              if (CurrentFrame.mvpMapPoints[i2]->Observations() > 0)
                continue;

            if (CurrentFrame.mvuRight[i2] > 0)
            {
              const float ur = u - CurrentFrame.mbf * invzc;
              const float er = fabs(ur - CurrentFrame.mvuRight[i2]);
              if (er > radius)
                continue;
            }

            const cv::Mat &d = CurrentFrame.mDescriptors.row(i2);

            const int dist = DescriptorDistance(dMP, d);

            if (dist < bestDist)
            {
              bestDist = dist;
              bestIdx2 = i2;
            }
          }

          if (bestDist <= TH_HIGH)
          {
            CurrentFrame.mvpMapPoints[bestIdx2] = pMP;
            nmatches++;

            if (mbCheckOrientation)
            {
              float rot = LastFrame.mvKeysUn[i].angle -
                          CurrentFrame.mvKeysUn[bestIdx2].angle;
              if (rot < 0.0)
                rot += 360.0f;
              int bin = round(rot * factor);
              if (bin == HISTO_LENGTH)
                bin = 0;
              assert(bin >= 0 && bin < HISTO_LENGTH);
              rotHist[bin].push_back(bestIdx2);
            }
          }
        }
      }
    }

    // Apply rotation consistency
    if (mbCheckOrientation)
    {
      int ind1 = -1;
      int ind2 = -1;
      int ind3 = -1;

      ComputeThreeMaxima(rotHist, HISTO_LENGTH, ind1, ind2, ind3);

      for (int i = 0; i < HISTO_LENGTH; i++)
      {
        if (i != ind1 && i != ind2 && i != ind3)
        {
          for (size_t j = 0, jend = rotHist[i].size(); j < jend; j++)
          {
            CurrentFrame.mvpMapPoints[rotHist[i][j]] =
                static_cast<MapPoint *>(NULL);
            nmatches--;
          }
        }
      }
    }

    return nmatches;
  }
} // namespace defSLAM
