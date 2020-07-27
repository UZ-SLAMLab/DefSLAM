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

#include "GroundTruthKeyFrame.h"
#include "GroundTruthCalculator.h"
#include "MapPoint.h"
#include <Modules/ToolsPCL/PCLNormalEstimator.h>
#include <math.h>
#include <mutex>
#include <opencv2/core/core.hpp>
#include "Modules/ToolsPCL/SmootherMLS.h"

namespace defSLAM
{
  /*****
   *  Constructor. Initialize the parameters for the groundtruth estimation.
   *  TODO: Implement with depth image as in GroundTruthFrame
   *****/

  GroundTruthKeyFrame::GroundTruthKeyFrame(Frame &F, Map *pMap,
                                           KeyFrameDatabase *pKFDB)
      : DefKeyFrame(F, pMap, pKFDB), imRight(F.imRight.clone()),
        StereoAvaliable(false)
  {
    if (F.StereoAvailable)
    {
      this->StereoAvaliable = true;
    }
  }
  /*****
   * estimateAngleErrorAndScale. This method estimates the 3D groundtruth of the
   * surface estimated for this keyframe. It use the keypoints of the left image 
   * with a normal and search for estimates the 3D in the right image. It uses the 
   * pcl library to determinate the normals of the point cloud and compares them 
   * with the estimated by the NRSfM and the SfN. NRSfM tends to be quite noisy.
   */
  float GroundTruthKeyFrame::estimateAngleErrorAndScale()
  {
    size_t nval = this->mvKeysUn.size();
    cv::Mat iLeft(this->imGray.clone());
    posMono_.clear();
    posStereo_.clear();
    std::vector<int> ptreal;
    std::vector<float> points_est;
    std::vector<float> points_gt;
    std::vector<float> normals_est;
    posMono_.reserve(nval);
    posStereo_.reserve(nval);
    std::vector<double *> covariances;
    std::vector<float> zs;
    std::vector<float> zs_est;

    for (uint i = 0; i < nval; i++)
    {
      MapPoint *pMP = this->GetMapPoint(i);
      if (!pMP)
        continue;
      if (pMP->isBad())
        continue;
      cv::Vec3f pos;

      this->surface->get3DSurfacePoint(i, pos);
      cv::Mat PosMat(4, 1, CV_32F);
      for (uint i(0); i < 3; i++)
        PosMat.at<float>(i) = pos(i);

      pos(3) = 1;

      cv::KeyPoint kp = this->mvKeysUn[i];

      std::vector<float> ps = GroundTruthTools::estimateGT(kp, iLeft, imRight, mbf, cx, cy, fx, fy);

      if (ps[2] < 0)
        continue;

      zs.push_back(ps[2]);
      zs_est.push_back(PosMat.at<float>(2));
      {
        cv::Vec3f normal;
        // For normal estimation
        if (this->surface->getNormalSurfacePoint(i, normal))
        {
          for (uint j(0); j < 3; j++)
            normals_est.push_back(normal(j));
          points_gt.push_back(ps[0]);
          points_gt.push_back(ps[1]);
          points_gt.push_back(ps[2]);
          points_est.push_back(PosMat.at<float>(0));
          points_est.push_back(PosMat.at<float>(1));
          points_est.push_back(PosMat.at<float>(2));
          std::vector<float> pm;
          pm.push_back(PosMat.at<float>(0));
          pm.push_back(PosMat.at<float>(1));
          pm.push_back(PosMat.at<float>(2));
          posMono_.push_back(pm);
          posStereo_.push_back(ps);
          covariances.push_back(pMP->covNorm);
          ptreal.push_back(i);
        }
      }
    }
    if (points_gt.size() / 3 < 20)
    {
      std::cout << "Points" << std::endl;
      return 1;
    }
    std::sort(zs.begin(), zs.end());
    std::sort(zs_est.begin(), zs_est.end());

    PCLNormalEstimator pne(points_gt, zs[zs.size() / 2] / 7);
    std::vector<float> normals_gt = pne.getNormals();
    PCLNormalEstimator pne2(points_est, zs_est[zs_est.size() / 2] / 7);
    std::vector<float> normals_sfn = pne2.getNormals();

    std::vector<float> ErrorAngleIso;
    std::vector<float> ErrorAngleSfN;

    double sumIso(0.0);
    double sumSfN(0.0);

    for (uint i(0); i < normals_gt.size() / 3; i++)
    {
      // std::cout << "Normal " << i << std::endl;
      Eigen::Vector3f nest;
      nest << normals_est[3 * i + 0], normals_est[3 * i + 1],
          normals_est[3 * i + 2];
      Eigen::Vector3f nreal;
      nreal << normals_gt[3 * i + 0], normals_gt[3 * i + 1],
          normals_gt[3 * i + 2];
      Eigen::Vector3f nsft;
      nsft << normals_sfn[3 * i + 0], normals_sfn[3 * i + 1],
          normals_sfn[3 * i + 2];
      nest.normalize();
      nreal.normalize();
      nsft.normalize();
      double aIso = std::acos(nest.transpose() * nreal) * 180 / M_PI;
      double aSfN = std::acos(nreal.transpose() * nsft) * 180 / M_PI;

      if (aIso > 90)
        aIso = 180 - aIso;

      if (aSfN > 90)
        aSfN = 180 - aSfN;

      if (aIso == aIso)
      {
        ErrorAngleIso.push_back(aIso);
        sumIso += aIso;
      }
      if (aSfN == aSfN)
      {
        ErrorAngleSfN.push_back(aSfN);
        sumSfN += aSfN;
      }
    }
    if (ErrorAngleIso.size() < 1)
      return 1;
    std::accumulate(ErrorAngleIso.begin(), ErrorAngleIso.end(), sumIso);
    std::sort(ErrorAngleIso.begin(), ErrorAngleIso.end());

    std::cout << "Mean angle Iso Error Keyframe : " << std::endl
              << "min : "
              << *std::min_element(ErrorAngleIso.begin(), ErrorAngleIso.end())
              << std::endl
              << "max : "
              << *std::max_element(ErrorAngleIso.begin(), ErrorAngleIso.end())
              << std::endl
              << "median : " << ErrorAngleIso[ErrorAngleIso.size() / 2]
              << std::endl
              << sumIso / ((double)ErrorAngleIso.size()) << " " << std::endl;

    static int is(0);
    std::ostringstream out;
    out << std::internal << std::setfill('0') << std::setw(5)
        << uint(this->mTimeStamp);
    std::string name("ErrorAngIso" + out.str() + "-" + std::to_string(is) +
                     ".txt");
    std::string name2("ErrorAngSfN" + out.str() + "-" + std::to_string(is) +
                      ".txt");
    is++;
    GroundTruthTools::saveResults(ErrorAngleIso, name);
    GroundTruthTools::saveResults(ErrorAngleSfN, name2);

    return 1;
  }
} // namespace defSLAM
