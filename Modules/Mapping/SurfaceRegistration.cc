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

#include "SurfaceRegistration.h"
#include "Converter.h"
#include "DefKeyFrame.h"
#include "DefMapPoint.h"
#include "DefOptimizer.h"
#include "Surface.h"
#include "SurfacePoint.h"
#include "GroundTruthCalculator.h"

namespace defSLAM
{
  using ORB_SLAM2::Converter;
  using ORB_SLAM2::MapPoint;

  SurfaceRegistration::SurfaceRegistration(KeyFrame *Keyframe, double chiLimit_,
                                           bool check_chi)
      : refKF(Keyframe), chiLimit_(chiLimit_), check_chi(check_chi) {}

  SurfaceRegistration::~SurfaceRegistration() { refKF = nullptr; }

  /**************
     * Function to align the point clouds. It returns true if the 
     * registration is correct. It use a Levenger-Marquard to 
     * recover the Sim(3) that align the surface estimated with the
     * point clouds register in the keyframe. 
     * ************/
  bool SurfaceRegistration::registerSurfaces()
  {
    std::vector<std::vector<float>> cloud1pc;
    std::vector<std::vector<float>> cloud2pc;
    cv::Mat Twc = refKF->GetPoseInverse();
    cv::Mat TcwKf = refKF->GetPose();

    cloud1pc.reserve(refKF->mvKeys.size());
    cloud2pc.reserve(refKF->mvKeys.size());

    uint CounterMP(0);
    // Set up the 3D points cloud to make the sim(3) alignment
    for (uint i(0); i < refKF->mvKeys.size(); i++)
    {
      MapPoint *pMP = refKF->GetMapPoint(i);
      if (pMP)
      {
        if (pMP->isBad())
          continue;

        if (static_cast<DefMapPoint *>(pMP)->getFacet())
        {
          if (!pMP->covNorm)
            continue;
          cv::Mat x3DMap(4, 1, CV_32F);
          cv::Mat x3Dy = static_cast<DefMapPoint *>(pMP)
                             ->PosesKeyframes[refKF]
                             .clone();
          if (x3Dy.empty())
            continue;
          std::vector<float> pt1;
          pt1.push_back(x3Dy.at<float>(0, 0));
          pt1.push_back(x3Dy.at<float>(1, 0));
          pt1.push_back(x3Dy.at<float>(2, 0));
          x3DMap.at<float>(0, 0) = x3Dy.at<float>(0, 0);
          x3DMap.at<float>(1, 0) = x3Dy.at<float>(1, 0);
          x3DMap.at<float>(2, 0) = x3Dy.at<float>(2, 0);
          x3DMap.at<float>(3, 0) = 1;
          cloud1pc.push_back(std::move(pt1));
          cv::Vec3f x3DsurfKf;
          static_cast<DefKeyFrame *>(refKF)->surface->get3DSurfacePoint(
              i, x3DsurfKf);
          cv::Mat x3DSurW(4, 1, CV_32F);
          x3DSurW.at<float>(0, 0) = x3DsurfKf(0);
          x3DSurW.at<float>(1, 0) = x3DsurfKf(1);
          x3DSurW.at<float>(2, 0) = x3DsurfKf(2);
          x3DSurW.at<float>(3, 0) = 1;
          cv::Mat x3DKf = (Twc)*x3DSurW;
          std::vector<float> pt2;
          pt2.push_back(x3DKf.at<float>(0, 0));
          pt2.push_back(x3DKf.at<float>(1, 0));
          pt2.push_back(x3DKf.at<float>(2, 0));
          cloud2pc.push_back(std::move(pt2));
          CounterMP++;
        }
      }
    }

    if (CounterMP < 15)
      return false;
    // Used to initialize the scale. Although is from GroundTruthTools
    // is not using anything of groundtruth to the scale computation.
    float scale = GroundTruthTools::scaleMinMedian(cloud2pc, cloud1pc);
    Eigen::Matrix4f TwcEigen;
    cv::cv2eigen(Twc, TwcEigen);

    Eigen::Matrix4f transform_;
    transform_ << 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1;

    Eigen::Matrix4d transform_d;
    for (uint i(0); i < 4; i++)
    {
      for (uint j(0); j < 4; j++)
      {
        transform_d(i, j) = transform_(i, j);
      }
    }

    g2o::Sim3 transf(transform_d.block(0, 0, 3, 3), transform_d.block(0, 3, 3, 1),
                     scale);

    double Huber(0.01);
    double chi(pow(chiLimit_, 2));
    bool aceptable =
        Optimizer::OptimizeHorn(cloud2pc, cloud1pc, transf, chi, Huber);

    if ((!aceptable) && (check_chi))
      return false;

    cv::Mat mScw = Converter::toCvMat(transf);
    Eigen::MatrixXf mScwEigen;
    cv::cv2eigen(mScw, mScwEigen);
    TwcEigen = mScwEigen * TwcEigen;

    Eigen::MatrixXf TT2 =
        TwcEigen.block(0, 0, 3, 3) * (TwcEigen.block(0, 0, 3, 3)).transpose();
    double s22 = std::sqrt(TT2(0, 0)); // Recover the scale
    static_cast<DefKeyFrame *>(refKF)->surface->applyScale(s22);

    TwcEigen.block(0, 0, 3, 3) = TwcEigen.block(0, 0, 3, 3) / (s22);
    Eigen::MatrixXf Tcw = TwcEigen.inverse();
    cv::eigen2cv(Tcw, mScw);
    refKF->SetPose(mScw); // Change the SE(3) pose of the camera.

    return true;
  }
} // namespace defSLAM
