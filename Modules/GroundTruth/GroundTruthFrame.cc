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

#include "GroundTruthFrame.h"
#include "GroundTruthCalculator.h"
#include "MinMedianFilter.h"
#include "Modules/ToolsPCL/SmootherMLS.h"
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include "set_MAC.h"

#ifndef ORBSLAM
#include "DefMapPoint.h"
#endif
#include "MapPoint.h"

#include "opencv2/calib3d.hpp"
#include "opencv2/core/core.hpp"

namespace defSLAM
{

  using ORB_SLAM2::MapPoint;

  /*****
   * Constructor. See Frame constructor for more details, here 
   * With this constructor does not initialize using the stereo pair.
  */
  GroundTruthFrame::GroundTruthFrame(const cv::Mat &imGray,
                                     const double &timeStamp,
                                     ORBextractor *extractor, ORBVocabulary *voc,
                                     cv::Mat &K, cv::Mat &distCoef,
                                     const float &bf, const float &thDepth,
                                     const cv::Mat &ImRGB, const cv::Mat &ImRight,
                                     cv::Mat _mask)
      : Frame(imGray, timeStamp, extractor, voc, K, distCoef, bf, thDepth, ImRGB,
              _mask)
  {
    this->imRight = (ImRight);
    // Do not use stereo to initialize
    this->StereoAvailable = (false);
    isDepth = false;
  }
  /*****
   * Constructor. See Frame constructor for more details.
   * With this constructor does initialize using the stereo pair.
  */
  GroundTruthFrame::GroundTruthFrame(
      const cv::Mat &imLeft, const cv::Mat &imRight, const double &timeStamp,
      ORB_SLAM2::ORBextractor *extractorLeft,
      ORB_SLAM2::ORBextractor *extractorRight, ORB_SLAM2::ORBVocabulary *voc,
      cv::Mat &K, cv::Mat &distCoef, const float &bf, const float &thDepth,
      const cv::Mat &ImRGB, cv::Mat _mask)
      : Frame(imLeft, imRight, timeStamp, extractorLeft, extractorRight, voc, K,
              distCoef, bf, thDepth, ImRGB, _mask)

  {
    this->imRight = (imRight);
    // Use stereo to initialize
    this->StereoAvailable = (true);
    isDepth = false;
  }
  /*****
   * Copy constructor. See Frame constructor for more details, here 
   * With this constructor does initialize using the stereo pair.
  */
  GroundTruthFrame::GroundTruthFrame(const GroundTruthFrame &frame)
      : Frame(frame),
        isDepth(frame.isDepth)
  {
    StereoAvailable = (frame.StereoAvailable);
  }
  /*****
   * Constructor. See Frame constructor for more details, here 
   * we set up the parameters to estimate the ground truth. 
   * With this constructor does initialize using the stereo pair.
   * And can use depth image. Usually with CT depth image
  */
  GroundTruthFrame::GroundTruthFrame(const cv::Mat &imGray,
                                     const double &timeStamp,
                                     ORBextractor *extractor, ORBVocabulary *voc,
                                     cv::Mat &K, cv::Mat &distCoef,
                                     const float &bf, const float &thDepth,
                                     const cv::Mat &ImRGB, const cv::Mat &imDepth,
                                     bool isDepth, cv::Mat _mask)
      : Frame(imGray, timeStamp, extractor, voc, K, distCoef, bf, thDepth, ImRGB,
              _mask),
        isDepth(isDepth),
        imDepth(imDepth.clone())
  {
    // Not to use stereo to initialize
    this->StereoAvailable = (false);
  }
  /*****
   * Estimate3DScale. This method estimates the 3D groundtruth of the
   * map points. It reprojects them into the left image and search for
   * estimates the 3D. It can use a depth image or use the right image 
   * and search by correlation. It removes the outliers and once done it
   * estimates the current scale with a minimal median procedure. 
  */
  double GroundTruthFrame::Estimate3DScale(Map *map)
  {
    const vector<MapPoint *> &vpMPs = map->GetAllMapPoints();

    set<MapPoint *> spRefMPs(vpMPs.begin(), vpMPs.end());

    posMono_.clear();
    posStereo_.clear();

    cv::Mat left_disp;
    std::vector<std::vector<float>> posMonoInit_;
    std::vector<std::vector<float>> posStereoInit_;
    posMono_.reserve(vpMPs.size());
    posStereo_.reserve(vpMPs.size());
    std::vector<float> zs;
    for (vector<MapPoint *>::const_iterator sit = vpMPs.begin(),
                                            send = vpMPs.end();
         sit != send; sit++)
    {
#ifndef ORBSLAM
      if (!static_cast<DefMapPoint *>(*sit)->getFacet())
      {
        continue;
      }
#endif
      if ((*sit)->isBad())
        continue;
      // OpenCV_Template Matching tutorial :
      // https://docs.opencv.org/master/de/da9/tutorial_template_matching.html
      cv::Mat pos = (*sit)->GetWorldPos();
      cv::KeyPoint kp = this->ProjectPoints(pos);
      if (kp.pt.x < 0)
        continue;
      if (isDepth)
      {
        float dpth = imDepth.at<float>(kp.pt.y, kp.pt.x);
        const cv::Mat Pc = mRcw * pos + mtcw;
        std::vector<float> pm;
        pm.push_back(Pc.at<float>(0));
        pm.push_back(Pc.at<float>(1));
        pm.push_back(Pc.at<float>(2));
        std::vector<float> ps;
        ps.push_back(dpth * (((float)kp.pt.x - cx) / fx));
        ps.push_back(dpth * (((float)kp.pt.y - cy) / fy));
        ps.push_back(dpth);
        std::unique_lock<std::mutex> lck(mutexPoints);
        posMono_.push_back(pm);
        posStereo_.push_back(ps);
        continue;
      }

      std::vector<float> ps = GroundTruthTools::estimateGT(kp, ImGray, imRight, mbf, cx, cy, fx, fy);

      if (ps[2] < 0)
        continue;

      const cv::Mat Pc = mRcw * pos + mtcw;
      std::vector<float> pm;
      pm.reserve(3);
      pm.push_back(Pc.at<float>(0));
      pm.push_back(Pc.at<float>(1));
      pm.push_back(Pc.at<float>(2));
      zs.push_back(ps[2]);
      std::unique_lock<std::mutex> lck(mutexPoints);
      posMonoInit_.push_back(pm);
      posStereoInit_.push_back(ps);
    }
    std::cout << "POINTS EVALUATED : " << posMono_.size() << std::endl;
    if (!isDepth)
    {
      if (posMonoInit_.size() < 20)
        return 1;
      double filtering_time = ((double)cv::getTickCount());
      std::sort(zs.begin(), zs.end());

      SmootherMLS smls(1, zs[zs.size() / 2] / 7);
      std::vector<int> notOutliers = smls.outlierRemovalRadius(posStereoInit_);
      posMono_.reserve(posStereoInit_.size());
      posStereo_.reserve(posStereoInit_.size());
      if (notOutliers.size() < 20)
        return 1;
      for (uint i(0); i < notOutliers.size(); i++)
      {
        posMono_.push_back(posMonoInit_[notOutliers[i]]);
        posStereo_.push_back(posStereoInit_[notOutliers[i]]);
      }
      if (posMono_.size() < 20)
        return 1;
      filtering_time =
          ((double)cv::getTickCount() - filtering_time) / cv::getTickFrequency();
      std::cout << "Time 2 " << filtering_time << "ms " << std::endl;
      std::cout << "POINTS EVALUATED : " << posMono_.size() << "/" << notOutliers.size() << "/"
                << vpMPs.size() << std::endl;
    }
    if (posMono_.size() < 50)
    {
      posMono_.clear();
      posStereo_.clear();
      return 1;
    }
    double Scale = GroundTruthTools::scaleMinMedian(posMono_, posStereo_);

    std::cout << "Scale Corr stereo : " << Scale << std::endl;
    return Scale;
  }
  /*****
   * Estimate3DError. Estimates the euclidean distance between the ground truth
   * and the scaled monocular map point. It also saves the results into a txt 
   * file
  */
  double GroundTruthFrame::Estimate3DError(Map *map, const double &s)
  {
    std::vector<float> Error;
    Error.reserve(posMono_.size());
    int count(0);
    double acc(0.0);
    mvLocalMapPoints.clear();
    mvStereoMapPoints.clear();
    std::unique_lock<std::mutex> lck(mutexPoints);

    for (size_t i(0); i < posMono_.size(); i++)
    {
      double er = sqrt(pow(posStereo_[i][0] - s * posMono_[i][0], 2) +
                       pow(posStereo_[i][1] - s * posMono_[i][1], 2) +
                       pow(posStereo_[i][2] - s * posMono_[i][2], 2));
      Error.push_back(er);
      acc += er;
      count++;
    }

    double sum(0.0);
    std::accumulate(Error.begin(), Error.end(), sum);
    double invc = 1 / ((double)count);
    std::cout << "Mean Error Surf : " << acc * invc << " " << Error.size() << " "
              << posMono_.size() << std::endl;
    std::ostringstream out;
    out << std::internal << std::setfill('0') << std::setw(5)
        << uint(this->mTimeStamp);
    std::string name("ErrorGTs" + out.str() + ".txt");
    GroundTruthTools::saveResults(Error, name);
    return acc * invc;
  }

  // return monocular map point 3D pose
  std::vector<std::vector<float>> GroundTruthFrame::getPosMono()
  {
    std::unique_lock<std::mutex> lck(mutexPoints);
    return posMono_;
  }

  // return groundtruth map point 3D pose
  std::vector<std::vector<float>> GroundTruthFrame::getPosStereo()
  {
    std::unique_lock<std::mutex> lck(mutexPoints);
    return posStereo_;
  }
} // namespace defSLAM
