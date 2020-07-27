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

#ifndef DEFTRACKING_H
#define DEFTRACKING_H

#include <opencv2/core/core.hpp>
#include <opencv2/features2d/features2d.hpp>

#include <Tracking.h>
#include <mutex>
#include <opencv2/core/eigen.hpp>

namespace ORB_SLAM2
{
  class Tracking;
  class ORBmatcher;
} // namespace ORB_SLAM2

namespace defSLAM
{
  using ORB_SLAM2::FrameDrawer;
  using ORB_SLAM2::KeyFrameDatabase;
  using ORB_SLAM2::Map;
  using ORB_SLAM2::MapDrawer;
  using ORB_SLAM2::ORBVocabulary;
  using ORB_SLAM2::System;
  using ORB_SLAM2::Tracking;

  class DefKeyFrame;
  class DefTracking : public Tracking
  {

  public:
    DefTracking(System *pSys, ORBVocabulary *pVoc,
                FrameDrawer *pFrameDrawer, MapDrawer *pMapDrawer,
                Map *pMap, KeyFrameDatabase *pKFDB,
                const string &strSettingPath,
                const int sensor = ORB_SLAM2::System::MONOCULAR,
                bool viewerOn = false);

  public:
    //Main function of tracking where the map is considered deformable.
    virtual bool TrackLocalMap();

    // Initial rigid tracking to get more matches and remove some outliers.
    virtual bool TrackWithMotionModel();

    // Run for image with ground truth through a stereo pair
    virtual cv::Mat GrabImageMonocularGT(const cv::Mat &im,
                                         const cv::Mat &imRight,
                                         const double &timestamp,
                                         cv::Mat _mask = cv::Mat());

    // Run for image with ground truth through a depth image. It is used for
    // the CT image of the phantom datasets, but it should be usable for rgbd image.
    virtual cv::Mat GrabImageMonocularCTGT(const cv::Mat &im,
                                           const cv::Mat &imDepth,
                                           const double &timestamp,
                                           cv::Mat _mask = cv::Mat());

  protected:
    // Main function of tracking.
    void Track() override;

    // Remove matches of the current keyframe.
    virtual void CleanMatches();

    // Not useful
    void EraseTemporalPoints();

    // Create new keyframe for deformable SLAM.
    void CreateNewKeyFrame() override;

    // Initialize scene with a plane.
    void MonocularInitialization() override;

    // Check local points with the covisible keyframes.
    // Check there are no repeated points.
    void UpdateLocalPoints() override;

    // Update data to last frame.
    void UpdateLastFrame() override;

  protected:
    KeyFrame *keyframe; // main keyframe
    uint LocalZone;
    ofstream myfile;
    bool saveResults;
  };

} // namespace defSLAM

#endif // DEFORMABLE TRACKING_H
