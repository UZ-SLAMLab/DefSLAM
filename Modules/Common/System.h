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

#ifndef SYSTEM_H
#define SYSTEM_H

#include <opencv2/core/core.hpp>
#include <string>
#include <thread>

#include "FrameDrawer.h"
#include "KeyFrameDatabase.h"
#include "LocalMapping.h"
#include "LoopClosing.h"
#include "Map.h"
#include "MapDrawer.h"
#include "ORBVocabulary.h"
#include "Tracking.h"
#include "Viewer.h"

namespace ORB_SLAM2
{
  class FrameDrawer;
  class Map;
  class Tracking;
  class LocalMapping;
  class LoopClosing;
  class Viewer;
} // namespace ORB_SLAM2

namespace defSLAM
{
  using ORB_SLAM2::FrameDrawer;
  using ORB_SLAM2::KeyFrameDatabase;
  using ORB_SLAM2::LocalMapping;
  using ORB_SLAM2::LoopClosing;
  using ORB_SLAM2::Map;
  using ORB_SLAM2::MapDrawer;
  using ORB_SLAM2::ORBextractor;
  using ORB_SLAM2::ORBVocabulary;
  using ORB_SLAM2::Tracking;
  using ORB_SLAM2::Viewer;

  class System
  {
  public:
    // Input sensor
    enum eSensor
    {
      MONOCULAR = 0,
      STEREO = 1,
      RGBD = 2
    };

  public:
    System() = default;

    // Initialize the SLAM system. It launches the Local Mapping, Loop Closing and
    // Viewer threads.
    System(const string &strVocFile, const string &strSettingsFile, const bool bUseViewer = true);

    // Proccess the given stereo frame. Images must be synchronized and rectified.
    // Input images: RGB (CV_8UC3) or grayscale (CV_8U). RGB is converted to
    // grayscale.
    // Returns the camera pose (empty if tracking fails).
    cv::Mat TrackStereo(const cv::Mat &imLeft, const cv::Mat &imRight,
                        const double &timestamp);

    // Process the given rgbd frame. Depthmap must be registered to the RGB frame.
    // Input image: RGB (CV_8UC3) or grayscale (CV_8U). RGB is converted to
    // grayscale.
    // Input depthmap: Float (CV_32F).
    // Returns the camera pose (empty if tracking fails).
    cv::Mat TrackRGBD(const cv::Mat &im, const cv::Mat &depthmap,
                      const double &timestamp);

    // Proccess the given monocular frame
    // Input images: RGB (CV_8UC3) or grayscale (CV_8U). RGB is converted to
    // grayscale.
    // Returns the camera pose (empty if tracking fails).
    virtual cv::Mat TrackMonocular(const cv::Mat &im, const double &timestamp,
                                   const cv::Mat _mask = cv::Mat());
    // virtual cv::Mat TrackMonocularwithOut(const cv::Mat &im,const cv::Mat
    // &imOut, const double &timestamp,const cv::Mat _mask = cv::Mat());

    // Using stereo for calculate the ground truth
    virtual cv::Mat TrackMonocularGT(const cv::Mat &im, const cv::Mat &imRight,
                                     const double &timestamp,
                                     const cv::Mat _mask = cv::Mat());

    // Using CT for depth image the ground truth (It may work with depth image too)
    cv::Mat TrackMonocularCTGT(const cv::Mat &im, const cv::Mat &CTdepth,
                               const double &timestamp, const cv::Mat _mask);

    // This stops local mapping thread (map building) and performs only camera
    // tracking.
    void ActivateLocalizationMode();

    // This resumes local mapping thread and performs SLAM again.
    void DeactivateLocalizationMode();

    // Returns true if there have been a big map change (loop closure, global BA)
    // since last call to this function
    bool MapChanged();

    // Reset the system (clear map)
    void Reset();

    // Restart with a different thickening or propagation band size
    void Restart(uint localzone, uint propagationzone);

    // All threads will be requested to finish.
    // It waits until all threads have finished.
    // This function must be called before saving the trajectory.
    virtual void Shutdown();

    // Information from most recent processed frame
    // You can call this right after TrackMonocular (or stereo or RGBD)
    int GetTrackingState();

  protected:
    // Input sensor
    eSensor mSensor;

    // ORB vocabulary used for place recognition and feature matching.
    ORBVocabulary *mpVocabulary;

    // KeyFrame database for place recognition (relocalization and loop
    // detection).
    KeyFrameDatabase *mpKeyFrameDatabase;

    // Map structure that stores the pointers to all KeyFrames and MapPoints.
    Map *mpMap;

    // Tracker. It receives a frame and computes the associated camera pose.
    // It also decides when to insert a new keyframe, create some new MapPoints
    // and
    // performs relocalization if tracking fails. // Non-Rigid Tracker. Class
    // inherited from Tracking that uses a laplacian
    // and inextensible method to ,once received a frame, compute the camera pose
    // and
    // the warp of the map.
    Tracking *mpTracker;

    // Local Mapper. It manages the local map and performs local bundle
    // adjustment.
    LocalMapping *mpLocalMapper;

    // Loop Closer. It searches loops with every new keyframe. If there is a loop
    // it performs
    // a pose graph optimization and full bundle adjustment (in a new thread)
    // afterwards.
    LoopClosing *mpLoopCloser;

    // The viewer draws the map and the current camera pose. It uses Pangolin.
    Viewer *mpViewer;

    FrameDrawer *mpFrameDrawer;
    MapDrawer *mpMapDrawer;

    // System threads: Local Mapping, Loop Closing, Viewer.
    // The Tracking thread "lives" in the main execution thread that creates the
    // System object.
    std::thread *mptLocalMapping;
    std::thread *mptLoopClosing;
    std::thread *mptViewer;

    // Reset flag
    std::mutex mMutexReset;
    bool mbReset;

    // Change mode flags
    std::mutex mMutexMode;
    bool mbActivateLocalizationMode;
    bool mbDeactivateLocalizationMode;

    // Tracking state
    int mTrackingState;
    std::vector<MapPoint *> mTrackedMapPoints;
    std::vector<cv::KeyPoint> mTrackedKeyPointsUn;

    std::mutex mMutexState;
    std::mutex mMutexdata;
  };

} // namespace defSLAM

#endif // SYSTEM_H
