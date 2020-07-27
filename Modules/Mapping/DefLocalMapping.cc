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

#include "DefLocalMapping.h"
#include "DefKeyFrame.h"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "GroundTruthKeyFrame.h"
#include "NormalEstimator.h"
#include "PolySolver.h"
#include "SchwarpDatabase.h"
#include "ShapeFromNormals.h"
#include "SurfaceRegistration.h"
#include <Converter.h>
#include <chrono>
#include <ctime>
#include <numeric>
#include <stdio.h>
#include <unistd.h>

namespace defSLAM
{
  /***********************************
   * Constructor of DefLocalMapping. It calls to the constructor to local mapping
   * and initializes the parameters of the deformable local mapping. we use
   * strSettingPath to initialize    *  
   *  pointsToTemplate_: keypoints in the unexplored area to create a new template
   *  chiLimit_ : Limit in residual error of the Surface registration to accept that
   *            a template has been correctly aligned.
   *  reg: Weight for the Schwarzian regularizer in Schwarp estimation.
   *  bendingReg_ : Weight for the beding regularizer in Shape-from-normals.
   *********************/
  DefLocalMapping::DefLocalMapping(Map *pMap,
                                   const string &strSettingPath)
      : LocalMapping(pMap, nullptr, 0.0),
        createTemplate_(false),
        pointsToTemplate_(100),
        chiLimit_(0.07)
  {
    cv::FileStorage fSettings(strSettingPath, cv::FileStorage::READ);
    pointsToTemplate_ = fSettings["LocalMapping.pointsToTemplate"];
    chiLimit_ = fSettings["LocalMapping.chiLimit"];
    double reg_ = fSettings["LocalMapping.Schwarp.Regularizer"];
    bendingReg_ = fSettings["LocalMapping.Bending"];
    saveResults_ = bool(int(fSettings["Viewer.SaveResults"]));
    warpDB_ = new SchwarpDatabase(reg_);
  }

  /***********************************
   * Destructor of DefLocalMapping. Just to remove the SchwarpDatabase.
   *********************/
  DefLocalMapping::~DefLocalMapping() { delete warpDB_; }

  /*********************************
   * Run. This function is to run in parallel the deformable 
   * mapping. It is an infinite loop that runs the deformable  
   * mapping in case of a new keyframe
   *********************************/
  void DefLocalMapping::Run()
  {

    mbFinished = false;

    while (1)
    {
      SetAcceptKeyFrames(false);
      insideTheLoop();
      if (Stop())
      {
        // Safe area to stop
        while (isStopped() && !CheckFinish())
        {
          usleep(3000);
        }
        if (CheckFinish())
          break;
      }

      ResetIfRequested();

      if (CheckFinish())
        break;

      SetAcceptKeyFrames(true);

      usleep(100000);
    }

    SetFinish();
  }

  /*********************************
   * insideTheLoop(). This function runs the deformable 
   * mapping. It creates the 
   *********************************/
  void DefLocalMapping::insideTheLoop()
  {
    // Tracking will see that Local Mapping is busy
    // Check if there are keyframes in the queue
    if (CheckNewKeyFrames())
    {
      // BoW conversion and insertion in Map
      this->ProcessNewKeyFrame();
      // Remove points duplicated and not visible
      this->MapPointCulling();
      // Use NRSfM to create or refine surfaces.
      this->NRSfM();

      mbAbortBA = false;
    }
  }

  /*********************************
   * UpdateTemplate(). This function is called in the deformable 
   * tracking and updates the template with the reference keyframe 
   * selected. It creates the new map points with the surface and 
   * create the templates 
   *********************************/
  bool DefLocalMapping::updateTemplate()
  {
    unique_lock<mutex> lock(mMutexReset);
    if ((createTemplate_) & (!stopRequested()))
    {
      std::unique_lock<std::mutex> M(
          static_cast<DefMap *>(mpMap)->MutexUpdating);
      static_cast<DefMap *>(mpMap)->clearTemplate();
      this->CreateNewMapPoints();
      static_cast<DefMap *>(mpMap)->createTemplate(referenceKF_);
      static_cast<DefKeyFrame *>(referenceKF_)->assignTemplate();
      createTemplate_ = false;
      return true;
    }
    return false;
  }

  /******
   * ProcessNewKeyframe. It add the keyframe to the warp database
   * to extract its warp with its covisible. In adition it add the
   * keyframe to a spanning covisibility tree as in ORBSLAM2:LocalMapping. 
   ******/
  void DefLocalMapping::ProcessNewKeyFrame()
  {
    LocalMapping::ProcessNewKeyFrame();
    warpDB_->add(mpCurrentKeyFrame);
  }

  /*************************
   * NRSfM. Given the warp and its derivates it estimates the normals.
   * With those normals it computes an up-to-scale surface. It align that
   * surface to recover the scale with the map points pose registered when
   * the keyframe was created
   *****************************/
  void DefLocalMapping::NRSfM()
  {
    // When there is more than two keyframes, the non-rigid reconstruction starts.

    {
      std::cout << "NORMAL ESTIMATOR IN - ";
      NormalEstimator NormalEstimator_(this->warpDB_);
      NormalEstimator_.ObtainK1K2();
      std::cout << " NORMAL ESTIMATOR OUT";
    }
    if (mpMap->GetAllKeyFrames().size() == 1)
    {
      referenceKF_ = mpMap->GetAllKeyFrames()[0];
    }

    auto newTemplate = needNewTemplate();
    KeyFrame *kfForTemplate;
    if (newTemplate)
    {
      printf("New template requested \n");
      kfForTemplate = mpCurrentKeyFrame;
    }
    else
    {
      kfForTemplate = selectKeyframe();
    }

    // Check if there are enough normals to do the shape-from-normals
    if (!static_cast<DefKeyFrame *>(kfForTemplate)
             ->surface->enoughNormals())
    {
      printf("Not enough normals \n");
      return;
    }
    {
      ShapeFromNormals SfN(kfForTemplate, bendingReg_);

      // Integrate the normals to recover the surface
      if (!SfN.estimate())
      {
        printf("ShapeFromNormals not sucessful \n");
        return;
      }
    }
    if (saveResults_)
    {
      float scale =
          static_cast<GroundTruthKeyFrame *>(kfForTemplate)->estimateAngleErrorAndScale();
      std::cout << "Scale Error Keyframe : " << scale << " " << std::endl;
    }
    if (kfForTemplate != mpMap->GetAllKeyFrames()[0])
    {
      SurfaceRegistration SurfaceRegistration_(kfForTemplate, chiLimit_, true);
      bool wellRegistered = SurfaceRegistration_.registerSurfaces();
      if (!wellRegistered)
      {
        printf("SurfaceRegistration not sucessful (Not enough points to align or chi2 too big \n");
        return;
      }
    }
    referenceKF_ = kfForTemplate;
    createTemplate_ = true;
  }

  /*********************************
   * Create the new map points. They are extracted from the surface 
   * estimated for the keyframe with the Isometric NRSfM.
   ********************************/
  void DefLocalMapping::CreateNewMapPoints()
  {
    cv::Mat Twc = referenceKF_->GetPoseInverse();
    size_t nval = referenceKF_->mvKeysUn.size();
    int cols = referenceKF_->imGray.cols;
    cv::Mat mask(referenceKF_->imGray.rows, referenceKF_->imGray.cols,
                 CV_8UC1, cv::Scalar(0));
    int const max_BINARY_value = 255;

    for (size_t i = 0; i < nval; i++)
    {
      MapPoint *pMP = referenceKF_->GetMapPoint(i);
      if (pMP)
      {
        if (pMP->isBad())
          continue;
        mask.at<char>(referenceKF_->mvKeysUn[i].pt.y,
                      referenceKF_->mvKeysUn[i].pt.x) = 255;
      }
    }
    cv::Mat kernel;
    int kernel_size = cols / 20;
    int ddepth = -1;
    cv::Point anchor(-1, -1);
    double delta;
    delta = 0;
    kernel = cv::Mat::ones(kernel_size, kernel_size, CV_32F);

    cv::filter2D(mask, mask, ddepth, kernel, anchor, delta, cv::BORDER_DEFAULT);
    double threshold_value = 1;
    //     1: Binary Inverted
    cv::threshold(mask, mask, threshold_value, max_BINARY_value, 0);
    uint newPoints(0);
    for (size_t i = 0; i < nval; i++)
    {
      MapPoint *pMP = referenceKF_->GetMapPoint(i);

      if (pMP)
      {
        if (pMP->isBad())
          continue;
        DefMapPoint *defMP = static_cast<DefMapPoint *>(pMP);
        cv::Vec3f x3c;

        static_cast<DefKeyFrame *>(referenceKF_)
            ->surface->get3DSurfacePoint(i, x3c);

        cv::Mat x3ch(4, 1, CV_32F);
        x3ch.at<float>(0, 0) = x3c(0);
        x3ch.at<float>(1, 0) = x3c(1);
        x3ch.at<float>(2, 0) = x3c(2);
        x3ch.at<float>(3, 0) = 1;

        cv::Mat x3wh(4, 1, CV_32F);
        x3wh = Twc * x3ch;
        cv::Mat x3w(3, 1, CV_32F);
        x3w.at<float>(0, 0) = x3wh.at<float>(0, 0);
        x3w.at<float>(1, 0) = x3wh.at<float>(1, 0);
        x3w.at<float>(2, 0) = x3wh.at<float>(2, 0);
        cv::Mat X3Do = static_cast<DefMapPoint *>(pMP)
                           ->PosesKeyframes[referenceKF_]
                           .clone();

        if (X3Do.empty())
        {
          pMP->SetWorldPos(x3w);
          defMP->lastincorporasion = false;
          continue;
        }
        pMP->SetWorldPos(x3w);
      }
      else
      {
        const auto &kpt = referenceKF_->mvKeysUn[i].pt;
        if (mask.at<char>(kpt.y, kpt.x))
        {
          continue;
        }
        cv::Vec3f x3c;
        static_cast<DefKeyFrame *>(referenceKF_)
            ->surface->get3DSurfacePoint(i, x3c);

        cv::Mat x3ch(4, 1, CV_32F);
        x3ch.at<float>(0, 0) = x3c(0);
        x3ch.at<float>(1, 0) = x3c(1);
        x3ch.at<float>(2, 0) = x3c(2);
        x3ch.at<float>(3, 0) = 1;

        cv::Mat x3wh(4, 1, CV_32F);
        x3wh = Twc * x3ch;
        cv::Mat x3w(3, 1, CV_32F);
        x3w.at<float>(0, 0) = x3wh.at<float>(0, 0);
        x3w.at<float>(1, 0) = x3wh.at<float>(1, 0);
        x3w.at<float>(2, 0) = x3wh.at<float>(2, 0);

        pMP = new DefMapPoint(x3w, referenceKF_, mpMap);

        pMP->AddObservation(referenceKF_, i);
        referenceKF_->addMapPoint(pMP, i);

        pMP->ComputeDistinctiveDescriptors();
        pMP->UpdateNormalAndDepth();
        mpMap->addMapPoint(pMP);
        mlpRecentAddedMapPoints.push_back(pMP);
      }
    }
    std::cout << "Points Created : " << newPoints << std::endl;
  }

  /*********************************
   * This function evaluates if the algorithm is visiting new zones. It
   * measures the ocupancy of the map points in the image with a mask.
   * Then, it counts how many keypoints are out of the masked zone to initialize.
   * If it goes over a certain threshold set up in the contructor, it is exploring.
   ********************************/
  bool DefLocalMapping::needNewTemplate()
  {
    int nval = this->mpCurrentKeyFrame->mvKeysUn.size();
    int cols = mpCurrentKeyFrame->imGray.cols;
    cv::Mat mask(mpCurrentKeyFrame->imGray.rows, mpCurrentKeyFrame->imGray.cols,
                 CV_8UC1, cv::Scalar(0));
    int const max_BINARY_value = 255;

    for (int i = 0; i < nval; i++)
    {
      MapPoint *pMP = mpCurrentKeyFrame->GetMapPoint(i);
      if (pMP)
      {
        if (pMP->isBad())
          continue;
        mask.at<char>(this->mpCurrentKeyFrame->mvKeysUn[i].pt.y,
                      this->mpCurrentKeyFrame->mvKeysUn[i].pt.x) = 255;
      }
    }
    cv::Mat kernel;
    int kernel_size = cols / 20;
    int ddepth = -1;
    cv::Point anchor(-1, -1);
    double delta;
    delta = 0;
    kernel = cv::Mat::ones(kernel_size, kernel_size, CV_32F);
    cv::filter2D(mask, mask, ddepth, kernel, anchor, delta, cv::BORDER_DEFAULT);
    double threshold_value = 1;
    cv::threshold(mask, mask, threshold_value, max_BINARY_value, 0);

    int newPoints(0);
    for (int i = 0; i < nval; i++)
    {
      MapPoint *pMP = mpCurrentKeyFrame->GetMapPoint(i);
      if (!pMP)
      {
        const auto &kpt = mpCurrentKeyFrame->mvKeysUn[i].pt;
        if (mask.at<char>(kpt.y, kpt.x))
        {
          continue;
        }
        newPoints++;
      }
    }

    bool createNewTemplate = (newPoints > pointsToTemplate_);
    std::cout << "Points potential : " << newPoints << "  " << pointsToTemplate_
              << std::endl;
    return createNewTemplate;
  }

  /***************************************
   *  This function pick up the reference keyframe in case that there is no 
   * exploration. From the observed points we select the keyframe with the highest number
   * of observed map points.
  ****************************************/
  ORB_SLAM2::KeyFrame *DefLocalMapping::selectKeyframe()
  {
    size_t nval = this->mpCurrentKeyFrame->mvKeysUn.size();
    std::unordered_map<KeyFrame *, int> countKFMatches;

    for (size_t i = 0; i < nval; i++)
    {
      MapPoint *pMP = mpCurrentKeyFrame->GetMapPoint(i);
      if (pMP)
      {
        if (pMP->isBad())
          continue;
        KeyFrame *refkfi = pMP->GetReferenceKeyFrame();
        if (countKFMatches.count(refkfi) == 0)
          countKFMatches[refkfi] = 0;
        countKFMatches[refkfi]++;
      }
    }
    KeyFrame *refkfMaxPoints = mpCurrentKeyFrame;
    int CountPoints(0);
    for (auto &k : countKFMatches)
    {
      if (CountPoints < k.second)
      {
        refkfMaxPoints = k.first;
        CountPoints = k.second;
      }
    }

    return refkfMaxPoints;
  }

  // Reset the algorithm if reset bottom is pushed in the Viewer.
  void DefLocalMapping::ResetIfRequested()
  {
    unique_lock<mutex> lock(mMutexReset);
    if (mbResetRequested)
    {
      static_cast<DefMap *>(mpMap)->clear();
      warpDB_->clear();
      createTemplate_ = false;
      mlNewKeyFrames.clear();
      mlpRecentAddedMapPoints.clear();
      mbResetRequested = false;
    }
  }
} // namespace defSLAM
