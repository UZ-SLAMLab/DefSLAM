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

#include "DefTracking.h"
#include "DefOptimizer.h"
#include "GroundTruthFrame.h"
#include "Optimizer.h"

#include "TriangularMesh.h"
#include "DefORBmatcher.h"
#include "GroundTruthKeyFrame.h"
#include "MinMedianFilter.h"
#include "PnPsolver.h"
#include <DefLocalMapping.h>
#include <GroundTruthCalculator.h>

#include <DefMap.h>
#include <DefMapDrawer.h>
#include <set_MAC.h>
namespace defSLAM
{
  class DefLocalMapping;
  class DefMap;
  class DefMapDrawer;
  class Node;

  // Constructor
  DefTracking::DefTracking(System *pSys, ORBVocabulary *pVoc,
                           FrameDrawer *pFrameDrawer,
                           MapDrawer *pMapDrawer, Map *pMap,
                           KeyFrameDatabase *pKFDB,
                           const string &strSettingPath,
                           const int sensor, bool viewerOn)
      : Tracking(pSys, pVoc, pFrameDrawer, pMapDrawer, pMap, pKFDB,
                 strSettingPath, sensor, viewerOn)
  {
    cv::FileStorage fSettings(strSettingPath, cv::FileStorage::READ);
    RegLap = fSettings["Regularizer.laplacian"];
    RegInex = fSettings["Regularizer.Inextensibility"];
    RegTemp = fSettings["Regularizer.temporal"];
    double a = fSettings["Regularizer.LocalZone"];
    double SaveResults = fSettings["Viewer.SaveResults"];
    saveResults = bool(uint(SaveResults));

    cout << endl
         << "Defomation tracking Parameters: " << endl;
    cout << "- Reg. Inextensibility: " << RegInex << endl;
    cout << "- Reg. Laplacian: " << RegLap << endl;
    cout << "- Reg. Temporal: " << RegTemp << endl;
    cout << "- Reg. LocalZone: " << a << endl;

    LocalZone = uint(a);

    ReliabilityThreshold = fSettings["Regularizer.Reliability"];
  }

  // Main function of tracking.
  void DefTracking::Track()
  {
    if (mState == NO_IMAGES_YET)
    {
      mState = NOT_INITIALIZED;
    }

    mLastProcessedState = mState;
    unique_lock<mutex> lock(mpMap->mMutexMapUpdate);

    if (mState == NOT_INITIALIZED)
    {
      this->MonocularInitialization();
      if (mState != OK)
        return;
    }
    else
    {
      // System is initialized. Track Frame.
      bool bOK;
      // Initial camera pose estimation using motion model or relocalization (if
      // tracking is lost)
      if (!mbOnlyTracking)
      {
        bOK = TrackWithMotionModel();
        mCurrentFrame->mpReferenceKF = mpReferenceKF;
        // If we have an initial estimation of the camera pose and matching. Track
        // the local map.
        if (bOK)
        {

          if ((static_cast<DefLocalMapping *>(mpLocalMapper)
                   ->updateTemplate()))
          {
            mpReferenceKF =
                static_cast<DefMap *>(mpMap)->GetTemplate()->kf;

            Optimizer::DefPoseOptimization(
                mCurrentFrame, mpMap, this->getRegLap(), this->getRegInex(), 0,
                LocalZone);

            if (viewerOn)
            {
              static_cast<DefMapDrawer *>(mpMapDrawer)
                  ->updateTemplateAtRest();
            }
          }
          bOK = TrackLocalMap();
        }
      }
      else
      {

        bOK = this->OnlyLocalisation();
        // If we have an initial estimation of the camera pose and matching. Track
        // the local map.
        mCurrentFrame->mpReferenceKF = mpReferenceKF;

        // mbVO true means that there are few matches to MapPoints in the map. We
        // cannot retrieve
        // a local map and therefore we do not perform TrackLocalMap(). Once the
        // system relocalizes
        // the camera we will use the local map again.
        if (bOK && !mbVO)
          bOK = TrackLocalMap();
      }

      if (bOK)
      {
        mState = OK;

        // Update motion model
        if (!mLastFrame.mTcw.empty())
        {
          cv::Mat LastTwc = cv::Mat::eye(4, 4, CV_32F);
          mLastFrame.GetRotationInverse().copyTo(
              LastTwc.rowRange(0, 3).colRange(0, 3));
          mLastFrame.GetCameraCenter().copyTo(LastTwc.rowRange(0, 3).col(3));
          mVelocity = mCurrentFrame->mTcw * LastTwc;
        }
        else
          mVelocity = cv::Mat();

        if (viewerOn)
        {
          mpMapDrawer->SetCurrentCameraPose(mCurrentFrame->mTcw);
          mpFrameDrawer->Update(this);
          mpMapDrawer->UpdatePoints(mCurrentFrame);
          static_cast<DefMapDrawer *>(mpMapDrawer)->updateTemplate();
        }

        // Clean VO matches
        CleanMatches();
        // Erase Temporal points;
        EraseTemporalPoints();

// Check if we need to insert a new keyframe
        if ((mCurrentFrame->mnId % 10) < 1)
        {
          this->CreateNewKeyFrame();
        }
        
        // We allow points with high innovation (considererd outliers by the Huber
        // Function) pass to the new keyframe, so that bundle adjustment will
        // finally decide if they are outliers or not. We don't want next frame to
        // estimate its position with those points so we discard them in the
        // frame.
        for (size_t i = 0; i < size_t(mCurrentFrame->N); i++)
        {
          if (mCurrentFrame->mvpMapPoints[i] && mCurrentFrame->mvbOutlier[i])
          {
            mCurrentFrame->mvpMapPoints[i] = static_cast<MapPoint *>(nullptr);
          }
        }
      }

      else
      {
        mState = LOST;
        this->status << mCurrentFrame->mTimeStamp << " " << 1 << std::endl;
        if (viewerOn)
        {
          mpFrameDrawer->Update(this);
          mpMapDrawer->UpdatePoints(mCurrentFrame);
        }

        cout << "Track lost soon after initialisation, reseting..." << endl;
        mpSystem->Reset();
        return;
      }

      if (!mCurrentFrame->mpReferenceKF)
        mCurrentFrame->mpReferenceKF = mpReferenceKF;
      mLastFrame = Frame(*mCurrentFrame);
    }

    if (!mCurrentFrame->mTcw.empty())
    {
      cv::Mat Tcr =
          mCurrentFrame->mTcw * mCurrentFrame->mpReferenceKF->GetPoseInverse();
      mlRelativeFramePoses.push_back(Tcr);
      mlpReferences.push_back(mpReferenceKF);
      mlFrameTimes.push_back(mCurrentFrame->mTimeStamp);
      mlbLost.push_back(mState == LOST);
    }
    else
    {
      // This can happen if tracking is lost
      mlRelativeFramePoses.push_back(mlRelativeFramePoses.back());
      mlpReferences.push_back(mlpReferences.back());
      mlFrameTimes.push_back(mlFrameTimes.back());
      mlbLost.push_back(mState == LOST);
    }
  }

  //Main function of tracking where the map is considered deformable.
  bool DefTracking::TrackLocalMap()
  {
    // We have an estimation of the camera pose and some map points tracked in the
    // frame. We retrieve the local map and try to find matches to points in the
    // local map.
    UpdateLocalMap();
    SearchLocalPoints();
    // Optimize with deformable function if there is a template.
    if (static_cast<DefMap *>(mpMap)->GetTemplate())
    {
      Optimizer::DefPoseOptimization(
          mCurrentFrame, mpMap, this->getRegLap(), this->getRegInex(),
          this->getRegTemp(), LocalZone);
    }
    else
    {
      ORB_SLAM2::Optimizer::poseOptimization(mCurrentFrame);
    }

    // Count inliers, outliers and make statistics for map point culling.
    mnMatchesInliers = 0;
    int mnMatchesOutliers(0);
    int DefnToMatchLOCAL(0);
    for (int i = 0; i < mCurrentFrame->N; i++)
    {
      if (mCurrentFrame->mvpMapPoints[i])
      {
        if (!mCurrentFrame->mvbOutlier[i])
        {
          mCurrentFrame->mvpMapPoints[i]->IncreaseFound();
          if (!mbOnlyTracking)
          {
            if (mCurrentFrame->mvpMapPoints[i]->Observations() > 0)
            {
              mnMatchesInliers++;
              if (static_cast<DefMapPoint *>(
                      mCurrentFrame->mvpMapPoints[i])
                      ->getFacet())
                DefnToMatchLOCAL++;
            }
          }
          else
            mnMatchesInliers++;
        }
        else
        {
          mnMatchesOutliers++;
        }
      }
    }
    auto points = mpMap->GetReferenceMapPoints();
    auto numberLocalMapPoints(0);
    for (auto pMP : points)
    {
      if (pMP)
      {
        if (pMP->isBad())
          continue;
        if (static_cast<DefMapPoint *>(pMP)->getFacet())
          if (mCurrentFrame->isInFrustum(pMP, 0.5))
          {
            numberLocalMapPoints++;
          }
      }
    }
    // Optimize Pose
    int observedFrame(0);
    int mI(0);
    int mO(0);
    for (int i = 0; i < mCurrentFrame->N; i++)
    {
      if (mCurrentFrame->mvpMapPoints[i])
      {
        if (mCurrentFrame->mvpMapPoints[i]->isBad())
          continue;
        observedFrame++;
        if (!mCurrentFrame->mvbOutlier[i])
        {
          mI++;
        }
        else
        {
          mO++;
        }
      }
    }

    // Save matching result
    std::ostringstream out;
    out << std::internal << std::setfill('0') << std::setw(5)
        << uint(mCurrentFrame->mTimeStamp);
    std::cout << out.str() << " " << mI << " " << mO << " "
              << numberLocalMapPoints << std::endl;
    this->matches << out.str() << " " << mI << " " << mO << " "
                  << numberLocalMapPoints << std::endl;
    // Decide if the tracking was succesful
    // More restrictive if there was a relocalization recently
    if (mCurrentFrame->mnId < mnLastRelocFrameId + mMaxFrames &&
        mnMatchesInliers < 20)
      return false;

    if (mnMatchesInliers < 10)
      return false;
    else
      return true;
  }

  // Initial tracking to locate rigidly the camera and discard outliers.
  bool DefTracking::TrackWithMotionModel()
  {
    DefORBmatcher Defmatcher(0.9, false);

    // Update last frame pose according to its reference keyframe
    // Create "visual odometry" points if in Localization Mode
    UpdateLastFrame();

    mCurrentFrame->SetPose(mLastFrame.mTcw);

    fill(mCurrentFrame->mvpMapPoints.begin(), mCurrentFrame->mvpMapPoints.end(),
         static_cast<MapPoint *>(nullptr));

    // Project points seen in previous frame
    int th(20);

    int nmatches = Defmatcher.SearchByProjection(*mCurrentFrame, mLastFrame, th,
                                                 mSensor == System::MONOCULAR);

    std::cout << "POINTS matched:" << nmatches << std::endl;

    // If few matches, uses a wider window search
    if (nmatches /* +nmatches2 */ < 20)
    {
      fill(mCurrentFrame->mvpMapPoints.begin(), mCurrentFrame->mvpMapPoints.end(),
           static_cast<MapPoint *>(nullptr));
      nmatches = Defmatcher.SearchByProjection(*mCurrentFrame, mLastFrame, th + 5,
                                               mSensor == System::MONOCULAR);
    }
    // std::cout << "mnMatches 2 Motion: " << nmatches << std::endl;

    if (nmatches < 15)
      return false;
    return true;

    /// Optimize frame pose with all matches with a rigid model to initialize the
    /// pose of the camera
    Optimizer::poseOptimization(mCurrentFrame, myfile);

    // Discard outliers
    int nmatchesMap = 0;
    for (int i = 0; i < mCurrentFrame->N; i++)
    {
      if (mCurrentFrame->mvpMapPoints[i])
      {
        if (mCurrentFrame->mvbOutlier[i])
        {
          MapPoint *pMP = mCurrentFrame->mvpMapPoints[i];
          mCurrentFrame->mvpMapPoints[i] = static_cast<MapPoint *>(nullptr);
          mCurrentFrame->mvbOutlier[i] = false;
          pMP->mbTrackInView = false;
          pMP->mnLastFrameSeen = mCurrentFrame->mnId;
          nmatches--;
        }
        else if (mCurrentFrame->mvpMapPoints[i]->Observations() > 0)
          nmatchesMap++;
      }
    }

    if (mbOnlyTracking)
    {
      mbVO = nmatchesMap < 10;
      return nmatches > 20;
    }

    return nmatchesMap >= 10;
  }

  // Update the last frame relative pose.
  void DefTracking::UpdateLastFrame()
  {
    // Update pose according to reference keyframe
    KeyFrame *pRef = mLastFrame.mpReferenceKF;
    cv::Mat Tlr = mlRelativeFramePoses.back();

    mLastFrame.SetPose(Tlr * pRef->GetPose());

    if (mnLastKeyFrameId == mLastFrame.mnId || mSensor == System::MONOCULAR ||
        !mbOnlyTracking)
      return;
  }

  // Check local points with the covisible keyframes.
  // Check there are no repeated points.
  void DefTracking::UpdateLocalPoints()
  {
    std::set<MapPoint *> mPall;
    for (vector<KeyFrame *>::const_iterator itKF = mvpLocalKeyFrames.begin(),
                                            itEndKF = mvpLocalKeyFrames.end();
         itKF != itEndKF; itKF++)
    {
      KeyFrame *pKF = *itKF;
      const vector<MapPoint *> vpMPs = pKF->GetMapPointMatches();
      for (vector<MapPoint *>::const_iterator itMP = vpMPs.begin(),
                                              itEndMP = vpMPs.end();
           itMP != itEndMP; itMP++)
      {
        MapPoint *pMP = *itMP;
        if (!pMP)
          continue;
        if (pMP->mnTrackReferenceForFrame == mCurrentFrame->mnId)
          continue;
        if (!pMP->isBad())
        {
          pMP->mnTrackReferenceForFrame = mCurrentFrame->mnId;
          mPall.insert(pMP);
        }
      }
    }
    mvpLocalMapPoints.clear();
    mvpLocalMapPoints.resize(mPall.size());
    std::copy(mPall.begin(), mPall.end(), mvpLocalMapPoints.begin());
  }

  // Run for image with ground truth through a stereo pair
  cv::Mat DefTracking::GrabImageMonocularGT(const cv::Mat &imRectLeft,
                                            const cv::Mat &imRectRight,
                                            const double &timestamp,
                                            cv::Mat _mask)
  {
    mImGray = imRectLeft.clone();
    cv::Mat imGrayRight = imRectRight;
    imRectLeft.copyTo(mImRGB);

    if (mImGray.channels() == 3)
    {
      if (mbRGB)
      {
        cv::cvtColor(mImGray, mImGray, cv::COLOR_RGB2GRAY);
        cv::cvtColor(imGrayRight, imGrayRight, cv::COLOR_RGB2GRAY);
      }
      else
      {
        cv::cvtColor(mImGray, mImGray, cv::COLOR_BGR2GRAY);
        cv::cvtColor(imGrayRight, imGrayRight, cv::COLOR_BGR2GRAY);
      }
    }
    else if (mImGray.channels() == 4)
    {
      if (mbRGB)
      {
        cv::cvtColor(mImGray, mImGray, cv::COLOR_RGBA2GRAY);
        cv::cvtColor(imGrayRight, imGrayRight, cv::COLOR_RGBA2GRAY);
      }
      else
      {
        cv::cvtColor(mImGray, mImGray, cv::COLOR_BGRA2GRAY);
        cv::cvtColor(imGrayRight, imGrayRight, cv::COLOR_BGRA2GRAY);
      }
    }
    else
    {
      cv::cvtColor(imRectLeft, mImRGB, cv::COLOR_GRAY2RGB);
    }

    mCurrentFrame = new GroundTruthFrame(
        mImGray, timestamp, mpORBextractorLeft, mpORBVocabulary, mK, mDistCoef,
        mbf, mThDepth, imRectLeft, imGrayRight, _mask);

    this->Track();

    if ((mState == eTrackingState::OK) && (saveResults))
    {
      float scale =
          static_cast<GroundTruthFrame *>(mCurrentFrame)->Estimate3DScale(mpMap);
      scalefile << mCurrentFrame->mTimeStamp << " " << scale << std::endl;
      double error = static_cast<GroundTruthFrame *>(mCurrentFrame)
                         ->Estimate3DError(mpMap, scale);

      if (viewerOn)
      {
        mpMapDrawer->UpdatePoints(mCurrentFrame, scale);
        this->mpFrameDrawer->SetError(error);
      }
    }

    return mCurrentFrame->mTcw.clone();
  }

  // Run for image with ground truth through a depth image. It is used for
  // the CT image of the phantom datasets, but it should be usable for rgbd image.
  cv::Mat DefTracking::GrabImageMonocularCTGT(const cv::Mat &imRectLeft,
                                              const cv::Mat &imDepth,
                                              const double &timestamp,
                                              cv::Mat _mask)
  {
    mImGray = imRectLeft.clone();
    imRectLeft.copyTo(mImRGB);

    if (mImGray.channels() == 3)
    {
      if (mbRGB)
      {
        cv::cvtColor(mImGray, mImGray, cv::COLOR_RGB2GRAY);
      }
      else
      {
        cv::cvtColor(mImGray, mImGray, cv::COLOR_BGR2GRAY);
      }
    }
    else if (mImGray.channels() == 4)
    {
      if (mbRGB)
      {
        cv::cvtColor(mImGray, mImGray, cv::COLOR_RGBA2GRAY);
      }
      else
      {
        cv::cvtColor(mImGray, mImGray, cv::COLOR_BGRA2GRAY);
      }
    }
    else
    {
      cv::cvtColor(imRectLeft, mImRGB, cv::COLOR_GRAY2RGB);
    }

    mCurrentFrame = new GroundTruthFrame(
        mImGray, timestamp, mpORBextractorLeft, mpORBVocabulary, mK, mDistCoef,
        mbf, mThDepth, imRectLeft, imDepth, true, _mask);

    this->Track();

    if ((mState == eTrackingState::OK) && (saveResults))
    {
      float scale =
          static_cast<GroundTruthFrame *>(mCurrentFrame)->Estimate3DScale(mpMap);
      scalefile << mCurrentFrame->mTimeStamp << " " << scale << std::endl;
      double error = static_cast<GroundTruthFrame *>(mCurrentFrame)
                         ->Estimate3DError(mpMap, scale);

      if (viewerOn)
      {
        mpMapDrawer->UpdatePoints(mCurrentFrame, scale);
        this->mpFrameDrawer->SetError(error);
      }
      // GroundTruthCalculator::CompareDefGT(&mCurrentFrame,timestamp);
    }
    return mCurrentFrame->mTcw.clone();
  }

  // Initialize scene with a plane.
  void DefTracking::MonocularInitialization()
  {
    /// Initialize the surface and the points in the surface considering a plane
    /// parallel to the camera plane
    if (mCurrentFrame->N > 100)
    {
      // Set Frame pose to the origin
      mCurrentFrame->SetPose(cv::Mat::eye(4, 4, CV_32F));

      // Create KeyFrame
      KeyFrame *pKFini =
          new GroundTruthKeyFrame(*mCurrentFrame, mpMap, mpKeyFrameDB);

      // Insert KeyFrame in the map
      mpMap->AddKeyFrame(pKFini);
      // Create MapPoints with inliers and associate to KeyFrame
      for (size_t i = 0; i < size_t(mCurrentFrame->N); i++)
      {
        cv::KeyPoint kp = mCurrentFrame->mvKeysUn[i];

        cv::Mat x3D = (cv::Mat_<float>(3, 1)
                           << (kp.pt.x - mCurrentFrame->cx) / mCurrentFrame->fx,
                       (kp.pt.y - mCurrentFrame->cy) / mCurrentFrame->fy, 1);
        MapPoint *pNewMP = new DefMapPoint(x3D, pKFini, mpMap);

        pNewMP->AddObservation(pKFini, i);
        pKFini->addMapPoint(pNewMP, i);
        pNewMP->ComputeDistinctiveDescriptors();
        pNewMP->UpdateNormalAndDepth();
        mpMap->addMapPoint(pNewMP);
        mCurrentFrame->mvpMapPoints[i] = pNewMP;
        pKFini->GetMapPoint(i);
      }

      double *Array;
      Array = new double[_NumberOfControlPointsU * _NumberOfControlPointsV];

      for (uint i(0); i < _NumberOfControlPointsU * _NumberOfControlPointsV;
           i++)
      {
        Array[i] = 1;
      }
      BBS::bbs_t bbs;
      auto defkf = static_cast<DefKeyFrame *>(pKFini);
      bbs.umin = defkf->umin;
      bbs.umax = defkf->umax;
      bbs.nptsu = _NumberOfControlPointsU;
      bbs.nptsv = _NumberOfControlPointsV;
      bbs.vmax = defkf->vmax;
      bbs.vmin = defkf->vmin;
      bbs.valdim = 1;
      defkf->surface->saveArray(Array, bbs);

      cout << "New map created with " << mpMap->MapPointsInMap() << " points"
           << endl;
      mLastFrame = Frame(*mCurrentFrame);
      mnLastKeyFrameId = uint(mCurrentFrame->mnId);
      mpLastKeyFrame = pKFini;
      mvpLocalKeyFrames.push_back(pKFini);
      mvpLocalMapPoints = mpMap->GetAllMapPoints();
      mpReferenceKF = pKFini;
      mCurrentFrame->mpReferenceKF = pKFini;
      mLastFrame.mpReferenceKF = pKFini;
      mpMap->SetReferenceMapPoints(mvpLocalMapPoints);
      mpMap->mvpKeyFrameOrigins.push_back(pKFini);
      mpLocalMapper->InsertKeyFrame(pKFini);
      cv::Mat Tcr =
          mCurrentFrame->mTcw * mCurrentFrame->mpReferenceKF->GetPoseInverse();
      mlRelativeFramePoses.push_back(Tcr);
      // Initialize the SLAM
      static_cast<DefKeyFrame *>(pKFini)->assignTemplate();
      static_cast<DefMap *>(mpMap)->createTemplate(pKFini);
      if (viewerOn)
      {
        mpMapDrawer->SetCurrentCameraPose(mCurrentFrame->mTcw);
        mpFrameDrawer->Update(this);
        mpMapDrawer->UpdatePoints(mCurrentFrame);
        static_cast<DefMapDrawer *>(mpMapDrawer)->updateTemplate();
        static_cast<DefMapDrawer *>(mpMapDrawer)->updateTemplateAtRest();
      }
      mState = OK;
    }
  }

  // Remove matches of the current keyframe.
  void DefTracking::CleanMatches()
  {
    for (int i = 0; i < mCurrentFrame->N; i++)
    {
      MapPoint *pMP = mCurrentFrame->mvpMapPoints[i];
      if (pMP)
        if (pMP->Observations() < 1)
        {
          mCurrentFrame->mvbOutlier[i] = false;
          mCurrentFrame->mvpMapPoints[i] = static_cast<MapPoint *>(nullptr);
        }
    }
  }

  // Erase temporal points.
  void DefTracking::EraseTemporalPoints()
  {
    // Delete temporal MapPoints
    for (list<MapPoint *>::iterator lit = mlpTemporalPoints.begin(),
                                    lend = mlpTemporalPoints.end();
         lit != lend; lit++)
    {
      MapPoint *pMP = *lit;
      delete pMP;
    }
    mlpTemporalPoints.clear();
  }

  // Create new keyframe for deformable SLAM.
  void DefTracking::CreateNewKeyFrame()
  {
    if (!mpLocalMapper->SetNotStop(true))
      return;

    KeyFrame *pKF = new GroundTruthKeyFrame(*mCurrentFrame, mpMap, mpKeyFrameDB);
    mCurrentFrame->mpReferenceKF = mpReferenceKF;
    mpLocalMapper->InsertKeyFrame(pKF);
    mpLocalMapper->SetNotStop(false);
    mnLastKeyFrameId = mCurrentFrame->mnId;
    mpLastKeyFrame = pKF;
  }
} // namespace defSLAM
