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

#include "DefKeyFrame.h"
#include "DefMapPoint.h"
#include <opencv2/core/core.hpp>
#include "bbs_MAC.h"
#include <mutex>

namespace defSLAM
{
  /**********************
 * The contructor register the current position of the map points 
 * to align the surface and recover the scale later. It also estimates
 * the retina coordinates.
 * 
 * Arg:
 *  Frame* F -- ptr to the current frame
 *  Map* pMap -- ptr to the map
 *  keyFrameDatabase -- to put the keyframe in the keyframe database. 
 * 
 * It also initialize the surface observed from the keyframe.
 * *******************/
  DefKeyFrame::DefKeyFrame(Frame &F, Map *pMap,
                           KeyFrameDatabase *pKFDB)
      : KeyFrame(F, pMap, pKFDB), Outliers(mvKeysUn.size(), false),
        KeyframesRelated(0), surface(static_cast<Surface *>(nullptr)), accMean(0),
        haveATemplate_(false)
  {
    // Initialize retina coordinates of the keypoints.
    this->umin = 0.75;
    this->umax = -0.75;
    this->vmin = 0.75;
    this->vmax = -0.75;
    this->NCu = _NumberOfControlPointsU;
    this->NCv = _NumberOfControlPointsV;
    this->valdim = 2;
    this->accMean = 1;
    this->NormaliseKeypoints();

    // It register the current position.
    for (size_t i(0); i < mvKeysUn.size(); i++)
    {
      MapPoint *pMP = this->GetMapPoint(i);
      if (pMP)
      {
        if (pMP->isBad())
          continue;

        if (static_cast<DefMapPoint *>(pMP)->getFacet())
        {
          static_cast<DefMapPoint *>(pMP)->PosesKeyframes[this] =
              pMP->GetWorldPos().clone();
        }
      }
    }
    // Initialize surface
    surface = new Surface(mvKeysUn.size());
  }

  //Destructor of the DefKeyFrame
  DefKeyFrame::~DefKeyFrame()
  {
    if (surface)
      delete surface;
    surface = nullptr;
    RGBimage.release();
    imGray.release();
    mpORBvocabulary = nullptr;
    mpKeyFrameDB = nullptr;
  }

  /**********************
 * Normalisation of the points. The keypoints are saved in retina coordinates.
 * *******************/
  void DefKeyFrame::NormaliseKeypoints()
  {
    // Calculate From One To Two
    Eigen::MatrixXf PtNormalizedkp1(3, mvKeysUn.size()),
        Pt2Dkp1(3, mvKeysUn.size());
    for (uint ikp = 0; ikp < mvKeysUn.size(); ikp++)
    {
      const cv::KeyPoint &kp1 = this->mvKeysUn[ikp];
      Pt2Dkp1(0, ikp) = kp1.pt.x;
      Pt2Dkp1(1, ikp) = kp1.pt.y;
      Pt2Dkp1(2, ikp) = 1.0;
    }
    const cv::Mat Kinvcv = this->mK.inv();
    Eigen::Matrix3f Kinv;
    cv::cv2eigen(Kinvcv, Kinv);
    PtNormalizedkp1 << Kinv * Pt2Dkp1;
    mpKeypointNorm.resize(this->mvKeysUn.size());
    for (uint ikp = 0; ikp < mvKeysUn.size(); ikp++)
    {
      mpKeypointNorm[ikp].pt.x = (PtNormalizedkp1(0, ikp));
      mpKeypointNorm[ikp].pt.y = (PtNormalizedkp1(1, ikp));
      const cv::KeyPoint &kp1 = this->mpKeypointNorm[ikp];
      if (kp1.pt.x < this->umin)
      {
        this->umin = kp1.pt.x - 0.10;
      }
      if (kp1.pt.x > this->umax)
      {
        this->umax = kp1.pt.x + 0.10;
      }
      if (kp1.pt.y < this->vmin)
      {
        this->vmin = kp1.pt.y - 0.10;
      }
      if (kp1.pt.y > this->vmax)
      {
        this->vmax = kp1.pt.y + 0.10;
      }
    }
  }
  /************
   * This function define a template as an anchor keyframe in which points 
   * and template are initialized.
   *************************/
  void DefKeyFrame::assignTemplate()
  {
    std::unique_lock<std::mutex> lk(mutex);
    this->haveATemplate_ = true;
  }
  /************
   * This function returns true if this keyframe has a template associated, 
   * i.e., if it is anchor keyframe. It is used just for visualization.
   *************************/
  bool DefKeyFrame::templateAssigned()
  {
    std::unique_lock<std::mutex> lk(mutex);
    return this->haveATemplate_;
  }
} // namespace defSLAM
