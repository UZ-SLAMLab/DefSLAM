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

#include "DefMapPoint.h"
#include "ORBmatcher.h"

#include <mutex>

namespace defSLAM
{
  class DiffProp;
  /// Constructor. See constructor Map Point
  DefMapPoint::DefMapPoint(const cv::Mat &Pos, KeyFrame *pRefKF,
                           Map *pMap)
      : MapPoint(Pos, pRefKF, pMap),
        facet(static_cast<Facet *>(nullptr))
  {
  }

  /// Constructor. See constructor Map Point
  DefMapPoint::DefMapPoint(const cv::Mat &Pos, Map *pMap,
                           Frame *pFrame, const int &idxF)
      : MapPoint(Pos, pMap, pFrame, idxF),
        facet(static_cast<Facet *>(nullptr)) {}

  /// Remove this map point from the template temp.
  void DefMapPoint::RemoveTemplate()
  {
    unique_lock<mutex> lock1(mMutexFeatures);
    unique_lock<mutex> lock2(mMutexPos);
    facet = static_cast<Facet *>(nullptr);
  }

  /// Get position of the point at shape-at-rest.
  cv::Mat DefMapPoint::GetWorldPosAtRest()
  {
    std::set<Node *> Nodes = facet->getNodes();
    std::vector<Node *> Nodesvector;
    for (std::set<Node *>::iterator itn = Nodes.begin(); itn != Nodes.end();
         itn++)
    {
      Nodesvector.push_back(*itn);
    }
    cv::Mat mWorldPosAtRest = mWorldPos.clone();
    mWorldPosAtRest.at<float>(0) = b1 * Nodesvector[0]->xO +
                                   b2 * Nodesvector[1]->xO +
                                   b3 * Nodesvector[2]->xO;
    mWorldPosAtRest.at<float>(1) = b1 * Nodesvector[0]->yO +
                                   b2 * Nodesvector[1]->yO +
                                   b3 * Nodesvector[2]->yO;
    mWorldPosAtRest.at<float>(2) = b1 * Nodesvector[0]->zO +
                                   b2 * Nodesvector[1]->zO +
                                   b3 * Nodesvector[2]->zO;
    return mWorldPosAtRest.clone();
  }

  /// Delete map point.
  void DefMapPoint::setBadFlag()
  {
    map<KeyFrame *, size_t> obs;
    {
      unique_lock<mutex> lock1(mMutexFeatures);
      unique_lock<mutex> lock2(mMutexPos);
      mbBad = true;
      obs = mObservations;
      mObservations.clear();
    }
    for (map<KeyFrame *, size_t>::iterator mit = obs.begin(), mend = obs.end();
         mit != mend; mit++)
    {
      KeyFrame *pKF = mit->first;
      pKF->EraseMapPointMatch(mit->second);
    }

    mpMap->eraseMapPoint(this);
  }

  /// Assign facet in the template.
  void DefMapPoint::SetFacet(Facet *face)
  {
    if (face)
      face->addMapPoint(this);
    facet = face;
    std::unique_lock<std::mutex> MutexNce(MutexNFace);
  }

  /// Get facet in the template.
  Facet *DefMapPoint::getFacet()
  {
    std::unique_lock<std::mutex> MutexNce(MutexFacet);
    return facet;
  }

  /// Set barycentric coordinates in the template.
  void DefMapPoint::SetCoordinates(double B1, double B2, double B3)
  {
    b1 = B1;
    b2 = B2;
    b3 = B3;
  }

  /// Once assigned the baricentric recalculate its position
  /// and its normals
  void DefMapPoint::Repose()
  {
    this->RecalculatePosition();
    this->UpdateNormalAndDepth();
  }

  // Set 3D from barycentric
  void DefMapPoint::RecalculatePosition()
  {
    std::set<Node *> Nodes = facet->getNodes();
    std::vector<Node *> Nodesvector;
    for (std::set<Node *>::iterator itn = Nodes.begin(); itn != Nodes.end();
         itn++)
    {
      Nodesvector.push_back(*itn);
    }

    unique_lock<mutex> lock(mMutexPos);

    mWorldPos.at<float>(0) =
        b1 * Nodesvector[0]->x + b2 * Nodesvector[1]->x + b3 * Nodesvector[2]->x;
    mWorldPos.at<float>(1) =
        b1 * Nodesvector[0]->y + b2 * Nodesvector[1]->y + b3 * Nodesvector[2]->y;
    mWorldPos.at<float>(2) =
        b1 * Nodesvector[0]->z + b2 * Nodesvector[1]->z + b3 * Nodesvector[2]->z;
  }

  // Set world position in keyframe kf.
  void DefMapPoint::SetWorldPos(cv::Mat x3d, KeyFrame *kf)
  {
    PosesKeyframes[kf] = x3d;
  }

} // namespace defSLAM
