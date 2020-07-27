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

#ifndef SCHWARPDATABASE_H
#define SCHWARPDATABASE_H
#pragma once

#include <list>
#include <set>
#include <vector>

#include "WarpDatabase.h"

#include "KeyFrame.h"
#include "MapPoint.h"

#include <Thirdparty/BBS/bbs_MAC.h>
#include <mutex>

namespace defSLAM
{
  using ORB_SLAM2::Frame;
  using ORB_SLAM2::KeyFrame;
  using ORB_SLAM2::MapPoint;

  class WarpDatabase;
  class DiffProp;

  class SchwarpDatabase : public WarpDatabase
  {
  public:
    // Constructor. lambda is set in the file .yalm
    SchwarpDatabase(double lambda) : lambda_(lambda) {}

    /*********** 
     * It insert the keyframe to the database. There we 
     * search for its covisible anchor keyframes. Once, 
     * we have them, if they have a certain number of matches
     * it estimates a initial warp between the anchor keyframes 
     * and the new keyframe. With that warp, it search for more 
     * matches. Once it is done it recalculate the warp and then
     * it save the correspondences and its derivatives
     * in mapPointsDB_ (WarpDatabase).
     * ******/
    void add(KeyFrame *kf) override;
    // Delete kf from database
    void erase(KeyFrame *kf) override;

    // clear the database
    void clear() override;

  protected:
    /************************
   *  Function that estimates the warps between two keyframes 
   * given the correspondences vMatchesIndices and 
   * saves their correspondances and its derivatives. It uses 
   * the Schwarzian penalization with a lambda weight. 
   * ******************/
    void calculateSchwarps(
        KeyFrame *KF, KeyFrame *KF2,
        vector<pair<size_t, size_t>> &vMatchedIndices,
        double (&x)[_NumberOfControlPointsU * _NumberOfControlPointsV * 2],
        double lambda_);

  protected:
    std::set<KeyFrame *> mpkeyframes; // Keyframes registered
    double lambda_;                   // Schwarzian regularization
  };

} // namespace defSLAM

#endif
