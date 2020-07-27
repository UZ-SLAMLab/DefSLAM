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

#ifndef WARPDATABASE_H
#define WARPDATABASE_H

#include "KeyFrame.h"
#include "MapPoint.h"
#include "diffProp.h"

#include <memory>
#include <mutex>

namespace defSLAM
{
  using ORB_SLAM2::KeyFrame;
  using ORB_SLAM2::MapPoint;
  typedef std::vector<std::shared_ptr<DiffProp>> kr2krdata;

  class DiffProp;
  class WarpDatabase
  {
    /*******
     * Base class for the class that estimate the warps and
     * each differential properties. DefSLAM implemented with 
     * Schwarzian warps (See SchwarpDatabase.h).
     * ******/
  public:
    // Constructor by default
    virtual ~WarpDatabase() = default;

    // It insert the keyframe to the database.
    virtual void add(KeyFrame *kf) = 0;

    // Delete kf from database
    virtual void erase(KeyFrame *kf) = 0;

    // clear the database
    virtual void clear() = 0;

    // Get DB of correspondencies and differential properties
    std::map<MapPoint *, kr2krdata> &getDiffDatabase() { return mapPointsDB_; }

    // Check if the correspondencies and differential properties
    // of a map has changed since it was used.
    std::map<MapPoint *, bool> &getToProccess() { return newInformation_; }

  protected:
    std::map<MapPoint *, bool> newInformation_; // Register of new information
                                                // for a map point.

    std::map<MapPoint *, kr2krdata> mapPointsDB_; // DB of correspondencies and differential properties
                                                  // for a pair of keyframes used in normal estimator
  };
} // namespace defSLAM

#endif
