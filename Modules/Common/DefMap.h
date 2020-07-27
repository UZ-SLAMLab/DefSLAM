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
#ifndef DEFMAP_H
#define DEFMAP_H

#include "Map.h"

#include "DefMapPoint.h"
#include "Template.h"
#include <set>

#include <mutex>

namespace ORB_SLAM2
{
  class Map;
  class Template;
  class KeyFrame;
  class MapPoint;
} // namespace ORB_SLAM2

namespace defSLAM
{
  using ORB_SLAM2::KeyFrame;
  using ORB_SLAM2::Map;

  class DefMap : public Map
  {
  public:
    DefMap();

    /***********
   * Constructor:
   *  Function that generates a template with the current points 
   *  and the estimated surface for the keyframe kf. (The surface
   *  must have been estimated.)
   *********/
    DefMap(std::vector<std::vector<double>> &,
           std::vector<std::vector<int>> &);

    /***********
   *  Function that generates a template with the estimated surface for the keyframe kf. 
   * (Kf->surface must have been initializated.)
   *********/
    virtual void createTemplate(KeyFrame *Kf);

    // return the template
    virtual void clear();

    // deletes the current template
    Template *GetTemplate();

    // clear the entire map
    void clearTemplate();

    // Return the reference keyframe or the keyframe used to
    // create the template
    KeyFrame *getLastKFTemplate();

    // Set the reference keyframe or the keyframe used to
    // create the template
    void setLastKFTemplate(KeyFrame *);

  public:
    bool ThereIsATemplate;

    std::mutex MutexUpdating;

  protected:
    Template *Mesh;
    KeyFrame *lastKfTemplate_;
    std::mutex MutexUsingMap;
    std::mutex MutexKf;
  };

} // namespace defSLAM

#endif // DEFORMABLE MAP_H
