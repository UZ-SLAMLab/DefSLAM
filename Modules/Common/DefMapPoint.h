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

#ifndef DEFMAPPOINT_H
#define DEFMAPPOINT_H
#include "MapPoint.h"

#include "DefMap.h"

#include "Facet.h"
#include <mutex>
#include <opencv2/core/core.hpp>

namespace ORB_SLAM2
{
  class KeyFrame;
  class Map;
  class Frame;
  class MapPoint;
} // namespace ORB_SLAM2

namespace defSLAM
{
  using ORB_SLAM2::Frame;
  using ORB_SLAM2::KeyFrame;
  using ORB_SLAM2::Map;
  using ORB_SLAM2::MapPoint;

  class DefMapPoint : public MapPoint
  {
  public:
    DefMapPoint() = delete;

    /// Constructor. See constructor Map Point
    DefMapPoint(const cv::Mat &Pos, KeyFrame *pRefKF, Map *pMap);

    /// Constructor. See constructor Map Point
    DefMapPoint(const cv::Mat &Pos, Map *pMap, Frame *pFrame,
                const int &idxF);

    /// Constructor by default
    virtual ~DefMapPoint() = default;

    /// Remove this map point from the template temp.
    void RemoveTemplate();

    /// Get position of the point at shape-at-rest.
    cv::Mat GetWorldPosAtRest();

    /// Delete map point.
    void setBadFlag() override;

    /// Assign facet in the template.
    void SetFacet(Facet *face);

    /// Get facet in the template.
    Facet *getFacet();

    /// Set barycentric coordinates in the template.
    void SetCoordinates(double, double, double);

    /// Once assigned the baricentric recalculate its position
    /// and its normals
    void Repose();

    // Set 3D from barycentric
    void RecalculatePosition();

    // Set world position in keyframe kf.
    void SetWorldPos(cv::Mat, KeyFrame *);

  public:
    // Register of the points 3D pose when the keyframe was initialized
    std::unordered_map<KeyFrame *, cv::Mat> PosesKeyframes;

    // Barycentrics
    double b1, b2, b3;

  protected:
    // Template
    Facet *facet;

  protected:
    std::mutex MutexFacet;
    std::mutex MutexNFace;
  };

} // namespace defSLAM

#endif // DEFORMATION MAPPOINT_H
