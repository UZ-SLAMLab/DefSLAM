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
#ifndef DEFKEYFRAME_H
#define DEFKEYFRAME_H
#include "MapPoint.h"

#include "KeyFrame.h"
#include "Surface.h"

#include <mutex>
#include <opencv2/core/core.hpp>
#include <opencv2/core/eigen.hpp>

namespace ORB_SLAM2
{
  class KeyFrame;
  class Map;
  class Frame;
  class MapPoint;
  class Node;
  class Template;
  class Facet;

} // namespace ORB_SLAM2

namespace defSLAM
{
  using ORB_SLAM2::Frame;
  using ORB_SLAM2::KeyFrame;
  using ORB_SLAM2::KeyFrameDatabase;
  using ORB_SLAM2::Map;
  using ORB_SLAM2::MapPoint;

  class DefKeyFrame : public KeyFrame
  {
  public:
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
    DefKeyFrame(Frame &F, Map *pMap, KeyFrameDatabase *pKFDB);

    //Destructor of the DefKeyFrame
    ~DefKeyFrame();
    /**********************
 * Normalisation of the points. The keypoints are saved in retina coordinates.
 * *******************/
    void NormaliseKeypoints();
    /************
   * This function define a template as an anchor keyframe in which points 
   * and template are initialized.
   *************************/
    void assignTemplate();
    /************
   * This function returns true if this keyframe has a template associated, 
   * i.e., if it is anchor keyframe. It is used just for visualization.
   *************************/
    bool templateAssigned();

  public:
    double umin;
    double umax;
    double vmin;
    double vmax;
    int NCu;
    int NCv;
    int valdim;
    std::vector<cv::KeyPoint> mpKeypointNorm;
    std::vector<bool> Outliers;
    int KeyframesRelated;
    Surface *surface;
    float accMean;

  private:
    std::mutex mutex;
    bool haveATemplate_;
  };

} // namespace defSLAM

#endif // DEFORMATION MAPPOINT_H
