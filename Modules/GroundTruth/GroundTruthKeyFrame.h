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

#ifndef GTKEYFRAME_H
#define GTKEYFRAME_H

#include "DefKeyFrame.h"

#include <mutex>
#include <opencv2/core/core.hpp>
#include <opencv2/core/eigen.hpp>

namespace defSLAM
{

  class DefKeyFrame;

  class GroundTruthKeyFrame : public DefKeyFrame
  {
  public:
    /*****
   *  Constructor. Not by default
   *****/
    GroundTruthKeyFrame() = delete;
    /*****
   *  Constructor. It initializes the parameters for the groundtruth estimation.
   *  TODO: Implement with depth image as in GroundTruthFrame
   *****/
    GroundTruthKeyFrame(Frame &F, Map *pMap, KeyFrameDatabase *pKFDB);

    //Default Destructor
    ~GroundTruthKeyFrame() = default;
    /*****
   * estimateAngleErrorAndScale. This method estimates the 3D groundtruth of the
   * surface estimated for this keyframe. It use the keypoints of the left image 
   * with a normal and search for estimates the 3D in the right image. It uses the 
   * pcl library to determinate the normals of the point cloud and compares them 
   * with the estimated by the NRSfM and the SfN. NRSfM tends to be quite noisy.
   */
    float estimateAngleErrorAndScale();

    // Drawer purposes. Not used now
    std::vector<cv::Point3f> mvLocalMapPoints, mvStereoMapPoints;

  private:
    cv::Mat imRight;

    bool StereoAvaliable;

  private:
    std::vector<std::vector<float>> posMono_;
    std::vector<std::vector<float>> posStereo_;

  private:
    std::mutex mutexPoints;
  };

} // namespace defSLAM

#endif // DEFORMATION MAPPOINT_H
