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
#ifndef DEFLOCALMAPPING_H
#define DEFLOCALMAPPING_H

#include "KeyFrame.h"
#include "LocalMapping.h"
#include "WarpDatabase.h"

namespace ORB_SLAM2
{
  class KeyFrame;
}
namespace defSLAM
{
  using ORB_SLAM2::KeyFrame;
  using ORB_SLAM2::Map;
  using ORB_SLAM2::MapDrawer;

  class DefLocalMapping : public ORB_SLAM2::LocalMapping
  {
  public:
    /***********************************
   * Constructor of DefLocalMapping. It calls to the constructor to local mapping
   * and initializes the parameters of the deformable local mapping. It initializes 
   *  pointsToTemplate_: keypoints in the unexplored area to create a new template
   *  chiLimit_ : Limit in residual error of the Surface registration to accept that
   *            a template has been correctly aligned.
   *  reg: Weight for the Schwarzian regularizer in Schwarp estimation.
   *  bendingReg_ : Weight for the beding regularizer in Shape-from-normals.
   *********************/
    DefLocalMapping(Map *pMap, const string &strSettingPath);

    /***********************************
   * Destructor of DefLocalMapping. Just to remove the SchwarpDatabase.
   *********************/
    ~DefLocalMapping();

    /*********************************
   * Run. This function is to run in parallel the deformable 
   * mapping. It is an infinite loop that runs the deformable  
   * mapping in case of a new keyframe
   *********************************/
    virtual void Run() override;

    /*********************************
   * insideTheLoop(). This function runs the deformable 
   * mapping. It creates the 
   *********************************/
    void insideTheLoop();

    /*********************************
   * UpdateTemplate(). This function is called in the deformable 
   * tracking and updates the template with the reference keyframe 
   * selected. It creates the new map points with the surface and 
   * create the templates 
   *********************************/
    bool updateTemplate();

    // Stop the internal processes and reset the system.
    void ResetIfRequested() override;

    /*********************************
   * Create the new map points. They are extracted from the surface 
   * estimated for the keyframe with the Isometric NRSfM.
   ********************************/
    void CreateNewMapPoints() override;

  protected:
    /******
   * ProcessNewKeyframe. It add the keyframe to the warp database
   * to extract its warp with its covisible. In adition it add the
   * keyframe to a spanning covisibility tree as in ORBSLAM2:LocalMapping. 
   ******/
    void ProcessNewKeyFrame() override;
    /*************************
   * NRSfM. Given the warp and its derivates it estimates the normals.
   * With those normals it computes an up-to-scale surface. It align that
   * surface to recover the scale with the map points pose registered when
   * the keyframe was created
   *****************************/
    void NRSfM();

    /***************************************
   *  This function pick up the reference keyframe in case that there is no 
   * exploration. From the observed points we select the keyframe with the highest number
   * of observed map points.
  ****************************************/
    bool needNewTemplate();

    /***************************************
   *  This function pick up the reference keyframe in case that there is no 
   * exploration. From the observed points we select the keyframe with the highest number
   * of observed map points.
  ****************************************/
    KeyFrame *selectKeyframe();

  protected:
    // createTemplate tells when the Template must be changed
    bool createTemplate_;

    // DB of warps for the NRSfM. It contains the diferential
    // properties between keyframes to do the NRSfM.
    WarpDatabase *warpDB_;

    // Reference keyframe
    KeyFrame *referenceKF_;

    // keypoints in the unexplored area to create a new template
    int pointsToTemplate_;

    // Limit in residual error of the Surface registration to accept that
    // a template has been correctly aligned.
    double chiLimit_;

    // Weight for the beding regularizer in Shape-from-normals
    double bendingReg_;
    // Save results
    bool saveResults_;
  };

} // namespace defSLAM

#endif // LOCALMAPPING_H
