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

#ifndef DEFVIEWER_H
#define DEFVIEWER_H
#include "Viewer.h"
#include "FrameDrawer.h"
#include "MapDrawer.h"
#include "Tracking.h"
#include "System.h"
#include <mutex>

namespace ORB_SLAM2
{
    class Tracking;
    class FrameDrawer;
    class MapDrawer;
    class System;
    class Viewer;
} // namespace ORB_SLAM2

namespace defSLAM
{
    using ORB_SLAM2::FrameDrawer;
    using ORB_SLAM2::MapDrawer;
    using ORB_SLAM2::System;
    using ORB_SLAM2::Tracking;

    class DefViewer : public ORB_SLAM2::Viewer
    {
    public:
        // Constructor.
        DefViewer(System *pSystem, FrameDrawer *pFrameDrawer, MapDrawer *pMapDrawer, Tracking *pTracking, const string &strSettingPath);

        // Main thread function. Draw points, keyframes, the current camera pose and the last processed
        // frame. Drawing is refreshed according to the camera fps. We use Pangolin.
        void Run() override;

    private:
        double RegInex, RegLap, RegTemp; // Regularizers for deformable tracking.
    };

} // namespace defSLAM

#endif // VIEWER_H
