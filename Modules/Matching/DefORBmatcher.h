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

#ifndef DEFORBMATCHER_H
#define DEFORBMATCHER_H
#pragma once
#include <ORBmatcher.h>
#include <Thirdparty/BBS/bbs_MAC.h>
#include <vector>

namespace defSLAM
{
    using ORB_SLAM2::Frame;
    using ORB_SLAM2::KeyFrame;
    using ORB_SLAM2::ORBmatcher;

    class DefORBmatcher : public ORBmatcher
    {
    public:
        // Constructor.
        DefORBmatcher(float nnratio = 0.6, bool checkOri = true);

        // Find more matches using a Schwarp. x is the array with the nodes of the warp.
        void findbyWarp(
            KeyFrame *Kf1, KeyFrame *Kf2,
            vector<pair<size_t, size_t>> &vMatchedIndices,
            double (&x)[_NumberOfControlPointsU * _NumberOfControlPointsV * 2],
            double lambda);

        /*********************
        * Initialize a warp and make a refinement with the Schwarzian equation
        * to search for map points. HERE IS INITIALIZED THE WARP. IF YOU DO NOT
        * CALL findbyWarp IN SchwarpDatabase.cc, YOU MUST INITIALIZE THE WARP
        * BEFORE OPTIMIZING THE SCHWARP. (To initialize a warp, use the function
        * Warps::Warp::initialize() from Mapping/Schwarp.h)
        **********************/
        void CalculateInitialSchwarp(
            KeyFrame *Kf1i, KeyFrame *Kf2i,
            vector<pair<size_t, size_t>> &vMatchedIndices,
            double (&x)[_NumberOfControlPointsU * _NumberOfControlPointsV * 2],
            double lambda);

        // It use the result of the schwarzian warp to search for map points.
        int searchBySchwarp(
            KeyFrame *pKF1, KeyFrame *pKF2,
            double (&x)[_NumberOfControlPointsU * _NumberOfControlPointsV * 2],
            std::vector<pair<size_t, size_t>> &vMatchedPairs);

        // Copy of function SearchByProjection of ORBmatcher.h. It is used in
        // in the tracking. It only get points embedded in the 3D template
        int SearchByProjection(Frame &CurrentFrame, const Frame &LastFrame,
                               const float th, const bool bMono);
    };

} // namespace defSLAM

#endif // ORBMATCHER_H
