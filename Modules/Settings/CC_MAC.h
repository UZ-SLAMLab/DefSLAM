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
// Parameters for the ground truth stereo with matching
// by cross correlation. Set up depending of the dataset
#ifndef CC_MAC
#define CC_MAC
const int TEMPX(15);          // Parameters for the template to match
const int TEMPY(15);          // Parameters for the template
const int MARGIN(2);          // margin in the epipolar line. The images must be rectified
const int SEARCHX(300);       // Search window in columns
const double THRESHOLD(0.99); // Threshold of acceptance.
#endif
