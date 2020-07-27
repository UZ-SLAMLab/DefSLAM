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

#ifndef MESHDRAWER_H
#define MESHDRAWER_H

#include <opencv2/core.hpp>

#include <iostream>
#include <vector>
#include <set>
#include <mutex>

namespace defSLAM
{
    class MeshDrawer
    {
        /*************
         * Class to draw meshes
         ************/
    public:
        // Constructor
        MeshDrawer() = default;

        // Destructor
        ~MeshDrawer() = default;

        /******************
         * In case the mesh was created from a keyframe. The 
         * image can be used as texture. 
         * If there is no texture fit it with a random image
         * (3 Channels).
         *************/
        void addTextureImage(cv::Mat im);

        /******************
        * It registers the edges of the mesh <e_{i},e_{j}> to draw them.
        *************/
        void addEdge(std::vector<int> edge);

        /******************
         * Add a node to draw and its projection. The position as 
         * a vector of doubles with 3 components {x,y,z} and the 
         * projection as a vector of 2 components. The role defines 
         * the color of the node. 
          *************/
        void addNode(std::vector<double> pos, std::vector<double> proj, int role);

        /******************
        * Add a facet as a vector of 3 integers 
        *************/
        void addFacet(std::vector<int> fac);

        /****************
         *  Function to draw the mesh.
         *      alpha: level of transparency
         *      drawedges: flag to draw the edges.
         *************/
        void drawMesh(double alpha = 1, bool drawedges = true);

    private:
        std::vector<int> Facets;    // One vector [V1_1 V1_2 V1_3 ..Vn_1 Vn_2 Vn_3 ... VN_1 VN_2 VN_3] Facets
        std::vector<int> Edges;     // One vector [V1_1 V1_2 ... Vn_1 Vn_2 ... VN_1 VN_2] Edges
        std::vector<double> Nodes;  // One vector [V1 V2 ... Vn ... VN] Vertex position
        std::vector<int> NodesRole; // One vector [R1 R2 ... Rn ... RN] Nodes role

        cv::Mat Texture;
        std::vector<float> NodesProjection; // One vector with the 2d projection of the vertex [v1 v2 ... vn ... vN]

    private:
        std::mutex mtemp;
    };

} // namespace defSLAM
#endif // MESHS_H
