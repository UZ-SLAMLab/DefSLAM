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

#include "DefKeyFrame.h"
#include <TriangularMesh.h>
#include <algorithm>
#include <iterator>

namespace defSLAM
{
  using ORB_SLAM2::MapPoint;

  // Constructor from mesh given nodes and facets.
  TriangularMesh::TriangularMesh(std::vector<std::vector<double>> &vertexes,
                                 std::vector<std::vector<int>> &indexes, Map *map)
      : Template(map)
  {
    std::vector<Node *> NodesIndex;
    for (uint j = 0; j < vertexes.size() - 1; j++)
    {
      Node *NewNode =
          new Node(vertexes[j][0], vertexes[j][1], vertexes[j][2], j, this);
      nodes_.insert(NewNode);
      NodesIndex.push_back(NewNode);
      nodeArray_.push_back(NewNode);
      NewNode = nullptr;
      numVertices++;
    }

    for (uint j = 0; j < indexes.size(); j++)
    {
      facets_.insert(new Facet(NodesIndex[size_t(indexes[j][0])],
                               NodesIndex[size_t(indexes[j][1])],
                               NodesIndex[size_t(indexes[j][2])], this));
    }
  }

  // Constructor for a mesh given a surface estimated in the kf.
  TriangularMesh::TriangularMesh(std::set<MapPoint *> &mspMapPoints, Map *map,
                                 KeyFrame *kf)
      : Template(mspMapPoints, map, kf)
  {
    Surface *RefSurface = static_cast<DefKeyFrame *>(kf)->surface;
    std::vector<cv::Mat> NodesSurface;
    int nodVer = 10;
    int nodHor = 10;
    RefSurface->getVertex(NodesSurface, nodVer, nodHor);
    auto facets = this->regularTriangulation(nodVer, nodHor);

    std::vector<std::vector<float>> vertexW;
    std::vector<std::vector<float>> vertexC;

    cv::Mat Twc = kf->GetPoseInverse();
    vertexW.reserve(NodesSurface.size());

    for (uint i(0); i < NodesSurface.size(); i++)
    {
      std::vector<float> ptW;
      ptW.reserve(3);
      cv::Mat x3wh(4, 1, CV_32F);
      x3wh = Twc * NodesSurface[i];
      ptW.push_back(x3wh.at<float>(0, 0));
      ptW.push_back(x3wh.at<float>(1, 0));
      ptW.push_back(x3wh.at<float>(2, 0));
      vertexW.push_back(ptW);
    }
    this->setNodes(vertexW, facets);
    this->calculateFeaturesCoordinates();
    this->texture = kf->RGBimage.clone();
    this->getFacetTexture(kf);
  }

  // Create a regular triangular mesh from the keyframe surface estimated.
  std::vector<std::vector<int>> TriangularMesh::regularTriangulation(int nodesVer, int nodesHor)
  {
    std::vector<std::vector<int>> facets;
    for (int j(0); j < nodesHor - 1; j++)
    {
      for (int i(0); i < nodesVer - 1; i++)
      {
        std::vector<int> facet1 = {i + nodesHor * j, i + nodesHor * j + 1, (nodesHor * (j + 1)) + i};
        std::vector<int> facet2 = {i + nodesHor * j + 1, (nodesHor * (j + 1)) + i, (nodesHor * (j + 1)) + i + 1};
        facets.push_back(facet1);
        facets.push_back(facet2);
      }
    }

    return facets;
  }

  // Create the nodes and the facets to define the triagular mesh.
  void TriangularMesh::setNodes(std::vector<std::vector<float>> &vertex,
                                std::vector<std::vector<int>> &facets)
  {
    std::vector<Node *> NodesIndex;
    NodesIndex.reserve(vertex.size());
    for (uint j = 0; j < vertex.size(); j++)
    {
      Node *NewNode = new Node(vertex[j][0], vertex[j][1], vertex[j][2], j, this);
      nodes_.insert(NewNode);
      NodesIndex.push_back(NewNode);
      NewNode = nullptr;
      numVertices++;
    }
    facets.reserve(facets.size());
    for (uint j = 0; j < facets.size(); j++)
    {
      Facet *facet = new Facet(NodesIndex[size_t(facets[j][0])],
                               NodesIndex[size_t(facets[j][1])],
                               NodesIndex[size_t(facets[j][2])], this);
      facets_.insert(facet);
    }
  }
  // Embed the points inside the facets by calculating its barycentric coordinates.
  void TriangularMesh::calculateFeaturesCoordinates()
  {
    std::for_each(
        mapPoints_.begin(), mapPoints_.end(),
        [this](MapPoint *const &it) {
          if ((it))
          {
            if (!it->isBad())
            {

              (*it).lastincorporasion = false;
              // Barycentric coordinates calculated with each KeyFrame
              Eigen::Vector3f MapPoint;

              cv::Mat MapPointPosition = (*it).GetWorldPos();

              MapPoint << MapPointPosition.at<float>(0),
                  MapPointPosition.at<float>(1), MapPointPosition.at<float>(2);

              Node *ClosestNode = static_cast<Node *>(nullptr);
              double bestdist(100);
              for (std::set<Node *>::iterator itNd = nodes_.begin();
                   itNd != nodes_.end(); itNd++)
              {
                double x, y, z;
                (*itNd)->getXYZ(x, y, z);
                double dist =
                    std::sqrt(pow(x - MapPoint(0), 2) + pow(y - MapPoint(1), 2) +
                              pow(z - MapPoint(2), 2));
                if (dist < bestdist)
                {
                  ClosestNode = (*itNd);
                  bestdist = dist;
                }
              }
              if (ClosestNode)
              {
                std::set<Facet *> FacetsNode = ClosestNode->GetFacets();
                for (std::set<Facet *>::iterator ita = FacetsNode.begin();
                     ita != FacetsNode.end(); ita++)
                {
                  std::set<Node *> node = (*ita)->getNodes();
                  std::vector<Node *> nods;
                  std::vector<Eigen::Vector3f> nodsEigen;
                  for (std::set<Node *>::iterator ite = node.begin();
                       ite != node.end(); ite++)
                  {
                    nods.push_back(*ite);
                    Eigen::Vector3f node;
                    node << (*ite)->x, (*ite)->y, (*ite)->z;
                    nodsEigen.push_back(node);
                  }
                  Eigen::Vector3f barycentric;
                  if (pointInTriangle(MapPoint, nodsEigen[0], nodsEigen[1],
                                      nodsEigen[2], barycentric))
                  {
                    static_cast<DefMapPoint *>(it)->SetCoordinates(
                        barycentric(0), barycentric(1), barycentric(2));
                    static_cast<DefMapPoint *>(it)->SetFacet(*ita);
                    static_cast<DefMapPoint *>(it)->Repose();
                    break;
                  }
                }
              }
            }
          }
        });
  }

  /******
   *  Estimate barycentric coordinates. It returns true if the point is inside the face.
   * Inspired in code presented in
   * https://math.stackexchange.com/questions/544946/determine-if-projection-of-3d-point-onto-plane-is-within-a-triangle
   * *******/
  bool TriangularMesh::pointInTriangle(const Eigen::Vector3f &query_point,
                                       const Eigen::Vector3f &triangle_vertex_0,
                                       const Eigen::Vector3f &triangle_vertex_1,
                                       const Eigen::Vector3f &triangle_vertex_2,
                                       Eigen::Vector3f &barycentric)
  {
    // u=P2−P1
    Eigen::Vector3f u = triangle_vertex_1 - triangle_vertex_0;
    // v=P3−P1
    Eigen::Vector3f v = triangle_vertex_2 - triangle_vertex_0;
    // n=u×v
    Eigen::Vector3f n = u.cross(v);
    // w=P−P1
    Eigen::Vector3f w = query_point - triangle_vertex_0;
    // Barycentric coordinates of the projection P′of P onto T:
    // γ=[(u×w)⋅n]/n²
    float gamma = u.cross(w).dot(n) / n.dot(n);
    // β=[(w×v)⋅n]/n²
    float beta = w.cross(v).dot(n) / n.dot(n);
    float alpha = 1 - gamma - beta;
    barycentric << alpha, beta, gamma;
    Eigen::Vector3f newPose = triangle_vertex_0 * alpha +
                              triangle_vertex_1 * beta +
                              triangle_vertex_2 * gamma;
    if ((newPose - query_point).squaredNorm() > 1E-1)
      return false;
    // The point P′ lies inside T if:
    return ((0 <= alpha) && (alpha <= 1) && (0 <= beta) && (beta <= 1) &&
            (0 <= gamma) && (gamma <= 1));
  }

  // Extract texture for each facet.
  void TriangularMesh::getFacetTexture(KeyFrame *KF)
  {
    ////// Loop the faces
    for (std::set<Facet *>::iterator itf = facets_.begin();
         itf != facets_.end(); itf++)
    {
      /// Extract texture of the image;
      (*itf)->getTextureCoordinates(KF);
    }
  }
} // namespace defSLAM
