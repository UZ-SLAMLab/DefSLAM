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

#ifndef TEMPLATE_H
#define TEMPLATE_H

#include "opencv2/core.hpp"
#include <Edge.h>
#include <Facet.h>
#include <MapPoint.h>
#include <Node.h>
#include <Surface.h>
#include <set>

namespace defSLAM
{
  class MapPoint;
  class Node;
  class Edge;
  class Facet;

  class Template
  {
  public:
    /*********
   *  Constructor for the Template. It is unique for each keyframe. In this point
   * it could be a template with any representation. The only requirements is that
   * it is going to embed map points.
   * *******/
    Template(std::set<MapPoint *> &mspMapPoints, Map *map, KeyFrame *kf);

    // Constructor just with the map.
    Template(Map *map);

    // Destructor.
    virtual ~Template();

    // Add map point to the template
    void addMapPoint(MapPoint *pMP);

    // Erase map point from the template
    void eraseMapPoint(MapPoint *);

    // Add node to the template
    void addNode(Node *pN);

    // Remove node from the template
    void eraseNode(Node *);

    // Add facet to the template
    void addFacet(Facet *);

    // Remove facet to the template
    void eraseFacet(Facet *);

    // Add edge to the template
    void addEdge(Edge *);

    // Remove edge from the template
    void EraseEdge(Edge *);

    // Get median size of the edge of the template
    double getEdgeMeanSize();

    // Get map points embedded in the template
    const std::set<MapPoint *> getPoints();

    // Get nodes.
    const std::set<Node *> getNodes();

    // Get facets.
    const std::set<Facet *> getFacets();

    // Get edges.
    const std::set<Edge *> getEdges();

    // Restart template
    virtual void restart();

    // Get texture. This texture is given by an image.
    cv::Mat getTexture();

  public:
    std::vector<Node *> nodeArray_; // Array of nodes of the template

    std::mutex mutexDrawer_; //

    KeyFrame *kf; // Keyframe of the template

  protected:
    Map *mMap; // Map.

    unsigned int numVertices;

    std::set<MapPoint *> mapPoints_;

    std::set<Node *> nodes_;

    std::set<Edge *> edges_;

    std::set<Facet *> facets_;

    std::mutex mMutexMap;

    std::mutex mMutexNodes;

    std::mutex mMutexEdges;

    std::mutex mMutexFaces;

    std::mutex MutexUsingMap;

    cv::Mat texture;

    bool UsingMap;

    bool Updating;
  };
} // namespace defSLAM
#endif
