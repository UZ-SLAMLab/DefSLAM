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

#include <Template.h>

#include <DefMapPoint.h>

namespace defSLAM
{
  /*********
   *  Constructor for the Template. It is unique for each keyframe. In this point
   * it could be a template with any representation. The only requirements is that
   * it is going to embed map points.
   * *******/
  Template::Template(std::set<MapPoint *> &mspMapPoints, Map *map, KeyFrame *kf)
      : kf(kf), mMap(map), numVertices(0)
  {
    for (std::set<MapPoint *>::iterator it = mspMapPoints.begin();
         it != mspMapPoints.end(); it++)
    {
      if (*it)
      {
        if (!(*it)->isBad())
        {
          this->addMapPoint(*it);
        }
      }
    }
    edges_.clear();
  }

  // Constructor for the map.
  Template::Template(Map *map) : mMap(map) {}

  // Destructor.
  Template::~Template()
  {
    std::unique_lock<std::mutex> lc(mutexDrawer_);
    texture.release();
    mMap = nullptr;
    for (std::set<MapPoint *>::iterator itmp = mapPoints_.begin();
         itmp != mapPoints_.end(); itmp++)
    {
      if (*itmp)
        if (!(*itmp)->isBad())
          static_cast<DefMapPoint *>(*itmp)->RemoveTemplate();
    }
    mapPoints_.clear();
    {
      std::unique_lock<std::mutex> a(mMutexEdges);

      for (std::set<Edge *>::iterator ite = edges_.begin();
           ite != edges_.end(); ite++)
      {
        delete (*ite);
      }
      edges_.clear();
    }
    {
      std::unique_lock<std::mutex> a(mMutexFaces);

      for (std::set<Facet *>::iterator itf = facets_.begin();
           itf != facets_.end(); itf++)
      {
        delete (*itf);
      }
      facets_.clear();
    }

    {
      std::unique_lock<std::mutex> a(mMutexNodes);
      for (std::set<Node *>::iterator itn = nodes_.begin();
           itn != nodes_.end(); itn++)
      {
        delete (*itn);
      }
      nodes_.clear();
    }

    Facet::TexturesDataset.clear();
  }

  // Add map point to the template
  void Template::addMapPoint(MapPoint *pMP)
  {
    unique_lock<mutex> lock(mMutexMap);
    mapPoints_.insert(pMP);
  }

  // Erase map point from the template
  void Template::eraseMapPoint(MapPoint *pMP)
  {
    unique_lock<mutex> lock(mMutexMap);
    mapPoints_.erase(pMP);
  }

  // Add node to the template
  void Template::addNode(Node *pN)
  {
    unique_lock<mutex> lock(mMutexNodes);
    nodes_.insert(pN);
  }

  // Remove node from the template
  void Template::eraseNode(Node *pN)
  {
    unique_lock<mutex> lock(mMutexNodes);
    nodes_.erase(pN);
  }

  // Add facet to the template
  void Template::addFacet(Facet *pN)
  {
    unique_lock<mutex> lock(mMutexFaces);
    facets_.insert(pN);
  }

  // Remove facet to the template
  void Template::eraseFacet(Facet *pN)
  {
    unique_lock<mutex> lock(mMutexFaces);
    facets_.erase(pN);
  }

  // Add edge to the template
  void Template::addEdge(Edge *pN)
  {
    unique_lock<mutex> lock(mMutexEdges);
    edges_.insert(pN);
  }

  // Remove edge from the template
  void Template::EraseEdge(Edge *pN)
  {
    unique_lock<mutex> lock(mMutexEdges);
    edges_.erase(pN);
  }

  // Get median size of the edge of the template
  double Template::getEdgeMeanSize()
  {
    unique_lock<mutex> lock(mMutexEdges);

    vector<double> dists;
    for (std::set<Edge *>::iterator ite = edges_.begin();
         ite != edges_.end(); ite++)
    {
      dists.push_back((*ite)->getDist());
    }
    std::sort(dists.begin(), dists.end());
    if (dists.size() > 0)
      return dists[dists.size() / 2];
    else
    {
      return 0.10;
    }
  }

  // Get map points embedded in the template
  const std::set<MapPoint *> Template::getPoints()
  {
    return this->mapPoints_;
  }

  // Get nodes.
  const std::set<Node *> Template::getNodes() { return this->nodes_; }

  // Get edges.
  const std::set<Edge *> Template::getEdges() { return this->edges_; }

  // Get facets.
  const std::set<Facet *> Template::getFacets() { return this->facets_; }

  // It takes the template to the shape-at-rest, i.e. its original form.
  void Template::restart()
  {
    for (std::set<Node *>::iterator it = nodes_.begin();
         it != nodes_.end(); it++)
    {
      (*it)->reset();
    }
  }

  // Get texture. This texture is given by an image.
  cv::Mat Template::getTexture() { return this->texture.clone(); }
} // namespace defSLAM
