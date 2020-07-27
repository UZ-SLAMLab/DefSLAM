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

#include "DefMap.h"
#include "TemplateGenerator.h"

#include <mutex>

namespace defSLAM
{
  /***********
   * Constructor:
   *  Mesh is the template (template is a keyword of c++)
   *  flag there is a template 
   * 
   *********/
  DefMap::DefMap()
      : Map(), ThereIsATemplate(false), Mesh(static_cast<Template *>(nullptr)) {}
  /***********
   * Constructor:
   *  It creates a map with a template with the vertex and facets pointed.
   *    std::vector<std::vector<double>> &vertex
   *    std::vector<std::vector<int>> &index
   *  flag there is a template 
   * 
   *********/
  DefMap::DefMap(std::vector<std::vector<double>> &vertex,
                 std::vector<std::vector<int>> &index)
      : Map(), Mesh(static_cast<Template *>(nullptr))
  {
    Mesh = TemplateGenerator::ChargeLaplacianMesh(vertex, index, this);
  }
  /***********
   *  Function that generates a template with the estimated surface for the keyframe kf. 
   * (Kf->surface must have been initializated.)
   *********/
  void DefMap::createTemplate(KeyFrame *Kf)
  {
    // It creates a template with the current mapped points
    Mesh = TemplateGenerator::LaplacianMeshCreate(mspMapPoints, this, Kf);

    ThereIsATemplate = true;
  }

  // return the template
  Template *DefMap::GetTemplate() { return Mesh; }

  // deletes the current template
  void DefMap::clearTemplate()
  {
    ThereIsATemplate = false;
    if (Mesh)
    {
      delete Mesh;
      Mesh = static_cast<Template *>(nullptr);
    }
    for (set<MapPoint *>::iterator sit = mspMapPoints.begin(),
                                   send = mspMapPoints.end();
         sit != send; sit++)
    {
      static_cast<DefMapPoint *>(*sit)->SetFacet(
          static_cast<Facet *>(nullptr));
    }
  }
  // clear the entire map
  void DefMap::clear()
  {
    this->clearTemplate();
    for (set<MapPoint *>::iterator sit = mspMapPoints.begin(),
                                   send = mspMapPoints.end();
         sit != send; sit++)
      delete *sit;

    for (set<KeyFrame *>::iterator sit = mspKeyFrames.begin(),
                                   send = mspKeyFrames.end();
         sit != send; sit++)
      delete *sit;

    mspMapPoints.clear();
    mspKeyFrames.clear();
    mnMaxKFid = 0;
    mvpReferenceMapPoints.clear();
    mvpKeyFrameOrigins.clear();
  }

  // Return the reference keyframe or the keyframe used to
  // create the template
  KeyFrame *DefMap::getLastKFTemplate()
  {
    std::unique_lock<std::mutex> lck(MutexKf);
    return lastKfTemplate_;
  }

  // Set the reference keyframe or the keyframe used to
  // create the template
  void DefMap::setLastKFTemplate(KeyFrame *kf)
  {
    std::unique_lock<std::mutex> lck(MutexKf);
    lastKfTemplate_ = kf;
  }
} // namespace defSLAM
