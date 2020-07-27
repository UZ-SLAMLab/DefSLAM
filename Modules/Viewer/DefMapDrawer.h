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

#ifndef DEFMAPDRAWER_H
#define DEFMAPDRAWER_H

#include "MapDrawer.h"
#include "MeshDrawer.h"

#include <pangolin/pangolin.h>

#include "Tracking.h"
#include <mutex>
#define MAX_NODES 300000
#define MAX_FACETS 100000

namespace ORB_SLAM2
{
  class MapDrawer;

}

namespace defSLAM
{
  class Node;
  class Edge;
  class Facet;
  class LaplacianMesh;
  using ORB_SLAM2::MapDrawer;

  class DefMapDrawer : public MapDrawer
  {
  public:
    // Constructor.
    DefMapDrawer(Map *, const string &);

    // Clear template history.
    void reset();

    // Draw the keyframes
    void DrawKeyFrames(const bool bDrawKF, const bool bDrawGraph) override;

    // Draw the template through meshdrawer.
    void DrawTemplate();

    // Draw the shape-at-rest of the template.
    void DrawTemplateAtRest(uint o);

    // Draw the templates through the execution.
    void DrawTemplatehist();

    // Updates from the tracking.
    void updateTemplate();

    // Update template at rest.
    void updateTemplateAtRest();

    // Draw template with texture.
    void ShowTexture(bool);

  protected:
    MeshDrawer *MeshDrawers;                                      // Mesh drawer of the template.
    std::map<KeyFrame *, unique_ptr<MeshDrawer>> MeshDrawershist; //Mesh drawers of the history of the template

  protected:
    bool Texture; // Flag for texture

    // Mutex
    std::mutex mTexture;
    std::mutex mTemplate;
    std::mutex mTemplatehist;
    std::mutex mPointsar;
  };

} // namespace defSLAM

#endif // DEFORMATION MAPDRAWER_H
