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

#include "LaplacianMesh.h"

namespace defSLAM
{
  using ORB_SLAM2::Map;
  using ORB_SLAM2::MapPoint;

  // Constructor from mesh.
  LaplacianMesh::LaplacianMesh(std::vector<std::vector<double>> &vertex,
                               std::vector<std::vector<int>> &index, Map *map)
      : TriangularMesh(vertex, index, map)
  {
    this->ExtractMeanCurvatures();
  }

  // Constructor from surface estimated for kF.
  LaplacianMesh::LaplacianMesh(std::set<MapPoint *> &mspMapPoints, Map *map,
                               KeyFrame *kF)
      : TriangularMesh(mspMapPoints, map, kF)
  {
    this->ExtractMeanCurvatures();
  }

  // Destructor
  LaplacianMesh::~LaplacianMesh() { LaplacianCoords.clear(); }

  /********
   * This function extract mean curvatures. It is thought for irregular meshes.
   * it use the method proposed in to extract Laplacian coordiantes that represent
   * the mean curvature.
   * ********/
  void LaplacianMesh::ExtractMeanCurvatures()
  {
    for (std::set<Node *>::iterator it = nodes_.begin();
         it != nodes_.end(); it++)
    {
      std::set<Node *> Neighbours = (*it)->GetNeighbours();
      Eigen::Vector3d Ni;
      Ni << (*it)->x, (*it)->y, (*it)->z;

      Eigen::Vector3d Laplacian;
      Laplacian.Zero();
      Eigen::Vector3d NeightMean;
      NeightMean.Zero();

      // Voronoi region of Ring of the node i
      for (std::set<Node *>::iterator ite = Neighbours.begin();
           ite != Neighbours.end(); ite++)
      {
        std::set<Node *> NeighboursofN = (*ite)->GetNeighbours();
        Eigen::Vector3d Nj;
        Nj << (*ite)->x, (*ite)->y, (*ite)->z;

        // Recover neighbours j+1 and j-1 of j
        std::vector<Node *> Nj1j_1;
        for (std::set<Node *>::iterator itn = NeighboursofN.begin();
             itn != NeighboursofN.end(); itn++)
        {
          if (Neighbours.count(*itn))
          {
            Nj1j_1.push_back(*itn);
          }
        }

        if (Nj1j_1.size() == 0)
        {
          (*ite)->setBadFlag();
        }
        else if (Nj1j_1.size() == 1)
        {
          (*ite)->setBoundary();
        }
        else
        {
          (*it)->NodesjJ_1J[*ite] =
              (std::pair<Node *, Node *>(Nj1j_1[0], Nj1j_1[1]));
          Eigen::Vector3d Nj1;
          Nj1 << (Nj1j_1[0])->x, (Nj1j_1[0])->y, (Nj1j_1[0])->z;

          Eigen::Vector3d Nj_1;
          Nj_1 << (Nj1j_1[1])->x, (Nj1j_1[1])->y, (Nj1j_1[1])->z;

          Eigen::Vector1d tnomega1;
          tnomega1 << ((Nj_1 - Ni).cross(Nj - Ni)).norm() /
                          ((Nj_1 - Ni).dot(Nj - Ni));
          Eigen::Vector1d tnomega2;
          tnomega2 << ((Nj1 - Ni).cross(Nj - Ni)).norm() /
                          ((Nj1 - Ni).dot(Nj - Ni));
          double wij = (tan(std::abs(atan(tnomega1(0))) / 2) +
                        tan(std::abs(atan(tnomega2(0))) / 2)) /
                       (Ni - Nj).norm();

          (*it)->weights[*ite] = (wij);
        }
      }
    }

    for (std::set<Node *>::iterator it = nodes_.begin();
         it != nodes_.end(); it++)
    {
      // Function not defined for boundaries.
      if (!(*it)->isBoundary())
      {
        std::set<Node *> Neighbours = (*it)->GetNeighbours();
        if (Neighbours.size() > 1)
        {
          Eigen::Vector3d Laplacian;
          Laplacian.Zero();
          Laplacian << 0, 0, 0;
          Eigen::Vector3d Ni;
          Ni << (*it)->x, (*it)->y, (*it)->z;
          double sumweights(0.0);
          std::set<Node *> Neighbours = (*it)->GetNeighbours();
          for (std::set<Node *>::iterator ite = Neighbours.begin();
               ite != Neighbours.end(); ite++)
          {
            Eigen::Vector3d Nj;
            Nj << (*ite)->x, (*ite)->y, (*ite)->z;
            Laplacian.noalias() = Laplacian + (*it)->weights[(*ite)] * (Nj);
            sumweights = sumweights + (*it)->weights[(*ite)];
          }

          LaplacianCoords[*it] = (Ni - (Laplacian / sumweights)).transpose();
        }
      }
    }
  }

  // Get the initial laplacian coordinate(vector of mean curvature) of the node
  Eigen::Vector3d const LaplacianMesh::GetLaplacianCoord(Node *n)
  {
    return LaplacianCoords.at(n);
  }

  // Get the initial mean curvature (norm of mean curvature)) of the node
  const Eigen::Vector1d LaplacianMesh::GetMeanCurvatureInitial(Node *n)
  {
    Eigen::Vector1d mc;
    mc << LaplacianCoords.at(n).norm();
    return mc;
  }

  // Get current mean curvature. It consider that the weights remain constant.
  Eigen::Vector1d const LaplacianMesh::GetMeanCurvature(Node *n)
  {
    Eigen::Vector3d Laplacian;
    Laplacian.Zero();
    Laplacian << 0, 0, 0;
    Eigen::Vector3d Ni;
    Ni << (n)->x, (n)->y, (n)->z;
    double sumweights(0.0);
    std::set<Node *> Neighbours = (n)->GetNeighbours();
    for (std::set<Node *>::iterator ite = Neighbours.begin();
         ite != Neighbours.end(); ite++)
    {
      Eigen::Vector3d Nj;
      Nj << (*ite)->x, (*ite)->y, (*ite)->z;

      Laplacian.noalias() = Laplacian + (n)->weights[(*ite)] * (Nj);
      sumweights =
          sumweights +
          (n)->weights[(
              *ite)];
    }
    Eigen::Vector1d mc;
    mc << (Laplacian.norm() / sumweights);

    return mc;
  }

  // Get current mean curvature as Laplacian vector
  Eigen::Vector3d const LaplacianMesh::GetMeanCurvatureVector(Node *n)
  {
    Eigen::Vector3d Laplacian;
    Laplacian.Zero();
    Laplacian << 0, 0, 0;
    Eigen::Vector3d Ni;
    Ni << (n)->x, (n)->y, (n)->z;
    double sumweights(0.0);
    std::set<Node *> Neighbours = (n)->GetNeighbours();
    for (std::set<Node *>::iterator ite = Neighbours.begin();
         ite != Neighbours.end(); ite++)
    {
      Eigen::Vector3d Nj;
      Nj << (*ite)->x, (*ite)->y, (*ite)->z;

      Laplacian.noalias() = Laplacian + (n)->weights[(*ite)] * (Nj);
      sumweights = sumweights + (n)->weights[(*ite)];
    }

    return Laplacian / sumweights;
  }

} // namespace defSLAM
