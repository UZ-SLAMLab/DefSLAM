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
*
* g2o_types defined for the defSLAM
*
**************************************/

#include <fstream>
#include <iostream>

// IF it fails, try including this

#include <../g2o/g2o/core/base_binary_edge.h>
#include <../g2o/g2o/core/base_unary_edge.h>
#include <../g2o/g2o/core/base_vertex.h>
#include <../g2o/g2o/core/eigen_types.h>

#include <../g2o/g2o/core/base_multi_edge.h>

namespace Eigen
{
  typedef Matrix<double, 1, 1, ColMajor> Vector1D;
}
namespace g2o
{
  using namespace std;

  class VertexDouble : public g2o::BaseVertex<1, Eigen::Vector1D>
  {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    VertexDouble() {}

    virtual bool read(std::istream & /*is*/)
    {
      cerr << __PRETTY_FUNCTION__ << " not implemented yet" << endl;
      return false;
    }

    virtual bool write(std::ostream & /*os*/) const
    {
      cerr << __PRETTY_FUNCTION__ << " not implemented yet" << endl;
      return false;
    }

    virtual void setToOriginImpl()
    {
      cerr << __PRETTY_FUNCTION__ << " not implemented yet" << endl;
    }

    virtual void oplusImpl(const double *update)
    {
      Eigen::Vector1D::ConstMapType v(update);
      _estimate.noalias() += v;
    }
  };

  class EdgeNodesCamera : public BaseMultiEdge<2, Eigen::Vector2d>
  {
    /****
     * Class to optimize the reprojection error using 
     * barycentric coordinates
     * *****/
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    EdgeNodesCamera() = default;

    virtual bool read(std::istream &is)
    {
      std::cout << "Not implemented" << std::endl;
      return false;
    }

    virtual bool write(std::ostream &os) const
    {
      std::cout << "Not implemented" << std::endl;
      return false;
    }

    void setBarycentric(Eigen::Vector3d const &Barycentrics)
    {
      Barycentric << Barycentrics;
    }

    void computeError()
    {
      const VertexSE3Expmap *Camera =
          static_cast<const VertexSE3Expmap *>(_vertices[0]);
      const VertexSBAPointXYZ *v1 =
          static_cast<const VertexSBAPointXYZ *>(_vertices[1]);
      const VertexSBAPointXYZ *v2 =
          static_cast<const VertexSBAPointXYZ *>(_vertices[2]);
      const VertexSBAPointXYZ *v3 =
          static_cast<const VertexSBAPointXYZ *>(_vertices[3]);

      // Eigen::Matrix3d vABC;
      auto Point3Dw = Barycentric(0) * v1->estimate() +
                      Barycentric(1) * v2->estimate() +
                      Barycentric(2) * v3->estimate();

      auto Point3Dc = Camera->estimate().map(Point3Dw);
      Vector2D Estimation;
      Estimation = this->cam_project_1(Point3Dc);
      _error << (_measurement - Estimation);
    }

    Vector2D cam_project_1(const Eigen::Vector3d &trans_xyz)
    {
      Vector2D res;
      Vector2D res1;
      res1(0) = trans_xyz(0) / trans_xyz(2);
      res1(1) = trans_xyz(1) / trans_xyz(2);
      res[0] = res1[0] * fx + cx;
      res[1] = res1[1] * fy + cy;
      return res;
    }

    double fx, fy, cx, cy;

    void linearizeOplus()
    {
      const VertexSE3Expmap *Camera =
          static_cast<const VertexSE3Expmap *>(_vertices[0]);
      const VertexSBAPointXYZ *v1 =
          static_cast<const VertexSBAPointXYZ *>(_vertices[1]);
      const VertexSBAPointXYZ *v2 =
          static_cast<const VertexSBAPointXYZ *>(_vertices[2]);
      const VertexSBAPointXYZ *v3 =
          static_cast<const VertexSBAPointXYZ *>(_vertices[3]);
      SE3Quat T(Camera->estimate());

      Eigen::Vector3d xyz1 = Camera->estimate().map(v1->estimate());
      Eigen::Vector3d xyz2 = Camera->estimate().map(v2->estimate());
      Eigen::Vector3d xyz3 = Camera->estimate().map(v3->estimate());

      Eigen::Vector3d xyz =
          xyz1 * Barycentric(0) + xyz2 * Barycentric(1) + xyz3 * Barycentric(2);

      auto Trot = T.rotation().toRotationMatrix();
      double x = xyz(0);
      double y = xyz(1);
      double z = xyz(2);
      double z_2 = z * z;

      _jacobianOplus[0](0, 0) = x * y / z_2 * fx;
      _jacobianOplus[0](0, 1) = -(1 + (x * x / z_2)) * fx;
      _jacobianOplus[0](0, 2) = y / z * fx;
      _jacobianOplus[0](0, 3) = -1. / z * fx;
      _jacobianOplus[0](0, 4) = 0;
      _jacobianOplus[0](0, 5) = x / z_2 * fx;

      _jacobianOplus[0](1, 0) = (1 + y * y / z_2) * fy;
      _jacobianOplus[0](1, 1) = -x * y / z_2 * fy;
      _jacobianOplus[0](1, 2) = -x / z * fy;
      _jacobianOplus[0](1, 3) = 0;
      _jacobianOplus[0](1, 4) = -1. / z * fy;
      _jacobianOplus[0](1, 5) = y / z_2 * fy;

      x = xyz1(0);
      y = xyz1(1);
      z = xyz1(2);

      Eigen::Matrix<double, 2, 3> tmp;
      tmp(0, 0) = fx;
      tmp(0, 1) = 0;
      tmp(0, 2) = -x / z * fx;

      tmp(1, 0) = 0;
      tmp(1, 1) = fy;
      tmp(1, 2) = -y / z * fy;

      _jacobianOplus[1] = -1. / z * tmp * Trot * Barycentric(0);

      x = xyz2(0);
      y = xyz2(1);
      z = xyz2(2);

      tmp(0, 2) = -x / z * fx;
      tmp(1, 2) = -y / z * fy;

      _jacobianOplus[2] = -1. / z * tmp * Trot * Barycentric(1);
      x = xyz3(0);
      y = xyz3(1);
      z = xyz3(2);
      tmp(0, 2) = -x / z * fx;
      tmp(1, 2) = -y / z * fy;

      _jacobianOplus[3] = -1. / z * tmp * Trot * Barycentric(2);
    }

  private:
    Eigen::Vector3d Barycentric;
  };

  class EdgeMeanCurvature : public BaseMultiEdge<1, Eigen::Vector1D>
  {
    /****
     * Class to optimize the changes in mean curvature. It
     * use fix weights pre-estimate wrt. the neightbours.
     * It use the class LaplacianMesh to estimate the initial 
     * curvature.
     * *****/
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    EdgeMeanCurvature() = default;

    virtual bool read(std::istream &is)
    {
      std::cout << "Not implemented" << std::endl;
      return false;
    }

    virtual bool write(std::ostream &os) const
    {
      std::cout << "Not implemented" << std::endl;
      return false;
    }

    void setBarycentric(Eigen::Vector3d const &Barycentrics)
    {
      initial_coordinate = Barycentrics;
    }

    void setWeights(std::vector<double> a) { weights = a; }

    void setDistanceEdges(double lenghtEdge) { lenghtEdge_ = lenghtEdge; }

    void SetNeighbourgEdge(uint indx) { Index_Number = indx; }

    void get_Neight()
    {
      for (uint i(1); i < _vertices.size(); i++)
      {
        const VertexSBAPointXYZ *v1 =
            static_cast<const VertexSBAPointXYZ *>(_vertices[i]);
        Vertex.push_back(v1->estimate());
      }
    }

    void computeError()
    {
      const VertexSBAPointXYZ *v =
          static_cast<const VertexSBAPointXYZ *>(_vertices[0]);
      std::vector<Eigen::Vector3d> Vertex2;
      Eigen::Vector3d Ni;
      Ni = v->estimate();

      for (uint i(1); i < _vertices.size(); i++)
      {
        const VertexSBAPointXYZ *v1 =
            static_cast<const VertexSBAPointXYZ *>(_vertices[i]);
        Vertex2.push_back(v1->estimate());
      }

      /// Calculate Zero Curvature Position
      Eigen::Vector3d Neightweight;
      Neightweight << 0, 0, 0;
      double sumweight_aux(0.0);

      for (uint i(0); i < Vertex2.size(); i++)
      {
        Neightweight.noalias() = Neightweight + weights[i] * (Vertex2[i]);
        sumweight_aux = sumweight_aux + weights[i];
      }
      sumWeights_ = sumweight_aux;
      Eigen::Vector3d ZeroCurv;
      ZeroCurv << Neightweight / sumWeights_;

      meanCurvature_ << 0, 0, 0;
      meanCurvature_ = Ni - ZeroCurv;
      meanCurvatureNorm_ = meanCurvature_.norm();
      _error << (meanCurvatureNorm_ - _measurement(0)) /
                    lenghtEdge_;
    }

    void linearizeOplus()
    {
      _jacobianOplus[0] = meanCurvature_.transpose();

      for (uint i(0); i < _jacobianOplus.size(); i++)
      {
        if (meanCurvatureNorm_ < 1E-15)
        {
          _jacobianOplus[i] << 0.0, 0.0, 0.0;
          continue;
        }
        if (i > 0)
        {
          double weight_aux = -(weights[i - 1] / sumWeights_);
          _jacobianOplus[i] = weight_aux * meanCurvature_.transpose();
        }
        _jacobianOplus[i] = _jacobianOplus[i] / (meanCurvatureNorm_ * lenghtEdge_);
      }
    }

    std::vector<JacobianType, Eigen::aligned_allocator<JacobianType>>
    getJacobian()
    {
      return _jacobianOplus;
    }

  private:
    Eigen::Vector3d initial_coordinate;
    std::vector<double> weights;
    uint Index_Number;
    std::vector<Eigen::Vector3d> Vertex;
    double lenghtEdge_;
    double sumWeights_;
    Eigen::Vector3d meanCurvature_;
    double meanCurvatureNorm_;
  };

  class EdgesStreching
      : public BaseBinaryEdge<1, Eigen::Vector1D, VertexSBAPointXYZ,
                              VertexSBAPointXYZ>
  {
    /****
     * Class to optimize the changes in length for each edge.
     * *****/
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    EdgesStreching() = default;

    virtual bool read(std::istream &is)
    {
      std::cout << "Not implemented" << std::endl;
      return false;
    }

    virtual bool write(std::ostream &os) const
    {
      std::cout << "Not implemented" << std::endl;
      return false;
    }

    void computeError()
    {
      const VertexSBAPointXYZ *v1 =
          static_cast<const VertexSBAPointXYZ *>(_vertices[0]);
      const VertexSBAPointXYZ *v2 =
          static_cast<const VertexSBAPointXYZ *>(_vertices[1]);
      _error =
          ((v1->estimate() - v2->estimate()).norm()) * _measurement.inverse() -
          Eigen::MatrixXd::Ones(1, 1);

      //  _error << (_error-_measurement);
    }
    void linearizeOplus()
    {
      const VertexSBAPointXYZ *v1 =
          static_cast<const VertexSBAPointXYZ *>(_vertices[0]);
      const VertexSBAPointXYZ *v2 =
          static_cast<const VertexSBAPointXYZ *>(_vertices[1]);
      Eigen::MatrixXd ddo =
          ((v1->estimate() - v2->estimate()).norm() * _measurement).inverse();

      Eigen::MatrixXd a = (v1->estimate() - v2->estimate()) * ddo;

      _jacobianOplusXi << a(0), a(1), a(2);
      _jacobianOplusXj << -a(0), -a(1), -a(2);
    }
  };

  class EdgesReference : public BaseUnaryEdge<3, Eigen::Vector3d, VertexSBAPointXYZ>
  {
    /*************************
     * It penalizes the different between the initial and the final 
     * position, decoupling the movements of the mesh and the camera.
     * ********************/
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    EdgesReference() = default;

    virtual bool read(std::istream &is)
    {
      std::cout << "Not implemented" << std::endl;
      return false;
    }

    virtual bool write(std::ostream &os) const
    {
      std::cout << "Not implemented" << std::endl;
      return false;
    }

    void computeError()
    {
      const VertexSBAPointXYZ *v1 =
          static_cast<const VertexSBAPointXYZ *>(vertex(0));
      _error << (v1->estimate() - _measurement);
    }

    void linearizeOplus() { _jacobianOplusXi.setIdentity(); }
  };
} // namespace g2o
