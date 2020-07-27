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

#include "DefOptimizer.h"
#include "LaplacianMesh.h"

#include "Converter.h"
#include "DefMap.h"
#include <Eigen/StdVector>
#include <fstream>
#include <iostream>
#include <mutex>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "Thirdparty/g2o/g2o/core/block_solver.h"
#include "Thirdparty/g2o/g2o/core/optimization_algorithm_levenberg.h"
#include "Thirdparty/g2o/g2o/core/robust_kernel_impl.h"
#include "Thirdparty/g2o/g2o/solvers/linear_solver_dense.h"
#include "Thirdparty/g2o/g2o/solvers/linear_solver_eigen.h"
#include "Thirdparty/g2o/g2o/types/sft_types.h"
#include "Thirdparty/g2o/g2o/types/types_seven_dof_expmap.h"
#include "Thirdparty/g2o/g2o/types/types_six_dof_expmap.h"
namespace defSLAM
{
  using ORB_SLAM2::Converter;

  namespace Optimizer
  {
    /******
     * Static functions for the optimizer.
     * ***/

    // Set nodes to optimize. This function numbers them from 1 to N+1 and
    // save the result in indx for each Node.
    void setMeshNodes(g2o::SparseOptimizer &optimizer, Map *mMap,
                      double fixed = false);

    // Updates node position.
    void updateNodes(g2o::SparseOptimizer &optimizer, Map *mMap);

    // Copy of pose Optimization from optimizer in ORBSLAM to increase thresholds
    // for deformable scenarios.
    int poseOptimization(Frame *pFrame)
    {
      g2o::SparseOptimizer optimizer;
      g2o::BlockSolver_6_3::LinearSolverType *linearSolver;

      linearSolver =
          new g2o::LinearSolverDense<g2o::BlockSolver_6_3::PoseMatrixType>();

      g2o::BlockSolver_6_3 *solver_ptr = new g2o::BlockSolver_6_3(linearSolver);

      g2o::OptimizationAlgorithmLevenberg *solver =
          new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
      optimizer.setAlgorithm(solver);

      int nInitialCorrespondences = 0;

      // Set Frame vertex
      g2o::VertexSE3Expmap *vSE3 = new g2o::VertexSE3Expmap();
      vSE3->setEstimate(Converter::toSE3Quat(pFrame->mTcw));
      vSE3->setId(0);
      vSE3->setFixed(false);
      optimizer.addVertex(vSE3);

      // Set MapPoint vertices
      const int N = pFrame->N;

      vector<g2o::EdgeSE3ProjectXYZOnlyPose *> vpEdgesMono;
      vector<size_t> vnIndexEdgeMono;
      vpEdgesMono.reserve(N);
      vnIndexEdgeMono.reserve(N);

      vector<g2o::EdgeStereoSE3ProjectXYZOnlyPose *> vpEdgesStereo;
      vector<size_t> vnIndexEdgeStereo;
      vpEdgesStereo.reserve(N);
      vnIndexEdgeStereo.reserve(N);

      const float deltaMono = sqrt(5.991);

      {
        unique_lock<mutex> lock(MapPoint::mGlobalMutex);

        for (int i = 0; i < N; i++)
        {
          MapPoint *pMP = pFrame->mvpMapPoints[i];
          if (pMP)
          {
            if (static_cast<DefMap *>(pMP->mpMap)->GetTemplate())
            {
              if (!static_cast<DefMapPoint *>(pMP)->getFacet())
              {
                continue;
              }

              nInitialCorrespondences++;
              // pFrame->mvbOutlier[i] = false;

              Eigen::Matrix<double, 2, 1> obs;
              const cv::KeyPoint &kpUn = pFrame->mvKeysUn[i];
              obs << kpUn.pt.x, kpUn.pt.y;

              g2o::EdgeSE3ProjectXYZOnlyPose *e =
                  new g2o::EdgeSE3ProjectXYZOnlyPose();

              e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex *>(
                                  optimizer.vertex(0)));
              e->setMeasurement(obs);
              const float invSigma2 = pFrame->mvInvLevelSigma2[kpUn.octave];
              e->setInformation(Eigen::Matrix2d::Identity() * invSigma2);

              g2o::RobustKernelHuber *rk = new g2o::RobustKernelHuber;
              e->setRobustKernel(rk);
              rk->setDelta(deltaMono);

              e->fx = pFrame->fx;
              e->fy = pFrame->fy;
              e->cx = pFrame->cx;
              e->cy = pFrame->cy;
              cv::Mat Xw = pMP->GetWorldPos();
              e->Xw[0] = Xw.at<float>(0);
              e->Xw[1] = Xw.at<float>(1);
              e->Xw[2] = Xw.at<float>(2);

              optimizer.addEdge(e);

              vpEdgesMono.push_back(e);
              vnIndexEdgeMono.push_back(i);
              e->computeError();
            }
          }
        }
      }

      if (nInitialCorrespondences < 3)
        return 0;

      // We perform 4 optimizations, after each optimization we classify observation
      // as inlier/outlier
      // At the next optimization, outliers are not included, but at the end they
      // can be classified as inliers again.
      const float chi2Mono[5] = {5.991, 5.991, 5.991, 5.991, 5.991};
      const int its[4] = {10, 10, 10, 10};

      int nBad = 0;
      for (uint ai = 0; ai < 2; ai++)
      {
        for (size_t it = 0; it < 4; it++)
        {
          vSE3->setEstimate(Converter::toSE3Quat(pFrame->mTcw));
          optimizer.initializeOptimization(ai);
          optimizer.optimize(its[it]);
          nBad = 0;
          if (ai == 0)
          {
            for (size_t i = 0, iend = vpEdgesMono.size(); i < iend; i++)
            {
              g2o::EdgeSE3ProjectXYZOnlyPose *e = vpEdgesMono[i];

              const size_t idx = vnIndexEdgeMono[i];

              if (pFrame->mvbOutlier[idx])
              {
                e->computeError();
              }

              const float chi2 = e->chi2();

              if (chi2 > chi2Mono[it])
              {
                pFrame->mvbOutlier[idx] = true;
                e->setLevel(1);
                nBad++;
              }
              else
              {
                pFrame->mvbOutlier[idx] = false;
                e->setLevel(0);
              }
            }
          }
          else
          {
            for (size_t i = 0, iend = vpEdgesMono.size(); i < iend; i++)
            {
              g2o::EdgeSE3ProjectXYZOnlyPose *e = vpEdgesMono[i];

              const size_t idx = vnIndexEdgeMono[i];

              if (pFrame->mvbOutlier[idx])
              {
                e->computeError();
              }

              const float chi2 = e->chi2();

              if (chi2 > chi2Mono[it])
              {
                pFrame->mvbOutlier[idx] = true;
                e->setLevel(2);
              }
              else
              {
                pFrame->mvbOutlier[idx] = false;
                nBad--;
              }
            }
          }
          if (optimizer.edges().size() < 10)
            break;
        }
      }

      vSE3->setEstimate(Converter::toSE3Quat(pFrame->mTcw));
      optimizer.initializeOptimization(0);
      optimizer.optimize(10);

      // Recover optimized pose and return number of inliers
      g2o::VertexSE3Expmap *vSE3_recov =
          static_cast<g2o::VertexSE3Expmap *>(optimizer.vertex(0));
      g2o::SE3Quat SE3quat_recov = vSE3_recov->estimate();
      cv::Mat pose = Converter::toCvMat(SE3quat_recov);
      pFrame->SetPose(pose);

      return nInitialCorrespondences - nBad;
    }

    // Shape-from-template with camera motion estimation. This function is used in
    // the deformation tracking to estimate the map deformation and camera pose each
    // frame.
    int DefPoseOptimization(Frame *pFrame, Map *mMap, double RegLap, double RegInex,
                            double RegTemp, uint NeighboursLayers)
    {
      g2o::SparseOptimizer optimizer;
      g2o::BlockSolverX::LinearSolverType *linearSolver;
      linearSolver =
          new g2o::LinearSolverDense<g2o::BlockSolverX::PoseMatrixType>();
      g2o::BlockSolverX *solver_ptr = new g2o::BlockSolverX(linearSolver);
      // std::cout << "OPTIMIZER " << PropagationZone << std::endl;
      g2o::OptimizationAlgorithmLevenberg *solver =
          new g2o::OptimizationAlgorithmLevenberg(solver_ptr);

      optimizer.setAlgorithm(solver);

      int nInitialCorrespondences = 0;

      // Set Frame vertex
      g2o::VertexSE3Expmap *vSE3 = new g2o::VertexSE3Expmap();
      vSE3->setEstimate(Converter::toSE3Quat(pFrame->mTcw));
      vSE3->setId(0);
      vSE3->setFixed(false);
      optimizer.addVertex(vSE3);
      setMeshNodes(optimizer, mMap);
      uint nBad(0);
      // Set MapPoint vertices
      const int N = pFrame->N;
      const int Ncorr = pFrame->mvpMapPointsCorr.size();

      std::set<Node *> ViewedNodes;
      vector<g2o::EdgeNodesCamera *> vpEdgesMono, vpEdgesMonocorr;
      vector<size_t> vnIndexEdgeMono, vnIndexEdgeMonocorr;
      vpEdgesMono.reserve(N);
      vnIndexEdgeMono.reserve(N);
      vnIndexEdgeMonocorr.reserve(Ncorr);
      /////////// OBSERVATIONS ////////////////////////
      const float deltaMono = sqrt(5.991);
      unique_lock<mutex> lock(MapPoint::mGlobalMutex);

      int Points_Obs(0);
      int Points_Facet(0);
      int PointsORB(0);

      for (int i = 0; i < N; i++)
      {
        if (!(pFrame->mvbOutlier[i]))
        {
          MapPoint *pMP = pFrame->mvpMapPoints[i];
          if ((pMP))
          {
            Points_Obs++;
            if ((pMP)->isBad())
              continue;
            if (static_cast<DefMapPoint *>(pMP)->getFacet())
            {
              // Monocular observation
              nInitialCorrespondences++;
              pFrame->mvbOutlier[i] = false;

              Eigen::Matrix<double, 2, 1> obs;
              const cv::KeyPoint &kpUn = pFrame->mvKeysUn[i];
              obs << kpUn.pt.x, kpUn.pt.y;

              g2o::EdgeNodesCamera *e = new g2o::EdgeNodesCamera();

              std::set<Node *> NodesMp =
                  static_cast<DefMapPoint *>(pMP)->getFacet()->getNodes();
              // std::cout << NodesMp.size() << std::endl;
              e->resize(NodesMp.size() + 1);
              e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex *>(
                                  optimizer.vertex(0)));
              uint ui(1);
              Eigen::Vector3d bary;
              bary << static_cast<DefMapPoint *>(pMP)->b1,
                  static_cast<DefMapPoint *>(pMP)->b2,
                  static_cast<DefMapPoint *>(pMP)->b3;
              for (std::set<Node *>::iterator it = NodesMp.begin();
                   it != NodesMp.end(); it++)
              {
                e->setVertex(ui, dynamic_cast<g2o::OptimizableGraph::Vertex *>(
                                     optimizer.vertex((*it)->getIndex())));
                ViewedNodes.insert(*it);
                (*it)->setViewed();

                ui++;
              }
              e->setBarycentric(bary);
              e->setMeasurement(obs);

              const float invSigma2 = pFrame->mvInvLevelSigma2[kpUn.octave];
              e->setInformation(Eigen::Matrix2d::Identity() * invSigma2 / N);

              g2o::RobustKernelHuber *rk = new g2o::RobustKernelHuber;
              e->setRobustKernel(rk);
              rk->setDelta(deltaMono);

              e->fx = pFrame->fx;
              e->fy = pFrame->fy;
              e->cx = pFrame->cx;
              e->cy = pFrame->cy;
              e->computeError();

              Points_Facet++;
              optimizer.addEdge(e);

              vpEdgesMono.push_back(e);
              vnIndexEdgeMono.push_back(i);
              PointsORB++;
            }
          }
        }
      }

      //// Temporal smoothness
      double m(1.0);
      if (static_cast<DefMap *>(mMap)->GetTemplate())
        m = (static_cast<DefMap *>(mMap)->GetTemplate())->getEdgeMeanSize();
      for (std::set<Node *>::iterator it = ViewedNodes.begin();
           it != ViewedNodes.end(); it++)
      {

        g2o::EdgesReference *e = new g2o::EdgesReference();
        e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex *>(
                            optimizer.vertex((*it)->getIndex())));
        double x, y, z;
        Eigen::Vector3d v;
        (*it)->getInitialPose(x, y, z);
        v << x, y, z;
        e->setMeasurement(v);
        e->setInformation(RegTemp * Eigen::Matrix3d::Identity() / pow(m, 2));
        optimizer.addEdge(e);
        e->computeError();
      }

      std::set<Node *> OptLap;
      std::set<Node *> NewOptNodes;
      OptLap = ViewedNodes;
      int a(0);
      for (uint n(0); n < (NeighboursLayers); n++)
      {
        for (std::set<Node *>::iterator it = ViewedNodes.begin();
             it != ViewedNodes.end(); it++)
        {
          std::set<Node *> Neighbours = (*it)->GetNeighbours();
          for (std::set<Node *>::iterator ite = Neighbours.begin();
               ite != Neighbours.end(); ite++)
          {
            NewOptNodes.insert(*ite);
            a++;
          }
        }
        for (std::set<Node *>::iterator it = NewOptNodes.begin();
             it != NewOptNodes.end(); it++)
        {
          OptLap.insert(*it);
        }
      }

      /// Laplacian Regularizer
      double meanlap(0.0);
      /// Implementation ECCV
      for (std::set<Node *>::iterator it = OptLap.begin(); it != OptLap.end();
           it++)
      {
        (dynamic_cast<g2o::OptimizableGraph::Vertex *>(
             optimizer.vertex((*it)->getIndex())))
            ->setFixed(false);
        (dynamic_cast<g2o::OptimizableGraph::Vertex *>(
             optimizer.vertex((*it)->getIndex())))
            ->setMarginalized(false);

        if (!(*it)->isBoundary())
        {
          std::set<Node *> Neighbours = (*it)->GetNeighbours();
          std::set<Edge *> edges = (*it)->getEdges();

          uint Index_Neigh(0);
          for (std::set<Edge *>::iterator ite = edges.begin();
               ite != edges.end(); ite++)
          {
            g2o::EdgeMeanCurvature *e = new g2o::EdgeMeanCurvature;
            e->resize(Neighbours.size() + 1);
            (*it)->setLocal();
            e->SetNeighbourgEdge(Index_Neigh);
            e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex *>(
                                optimizer.vertex((*it)->getIndex())));
            uint j(0);
            std::vector<double> weights;
            for (std::set<Node *>::iterator ite2 = Neighbours.begin();
                 ite2 != Neighbours.end(); ite2++)
            {
              e->setVertex(j + 1, dynamic_cast<g2o::OptimizableGraph::Vertex *>(
                                      optimizer.vertex((*ite2)->getIndex())));
              weights.push_back((*it)->weights[*ite2]);
              j++;
            }
            e->setDistanceEdges((*ite)->getDist());
            e->setWeights(weights);
            Eigen::Vector1D InitialMeanCurvature;
            InitialMeanCurvature =
                static_cast<LaplacianMesh *>(
                    static_cast<DefMap *>(mMap)->GetTemplate())
                    ->GetMeanCurvatureInitial(*it);
            e->setMeasurement(InitialMeanCurvature);
            e->computeError();
            auto a = e->errorData();
            double er = std::abs(a[0]);
            meanlap += er / OptLap.size();
            e->setInformation(RegLap * Eigen::Vector1D::Identity() / OptLap.size());
            optimizer.addEdge(e);
            Index_Neigh++;
          }
        }
      }

      // Inextensibility
      std::vector<double> err_inex;

      std::set<Edge *> medges;
      for (std::set<Node *>::iterator itn = OptLap.begin(); itn != OptLap.end();
           itn++)
      {
        std::set<Edge *> medgesnode = (*itn)->getEdges();
        for (std::set<Edge *>::iterator ite = medgesnode.begin();
             ite != medgesnode.end(); ite++)
        {
          medges.insert(*ite);
        }
      }
      //   std::cout << "Cauchy Strain :" << std::endl;
      double meaninex(0.0);
      std::vector<g2o::EdgesStreching *> vpInextensibilityEdges;
      for (std::set<Edge *>::iterator ite = medges.begin(); ite != medges.end();
           ite++)
      {
        g2o::EdgesStreching *e = new g2o::EdgesStreching;
        // e->resize(2);
        uint is(0);
        std::set<Node *> medgesnode = (*ite)->getNodes();
        for (std::set<Node *>::iterator itn = medgesnode.begin();
             itn != medgesnode.end(); itn++)
        {
          e->setVertex(is, dynamic_cast<g2o::OptimizableGraph::Vertex *>(
                               optimizer.vertex((*itn)->getIndex())));
          is++;
        }
        Eigen::Vector1D s;
        s << (*ite)->getDist();
        e->setMeasurement(s);
        e->setInformation(RegInex * Eigen::Vector1D::Identity() / medges.size());
        e->computeError();
        vpInextensibilityEdges.push_back(e);
        auto a = e->errorData();
        double er = sqrt(pow(a[0], 2));
        meaninex += (er / medges.size());
        err_inex.push_back(er);
        optimizer.addEdge(e);
      }

      optimizer.setVerbose(false);

      optimizer.initializeOptimization(0);

      optimizer.optimize(50);

      for (size_t i = 0, iend = vpEdgesMono.size(); i < iend; i++)
      {
        g2o::EdgeNodesCamera *e = vpEdgesMono[i];

        const size_t idx = vnIndexEdgeMono[i];

        if (pFrame->mvbOutlier[idx])
        {
          e->computeError();
        }

        const float chi2 = e->chi2();

        if (chi2 > 5.991)
        {
          pFrame->mvbOutlier[idx] = true;
          nBad++;
        }
        else
        {
          pFrame->mvbOutlier[idx] = false;
        }
      }
      double sumError(0.0);
      std::vector<float> vectorError;
      uint n(0);
      for (size_t i = 0, iend = vpEdgesMono.size(); i < iend; i++)
      {
        g2o::EdgeNodesCamera *e = vpEdgesMono[i];

        const size_t idx = vnIndexEdgeMono[i];

        if (!pFrame->mvbOutlier[idx])
        {
          e->computeError();
          auto a = e->errorData();
          double er = sqrt(pow(a[0], 2) + pow(a[1], 2));
          vectorError.push_back(er);
          sumError += er;
          n++;
        }
      }
      cout << "Reprojection error: " << sumError / n << endl;
      cout << "Points considered: " << n << endl;
      pFrame->repError = sumError / n;

      // updateNodes(optimizer,mMap);
      g2o::VertexSE3Expmap *vSE3_recov =
          static_cast<g2o::VertexSE3Expmap *>(optimizer.vertex(0));
      g2o::SE3Quat SE3quat_recov = vSE3_recov->estimate();
      cv::Mat pose = Converter::toCvMat(SE3quat_recov);
      pFrame->SetPose(pose);

      std::vector<MapPoint *> mPoints = mMap->GetAllMapPoints();

      updateNodes(optimizer, mMap);
      for (std::vector<MapPoint *>::iterator pMP = mPoints.begin();
           pMP != mPoints.end(); pMP++)
      {
        if (static_cast<DefMapPoint *>(*pMP)->getFacet())
          static_cast<DefMapPoint *>(*pMP)->RecalculatePosition();
      }
      return nInitialCorrespondences - nBad;
    }

    // Shape-from-template w/o camera motion. This function was used
    // for testing shape-from-template with conventional dataset.
    int DefPoseOptimization(const std::vector<std::vector<double>> &matches,
                            Frame *pFrame, Map *mMap, std::vector<bool> &outlier,
                            double RegLap, double RegInex, double RegTemp)
    {

      g2o::SparseOptimizer optimizer;
      g2o::BlockSolverX::LinearSolverType *linearSolver;
      linearSolver =
          new g2o::LinearSolverDense<g2o::BlockSolverX::PoseMatrixType>();
      g2o::BlockSolverX *solver_ptr = new g2o::BlockSolverX(linearSolver);

      g2o::OptimizationAlgorithmLevenberg *solver =
          new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
      // g2o::OptimizationAlgorithmDogleg* solver = new
      // g2o::OptimizationAlgorithmDogleg(solver_ptr);

      optimizer.setAlgorithm(solver);

      int nInitialCorrespondences = 0;

      // Set Frame vertex
      g2o::VertexSE3Expmap *vSE3 = new g2o::VertexSE3Expmap();
      vSE3->setEstimate(Converter::toSE3Quat(pFrame->mTcw));

      vSE3->setId(0);
      vSE3->setFixed(false);
      optimizer.addVertex(vSE3);

      setMeshNodes(optimizer, mMap, 0);

      // Set MapPoint vertices
      std::set<Node *> ViewedNodes;
      vector<g2o::EdgeNodesCamera *> vpEdgesMono;

      std::vector<Node *> Nodes =
          static_cast<DefMap *>(mMap)->GetTemplate()->nodeArray_;
      for (std::vector<Node *>::iterator it = Nodes.begin(); it != Nodes.end();
           it++)
      {
        optimizer.vertex((*it)->getIndex())->setFixed(false);
        optimizer.vertex((*it)->getIndex())->setMarginalized(false);
      }

      /////////// OBSERVATIONS ////////////////////////
      const float deltaMono = 0.5; // sqrt(5.991);
      unique_lock<mutex> lock(MapPoint::mGlobalMutex);
      for (uint i = 0; i < matches.size(); i++)
      {
        // Monocular observation
        nInitialCorrespondences++;

        Eigen::Matrix<double, 2, 1> obs;
        obs << matches[i][6], matches[i][7];
        g2o::EdgeNodesCamera *e = new g2o::EdgeNodesCamera();
        e->resize(4);
        e->setVertex(
            0, dynamic_cast<g2o::OptimizableGraph::Vertex *>(optimizer.vertex(0)));
        Eigen::Vector3d bary;
        bary << matches[i][3], matches[i][4], matches[i][5];

        for (int ui = 1; ui < 4; ui++)
        {
          e->setVertex(
              ui, dynamic_cast<g2o::OptimizableGraph::Vertex *>(
                      optimizer.vertex(Nodes[matches[i][ui - 1]]->getIndex())));
          static_cast<g2o::VertexSBAPointXYZ *>(
              optimizer.vertex(Nodes[matches[i][ui - 1]]->getIndex()))
              ->setMarginalized(false);
          static_cast<g2o::VertexSBAPointXYZ *>(
              optimizer.vertex(Nodes[matches[i][ui - 1]]->getIndex()))
              ->setFixed(false);
          ViewedNodes.insert(Nodes[matches[i][ui - 1]]);
          (Nodes[matches[i][ui - 1]])->setViewed();
        }

        e->setBarycentric(bary);
        e->setMeasurement(obs);

        e->setInformation(Eigen::Matrix2d::Identity() / double(matches.size()));

        g2o::RobustKernelHuber *rk = new g2o::RobustKernelHuber;
        e->setRobustKernel(rk);
        rk->setDelta(deltaMono);
        e->fx = pFrame->fx;
        e->fy = pFrame->fy;
        e->cx = pFrame->cx;
        e->cy = pFrame->cy;

        e->computeError();
        vpEdgesMono.push_back(e);
        optimizer.addEdge(e);
      }

      //// Temporal smoothness. It is regularized with the mean distance of the edges.
      double m(1.0);
      if (static_cast<DefMap *>(mMap)->GetTemplate())
        m = (static_cast<DefMap *>(mMap)->GetTemplate())->getEdgeMeanSize();

      for (std::set<Node *>::iterator it = ViewedNodes.begin();
           it != ViewedNodes.end(); it++)
      {
        g2o::EdgesReference *e = new g2o::EdgesReference();
        e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex *>(
                            optimizer.vertex((*it)->getIndex())));
        double x, y, z;
        Eigen::Vector3d v;
        (*it)->getXYZ(x, y, z);
        v << x, y, z;
        e->setMeasurement(v);
        e->setInformation(RegTemp * Eigen::Matrix3d::Identity() / pow(m, 2));
        // optimizer.addEdge(e);
        e->computeError();
      }

      ///////////
      //////////// LOCAL+AFFECTED DEFORMATION
      std::set<Node *> OptLap = ViewedNodes;

      /// Laplacian Regularizer
      //  std::cout << "Mean curvature :" << std::endl;
      double meanlap(0.0);
      // g2o::EdgeMeanCurvature* ejem = new g2o::EdgeMeanCurvature();

      /// Implementation ECCV
      for (std::set<Node *>::iterator it = OptLap.begin(); it != OptLap.end();
           it++)
      {
        if (!(*it)->isBoundary())
        {
          std::set<Node *> Neighbours = (*it)->GetNeighbours();
          std::map<Node *, uint> DictionaryNodeNumber;
          uint ind(0);
          for (std::set<Node *>::iterator ite = Neighbours.begin();
               ite != Neighbours.end(); ite++)
          {
            DictionaryNodeNumber[*ite] = ind;
            ind++;
          }
          uint Index_Neigh(0);
          for (std::set<Node *>::iterator ite = Neighbours.begin();
               ite != Neighbours.end(); ite++)
          {
            g2o::EdgeMeanCurvature *e = new g2o::EdgeMeanCurvature;
            e->resize(Neighbours.size() + 1);
            (*it)->setLocal();
            e->SetNeighbourgEdge(Index_Neigh);
            e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex *>(
                                optimizer.vertex((*it)->getIndex())));
            uint j(0);
            std::vector<double> weights;
            std::vector<std::pair<uint, uint>> NeightWeights;

            for (std::set<Node *>::iterator ite2 = Neighbours.begin();
                 ite2 != Neighbours.end(); ite2++)
            {
              e->setVertex(j + 1, dynamic_cast<g2o::OptimizableGraph::Vertex *>(
                                      optimizer.vertex((*ite2)->getIndex())));
              (dynamic_cast<g2o::OptimizableGraph::Vertex *>(
                   optimizer.vertex((*it)->getIndex())))
                  ->setFixed(false);
              (dynamic_cast<g2o::OptimizableGraph::Vertex *>(
                   optimizer.vertex((*it)->getIndex())))
                  ->setMarginalized(false);
              NeightWeights.push_back(std::pair<uint, uint>(
                  DictionaryNodeNumber[(*it)->NodesjJ_1J[(*ite2)].first],
                  DictionaryNodeNumber[(*it)->NodesjJ_1J[(*ite2)].second]));
              weights.push_back((*it)->weights[*ite2]);
              j++;
            }
            e->setWeights(weights);
            Eigen::Vector1D InitialMeanCurvature;
            InitialMeanCurvature =
                static_cast<LaplacianMesh *>(
                    static_cast<DefMap *>(mMap)->GetTemplate())
                    ->GetMeanCurvatureInitial(*it);
            //   e->setBarycentric(static_cast<LaplacianMesh*>(static_cast<DefMap*>(mMap)->GetTemplate())->GetLaplacianCoord(*it));
            e->setMeasurement(InitialMeanCurvature);
            e->computeError();
            auto a = e->errorData();
            double er = std::abs(a[0]); //+pow(a[1],2)+pow(a[2],2));
            meanlap += er / OptLap.size();
            e->setInformation(RegLap * Eigen::Vector1D::Identity() / OptLap.size());
            // ejem=e;

            optimizer.addEdge(e);
            Index_Neigh++;
          }
        }
      }

      // Inextensibility
      std::vector<double> err_inex;

      // std::cout << "Cauchy Strain :" << std::endl;
      double meaninex(0.0);
      std::set<Edge *> medges =
          static_cast<LaplacianMesh *>(
              static_cast<DefMap *>(mMap)->GetTemplate())
              ->getEdges();
      std::cout << medges.size() << std::endl;

      std::vector<g2o::EdgesStreching *> vpInextensibilityEdges;
      for (std::set<Edge *>::iterator ite = medges.begin(); ite != medges.end();
           ite++)
      {
        g2o::EdgesStreching *e = new g2o::EdgesStreching;
        // e->resize(2);
        uint is(0);
        std::set<Node *> medgesnode = (*ite)->getNodes();
        for (std::set<Node *>::iterator itn = medgesnode.begin();
             itn != medgesnode.end(); itn++)
        {
          e->setVertex(is, dynamic_cast<g2o::OptimizableGraph::Vertex *>(
                               optimizer.vertex((*itn)->getIndex())));
          is++;
        }
        Eigen::Vector1D s;
        s << (*ite)->getDist();
        e->setMeasurement(s);
        e->setInformation(RegInex * Eigen::Vector1D::Identity() / medges.size());
        e->computeError();
        vpInextensibilityEdges.push_back(e);
        auto a = e->errorData();
        double er = sqrt(pow(a[0], 2));
        meaninex += (er / medges.size());
        err_inex.push_back(er);
        optimizer.addEdge(e);
      }
      // std::cout << "RegInex : " << RegInex << " " <<    meaninex<< std::endl;

      /// Optimization Block
      /// ///
      optimizer.setVerbose(false);
      optimizer.initializeOptimization(0);
      optimizer.optimize(50);
      uint i(0);
      for (auto edgemono : vpEdgesMono)
      {
        auto a = edgemono->errorData();
        if (deltaMono < sqrt(pow(a[0], 2) + pow(a[1], 2)))
        {
          outlier.push_back(true);
        }
        else
        {
          outlier.push_back(false);
        }
        i++;
      }
      /////////////////// PROPAGATION OF THE SOLUTION
      ///

      updateNodes(optimizer, mMap);

      return 0;
    }

    // Function to Sim(3) alignment of two point clouds. Used in surface registration.
    bool OptimizeHorn(std::vector<std::vector<float>> &pts1,
                      std::vector<std::vector<float>> &pts2, g2o::Sim3 &g2oS12,
                      double chi, double huber)
    {
      g2o::SparseOptimizer optimizer;
      g2o::BlockSolverX::LinearSolverType *linearSolver;

      linearSolver =
          new g2o::LinearSolverDense<g2o::BlockSolverX::PoseMatrixType>();

      g2o::BlockSolverX *solver_ptr = new g2o::BlockSolverX(linearSolver);

      g2o::OptimizationAlgorithmLevenberg *solver =
          new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
      optimizer.setAlgorithm(solver);
      g2o::Sim3 g2oasdS12(g2oS12);
      g2o::VertexSim3ExpmapNoProj *vert0 = new g2o::VertexSim3ExpmapNoProj;
      vert0->setEstimate(g2oasdS12);
      vert0->setId(0);

      optimizer.addVertex(vert0);
      // Set MapPoint vertices
      const int N = pts1.size();
      // const vector<MapPoint*> vpMapPoints1 = pKF1->GetMapPointMatches();

      const float deltaHuber = sqrt(huber);

      std::vector<g2o::EdgeSim3Simple *> edgesSimple;
      edgesSimple.reserve(N);
      for (int i = 0; i < N; i++)
      {
        g2o::EdgeSim3Simple *simple = new g2o::EdgeSim3Simple;
        simple->setVertex(
            0, dynamic_cast<g2o::OptimizableGraph::Vertex *>(optimizer.vertex(0)));

        Eigen::Matrix<double, 3, 1> v1;
        v1 << pts1[i][0], pts1[i][1], pts1[i][2];

        Eigen::Matrix<double, 3, 1> v2;
        v2 << pts2[i][0], pts2[i][1], pts2[i][2];

        simple->setPoints(v1, v2);
        simple->setInformation(Eigen::Matrix3d::Identity());

        g2o::RobustKernelHuber *rk = new g2o::RobustKernelHuber;
        simple->setRobustKernel(rk);
        rk->setDelta(deltaHuber);
        optimizer.addEdge(simple);
        edgesSimple.push_back(simple);
        simple->setLevel(0);
      }
      // Optimize!
      optimizer.setVerbose(false);
      optimizer.initializeOptimization(0);
      optimizer.optimize(50);

      // Check inliers
      g2oS12 = vert0->estimate();

      int count(0);
      for (int i = 0; i < N; i++)
      {
        if (edgesSimple[i]->chi2() > chi)
        {
          // edgesSimple[i]->setLevel(1);
        }
        else
        {
          count++;
        }
      }
      optimizer.initializeOptimization(0);
      optimizer.optimize(50);

      std::cout << "chi " << optimizer.chi2() / count << "  " << chi << "  "
                << count << std::endl;
      if (std::isnan(optimizer.chi2()))
        return false;
      else if (std::isinf(optimizer.chi2()))
        return false;
      else
        return (optimizer.chi2() / count < chi);
    }

    // Set nodes to optimize. This function numbers them from 1 to N+1 and
    // save the result in indx for each Node.
    void setMeshNodes(g2o::SparseOptimizer &optimizer, Map *mMap, double prop)
    {
      // Add Nodes Vertex
      // Add Nodes from 2 to N nodes
      // The two first index are for the temporal constrains
      std::set<Node *> Nodes =
          static_cast<DefMap *>(mMap)->GetTemplate()->getNodes();
      uint j(1);
      if (prop)
        j = 0;

      for (std::set<Node *>::iterator it = Nodes.begin(); it != Nodes.end(); it++)
      {
        g2o::VertexSBAPointXYZ *vNode = new g2o::VertexSBAPointXYZ();
        Eigen::Vector3d v;
        double x, y, z;
        (*it)->getXYZ(x, y, z);
        v << x, y, z;
        vNode->setEstimate(v);
        vNode->setFixed(true);
        vNode->setMarginalized(true);
        vNode->setId(j);
        (*it)->setIndex(j);
        optimizer.addVertex(vNode);
        j++;
      };
    }

    // Updates node position.
    void updateNodes(g2o::SparseOptimizer &optimizer, Map *mMap)
    {
      std::set<Node *> Nodes =
          static_cast<DefMap *>(mMap)->GetTemplate()->getNodes();
      for (std::set<Node *>::iterator it = Nodes.begin(); it != Nodes.end(); it++)
      {
        (*it)->update();
        (*it)->resetRole();
        Eigen::Vector3d v = static_cast<const g2o::VertexSBAPointXYZ *>(
                                optimizer.vertex((*it)->getIndex()))
                                ->estimate();
        (*it)->setXYZ(v(0), v(1), v(2));
      }
    }

  } // namespace Optimizer
} // namespace defSLAM
