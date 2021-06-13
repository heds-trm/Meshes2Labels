/*
* This file is part of the m2l application (https://github.com/heds-trm/Meshes2Labels).
* Copyright (c) 2012-2021 Jerome Schmid.
* Geneva School of Health Sciences (HEdS) 
* University of Applied Sciences and Arts Western Switzerland (HES-SO)
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, version 3.
*
* This program is distributed in the hope that it will be useful, but
* WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
* General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <map>
#include <limits>
#include <tclap/CmdLine.h>
#include "vtkPolyDataReader.h"
#include "vtkPolyDataWriter.h"
#include "vtkPolyData.h"
#include "vtkStructuredPoints.h"
#include "vtkCell.h"
#include "vtkImageSeedConnectivity.h"
#include "vtkStructuredPointsReader.h"
#include "vtkStructuredPointsWriter.h"
#include "vtkDataSetWriter.h"
#include "vtkUnstructuredGridReader.h"
#include "vtkSTLReader.h"
#include "vtkUnstructuredGrid.h"
#include <itkVTKImageIO.h>
#include "vtkPolyDataNormals.h"
#include "vtkMatrix4x4.h"
#include "vtkMatrixToHomogeneousTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkTriangle.h"
#include "vtkCommand.h"
#include <vtkOBBTree.h>
#include <vtkSmartPointer.h>
#include <vtkAppendPolyData.h>
#include <itkVTKImageToImageFilter.h>
#include <itkImageToVTKImageFilter.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkImageStencil.h>
#include <vtkImageToStructuredPoints.h>

#include <itkChangeInformationImageFilter.h>
#include "itkVector.h"
#include "itkOffset.h"
#include "itkIndex.h"
#include "gdiam.h"
#include "io.h"
#include <vector>

// global vars ugly but so useful
double origin[3];
double spacing[3];
int size[3];
double directionx[3];
double directiony[3];
double directionz[3];
vtkStructuredPoints* oImage = nullptr;
bool bUchar = false;
double padding;
double glob_model2obbox_M[16];
double glob_obbox2model_M[16];
bool bForceUseOfOBBTree = false;
bool bTestMeshBboxInImage = true;
bool bUseVtkStencil = true;

template <typename T>
T norm(T u[3]) { return (T)sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]); }


void SafeDelete(vtkObjectBase* obj)
{
  if (obj) obj->Delete();
  obj = nullptr;
}


template <typename T>
void normalizeVec(T a[3])
{
  T n = norm<T>(a);
  if (n > 0.00000001)
  {
    a[0] /= n;
    a[1] /= n;
    a[2] /= n;
  }
}

template <typename T>
void Identity(T M[16]) { M[0] = 1; M[1] = 0; M[2] = 0; M[3] = 0; M[4] = 0; M[5] = 1; M[6] = 0; M[7] = 0; M[8] = 0; M[9] = 0; M[10] = 1; M[11] = 0; M[12] = 0; M[13] = 0; M[14] = 0; M[15] = 1; }

template <typename T>
void Invert_M(T M[16], T M_inv[16])
{
  int i, j;
  double** M_d = new T*[4];  T** M_inv_d = new T*[4];
  for (i = 0; i < 4; i++) { M_d[i] = new T[4]; M_inv_d[i] = new T[4]; for (j = 0; j < 4; j++) M_d[i][j] = M[j + i * 4]; }

  vtkMath::InvertMatrix(M_d, M_inv_d, 4);

  for (i = 0; i < 4; i++) { for (j = 0; j < 4; j++) M_inv[j + i * 4] = M_inv_d[i][j]; delete[] M_d[i]; delete[] M_inv_d[i]; }
  delete[] M_d; delete[] M_inv_d;
}

template <typename T>
void Transform(T pin[3], T pout[3], T M[16]) { pout[0] = M[0] * pin[0] + M[1] * pin[1] + M[2] * pin[2] + M[3]; pout[1] = M[4] * pin[0] + M[5] * pin[1] + M[6] * pin[2] + M[7]; pout[2] = M[8] * pin[0] + M[9] * pin[1] + M[10] * pin[2] + M[11]; }

void Transform(vtkPolyData* model, vtkPolyData* model_out, double M[16])
{
  vtkMatrix4x4* M44 = vtkMatrix4x4::New();
  for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) M44->SetElement(i, j, (double)M[4 * i + j]);

  vtkMatrixToHomogeneousTransform* transform = vtkMatrixToHomogeneousTransform::New();
  transform->SetInput(M44);
  transform->Update();
  vtkTransformPolyDataFilter* transf = vtkTransformPolyDataFilter::New();
  transf->SetTransform(transform);
  transf->SetInputData(model);
  transf->Update();
  model_out->DeepCopy(transf->GetOutput());
  SafeDelete(transform);
  SafeDelete(transf);
  SafeDelete(M44);
}

bool SaveHomogeneousMatrixToFile(const std::string& filename, double M[16])
{
  std::ofstream ofs(filename.c_str());

  if (ofs.is_open())
  {
    ofs << std::setprecision(16);
    ofs << M[0] << " " << M[1] << " " << M[2] << " ";
    ofs << M[4] << " " << M[5] << " " << M[6] << " ";
    ofs << M[8] << " " << M[9] << " " << M[10] << std::endl;
    ofs << M[3] << " " << M[7] << " " << M[11];
    ofs.close();
    return true;
  }

  return false;
}

template <unsigned int VDimension>
class BresenhamLine // ripped from itk, current version too old 
{
public:
  typedef BresenhamLine                      Self;
  // This defines the line direction
  typedef itk::Vector<float, VDimension>          LType;
  typedef itk::Offset<VDimension>                 OffsetType;
  typedef itk::Index<VDimension>                  IndexType;
  typedef std::vector<OffsetType>            OffsetArray;

  typedef typename IndexType::IndexValueType IndexValueType;

  // constructors
  BresenhamLine() {}
  ~BresenhamLine() {}

  OffsetArray BuildLine(LType Direction, unsigned int length)
  {
    // copied from the line iterator
    // The dimension with the largest difference between start and end
    unsigned int m_MainDirection;

    // Accumulated error for the other dimensions
    IndexType m_AccumulateError;

    // Increment for the error for each step. Two times the difference between
    // start and end
    IndexType m_IncrementError;

    // If enough is accumulated for a dimension, the index has to be
    // incremented. Will be the number of pixels in the line
    IndexType m_MaximalError;

    // Direction of increment. -1 or 1
    IndexType m_OverflowIncrement;

    // After an overflow, the accumulated error is reduced again. Will be
    // two times the number of pixels in the line
    IndexType m_ReduceErrorAfterIncrement;

    OffsetArray result(length);

    IndexType m_CurrentImageIndex, StartIndex, LastIndex;
    Direction.Normalize();
    // we are going to start at 0
    m_CurrentImageIndex.Fill(0);
    StartIndex.Fill(0);
    for (unsigned i = 0; i < VDimension; i++)
    {
      LastIndex[i] = (IndexValueType)(length*Direction[i]);
    }
    // Find the dominant direction
    IndexValueType maxDistance = 0;
    unsigned int maxDistanceDimension = 0;
    for (unsigned i = 0; i < VDimension; i++)
    {
      IndexValueType distance = abs(LastIndex[i]);
      if (distance > maxDistance)
      {
        maxDistance = distance;
        maxDistanceDimension = i;
      }
      m_IncrementError[i] = 2 * distance;
      m_OverflowIncrement[i] = (LastIndex[i] < 0 ? -1 : 1);
    }
    m_MainDirection = maxDistanceDimension;
    m_MaximalError.Fill(maxDistance);
    m_ReduceErrorAfterIncrement.Fill(2 * maxDistance);
    m_AccumulateError.Fill(0);
    unsigned int steps = 1;
    result[0] = m_CurrentImageIndex - StartIndex;
    while (steps < length)
    {
      // This part is from ++ in LineConstIterator
      // We need to modify m_AccumulateError, m_CurrentImageIndex, m_IsAtEnd
      for (unsigned int i = 0; i < VDimension; ++i)
      {
        if (i == m_MainDirection)
        {

          m_CurrentImageIndex[i] += m_OverflowIncrement[i];
        }
        else
        {
          m_AccumulateError[i] += m_IncrementError[i];
          if (m_AccumulateError[i] >= m_MaximalError[i])
          {
            m_CurrentImageIndex[i] += m_OverflowIncrement[i];
            m_AccumulateError[i] -= m_ReduceErrorAfterIncrement[i];
          }
        }
      }

      result[steps] = m_CurrentImageIndex - StartIndex; // produce an offset

      ++steps;
    }
    return(result);
  }

};

void Line3D_UC
(
  vtkStructuredPoints* vol,
  double iX0, double iY0, double iZ0,
  double iX1, double iY1, double iZ1)
{
  unsigned char* ptr = (unsigned char*)vol->GetScalarPointer();
  // starting point of line
  double iX = iX0, iY = iY0, iZ = iZ0;

  // direction of line
  double iDx = iX1 - iX0, iDy = iY1 - iY0, iDz = iZ1 - iZ0;

  // increment or decrement depending on direction of line
  double iSx = (iDx > 0 ? 0.05 : (iDx < 0 ? -0.05 : 0));
  double iSy = (iDy > 0 ? 0.05 : (iDy < 0 ? -0.05 : 0));
  double iSz = (iDz > 0 ? 0.05 : (iDz < 0 ? -0.05 : 0));

  // decision parameters for voxel selection
  if (iDx < 0) iDx = -iDx;
  if (iDy < 0) iDy = -iDy;
  if (iDz < 0) iDz = -iDz;
  double iAx = 2 * iDx, iAy = 2 * iDy, iAz = 2 * iDz;
  double iDecX, iDecY, iDecZ;

  // determine largest direction component, single-step related variable
  double iMax = iDx;
  int iVar = 0;
  if (iDy > iMax) { iMax = iDy; iVar = 1; }
  if (iDz > iMax) { iVar = 2; }

  // traverse Bresenham line
  switch (iVar)
  {
  case 0:  // single-step in iX-direction
    iDecY = iAy - iDx;
    iDecZ = iAz - iDx;
    for (/**/; /**/; iX += iSx, iDecY += iAy, iDecZ += iAz)
    {
      // process voxel
      int ptId = vol->FindPoint(iX, iY, iZ);
      if (ptId != -1) *(ptr + ptId) = 255;

      // take Bresenham step
      if (iX + 0.05 > iX1 && iX - 0.05 < iX1) break;
      if (iDecY >= 0) { iDecY -= iAx; iY += iSy; }
      if (iDecZ >= 0) { iDecZ -= iAx; iZ += iSz; }
    }
    break;
  case 1:  // single-step in iY-direction
    iDecX = iAx - iDy;
    iDecZ = iAz - iDy;
    for (/**/; /**/; iY += iSy, iDecX += iAx, iDecZ += iAz)
    {
      // process voxel
      int ptId = vol->FindPoint(iX, iY, iZ);
      if (ptId != -1) *(ptr + ptId) = 255;

      // take Bresenham step
      if (iY + 0.05 > iY1 && iY - 0.05 < iY1) break;
      if (iDecX >= 0) { iDecX -= iAy; iX += iSx; }
      if (iDecZ >= 0) { iDecZ -= iAy; iZ += iSz; }
    }
    break;
  case 2:  // single-step in iZ-direction
    iDecX = iAx - iDz;
    iDecY = iAy - iDz;
    for (/**/; /**/; iZ += iSz, iDecX += iAx, iDecY += iAy)
    {
      // process voxel
      int ptId = vol->FindPoint(iX, iY, iZ);
      if (ptId != -1) *(ptr + ptId) = 255;

      // take Bresenham step
      if (iZ + 0.05 > iZ1 && iZ - 0.05 < iZ1) break;
      if (iDecX >= 0) { iDecX -= iAz; iX += iSx; }
      if (iDecY >= 0) { iDecY -= iAz; iY += iSy; }
    }
    break;
  }
}

void GetModelBBox(double tol, vtkPolyData* model, double min[3], double max[3])
{
  double dblp[3];
  min[0] = std::numeric_limits< double >::max();
  min[1] = std::numeric_limits< double >::max();
  min[2] = std::numeric_limits< double >::max();

  max[0] = -std::numeric_limits< double >::max();
  max[1] = -std::numeric_limits< double >::max();
  max[2] = -std::numeric_limits< double >::max();

  for (int i = 0; i < model->GetNumberOfPoints(); i++)
  {
    model->GetPoint(i, dblp);
    if (dblp[0] < min[0]) min[0] = dblp[0]; if (dblp[1] < min[1]) min[1] = dblp[1]; if (dblp[2] < min[2]) min[2] = dblp[2];
    if (dblp[0] > max[0]) max[0] = dblp[0]; if (dblp[1] > max[1]) max[1] = dblp[1]; if (dblp[2] > max[2]) max[2] = dblp[2];
  }

  min[0] -= tol; min[1] -= tol; min[2] -= tol;
  max[0] += tol; max[1] += tol; max[2] += tol;
}

bool GetModelFittedOrientedBBox(
  double tol,
  double epsilon,
  double* points,
  unsigned int numberOfPoints,
  double minb[3], double maxb[3], // min and max of obbox in *local* space
  double model2obbox_M[16], // model to local obbox space coordinates system transform 
  double obbox2model_M[16]) // local obbox to model space coordinates system transform 
{
  double **gpoints_pters = new double*[numberOfPoints];
  for (unsigned int i = 0; i < numberOfPoints; ++i)
  {
    gpoints_pters[i] = points + 3 * i;
  }

  bool bError = true;
  gdiam_bbox  bb = gdiam_approx_mvbb(gpoints_pters, numberOfPoints, epsilon, bError);

  delete[] gpoints_pters;

  if (bError) return false;

  gdiam_point_t gpt;

  double bounds[8][3];

  bb.get_vertex(0, 0, 0, gpt); bounds[0][0] = gpt[0]; bounds[0][1] = gpt[1]; bounds[0][2] = gpt[2];
  bb.get_vertex(1, 0, 0, gpt); bounds[1][0] = gpt[0]; bounds[1][1] = gpt[1]; bounds[1][2] = gpt[2];
  bb.get_vertex(1, 1, 0, gpt); bounds[2][0] = gpt[0]; bounds[2][1] = gpt[1]; bounds[2][2] = gpt[2];
  bb.get_vertex(0, 1, 0, gpt); bounds[3][0] = gpt[0]; bounds[3][1] = gpt[1]; bounds[3][2] = gpt[2];
  bb.get_vertex(0, 0, 1, gpt); bounds[4][0] = gpt[0]; bounds[4][1] = gpt[1]; bounds[4][2] = gpt[2];
  bb.get_vertex(1, 0, 1, gpt); bounds[5][0] = gpt[0]; bounds[5][1] = gpt[1]; bounds[5][2] = gpt[2];
  bb.get_vertex(1, 1, 1, gpt); bounds[6][0] = gpt[0]; bounds[6][1] = gpt[1]; bounds[6][2] = gpt[2];
  bb.get_vertex(0, 1, 1, gpt); bounds[7][0] = gpt[0]; bounds[7][1] = gpt[1]; bounds[7][2] = gpt[2];


  double veci[3] = { bounds[1][0] - bounds[0][0], bounds[1][1] - bounds[0][1], bounds[1][2] - bounds[0][2] };
  double vecj[3] = { bounds[3][0] - bounds[0][0], bounds[3][1] - bounds[0][1], bounds[3][2] - bounds[0][2] };
  double veck[3] = { bounds[4][0] - bounds[0][0], bounds[4][1] - bounds[0][1], bounds[4][2] - bounds[0][2] };

  double li = norm<double>(veci) + 2.0*tol; double lj = norm<double>(vecj) + 2.0*tol; double lk = norm<double>(veck) + 2.0*tol;
  normalizeVec<double>(veci); normalizeVec<double>(vecj); normalizeVec<double>(veck);

  minb[0] = 0;
  minb[1] = 0;
  minb[2] = 0;

  maxb[0] = li;
  maxb[1] = lj;
  maxb[2] = lk;

  Identity<double>(obbox2model_M); Identity<double>(model2obbox_M);
  obbox2model_M[0] = veci[0]; obbox2model_M[4] = veci[1]; obbox2model_M[8] = veci[2];
  obbox2model_M[1] = vecj[0]; obbox2model_M[5] = vecj[1]; obbox2model_M[9] = vecj[2];
  obbox2model_M[2] = veck[0]; obbox2model_M[6] = veck[1]; obbox2model_M[10] = veck[2];
  obbox2model_M[3] = bounds[0][0] - veci[0] * tol - vecj[0] * tol - veck[0] * tol;
  obbox2model_M[7] = bounds[0][1] - veci[1] * tol - vecj[1] * tol - veck[1] * tol;
  obbox2model_M[11] = bounds[0][2] - veci[2] * tol - vecj[2] * tol - veck[2] * tol;

  // test matrix det before inverting it
  vtkSmartPointer< vtkMatrix4x4 > mat = vtkSmartPointer< vtkMatrix4x4 >::New();
  mat->DeepCopy(obbox2model_M);
  double det = mat->Determinant();

  if (abs(det) < 0.000001) return false; // non invertible!

  Invert_M<double>(obbox2model_M, model2obbox_M);

  return true;
}

bool GetModelFittedOrientedBBoxWithVtkOBBTree(
  double tol, // padding
  vtkPolyData* model,
  double minb[3], double maxb[3], // min and max of obbox in *local* space
  double model2obbox_M[16], // model to local obbox space coordinates system transform 
  double obbox2model_M[16]) // local obbox to model space coordinates system transform 
{
  vtkSmartPointer<vtkOBBTree> tree =
    vtkSmartPointer<vtkOBBTree>::New();
  tree->SetDataSet(model);
  tree->BuildLocator();

  double veci[3];
  double vecj[3];
  double veck[3];

  double corner[3]; double size[3];
  tree->ComputeOBB(model, corner, veci, vecj, veck, size);

  double li = norm< double >(veci) + 2.0*tol; normalizeVec< double >(veci);
  double lj = norm< double >(vecj) + 2.0*tol; normalizeVec< double >(vecj);
  double lk = norm< double >(veck) + 2.0*tol; normalizeVec< double >(veck);

  minb[0] = 0;
  minb[1] = 0;
  minb[2] = 0;

  maxb[0] = li;
  maxb[1] = lj;
  maxb[2] = lk;

  Identity<double>(obbox2model_M); Identity<double>(model2obbox_M);
  obbox2model_M[0] = veci[0]; obbox2model_M[4] = veci[1]; obbox2model_M[8] = veci[2];
  obbox2model_M[1] = vecj[0]; obbox2model_M[5] = vecj[1]; obbox2model_M[9] = vecj[2];
  obbox2model_M[2] = veck[0]; obbox2model_M[6] = veck[1]; obbox2model_M[10] = veck[2];
  obbox2model_M[3] = corner[0] - veci[0] * tol - vecj[0] * tol - veck[0] * tol;
  obbox2model_M[7] = corner[1] - veci[1] * tol - vecj[1] * tol - veck[1] * tol;
  obbox2model_M[11] = corner[2] - veci[2] * tol - vecj[2] * tol - veck[2] * tol;

  // test matrix det before inverting it
  vtkSmartPointer< vtkMatrix4x4 > mat = vtkSmartPointer< vtkMatrix4x4 >::New();
  mat->DeepCopy(obbox2model_M);
  double det = mat->Determinant();

  if (abs(det) < 0.000001) return false; // non invertible!

  Invert_M<double>(obbox2model_M, model2obbox_M);

  return true;
}

bool GetModelFittedOrientedBBox(
  double tol,
  double epsilon,
  vtkPolyData* model,
  double minb[3], double maxb[3], // min and max of obbox in *local* space
  double model2obbox_M[16], // model to local obbox space coordinates system transform 
  double obbox2model_M[16]) // local obbox to model space coordinates system transform 
{
  // extract oriented fitted bbox
  double* gpoints = new double[model->GetNumberOfPoints() * 3];
  unsigned int cter = 0;
  for (unsigned int k = 0; k < model->GetNumberOfPoints(); ++k)
  {
    double* p = model->GetPoint(k);
    gpoints[cter] = p[0]; ++cter;
    gpoints[cter] = p[1]; ++cter;
    gpoints[cter] = p[2]; ++cter;
  }

  // First attempt with gdiam
  bool b = false;

  if (!bForceUseOfOBBTree)
    b = GetModelFittedOrientedBBox(tol, epsilon, gpoints, model->GetNumberOfPoints(), minb, maxb, model2obbox_M, obbox2model_M);

  if (!b)
  {
    if (!bForceUseOfOBBTree)
      std::cerr << "Unable to compute oriented bounding box. Attempting to fix this by using an alternative method...";
    b = GetModelFittedOrientedBBoxWithVtkOBBTree(tol, model, minb, maxb, model2obbox_M, obbox2model_M);

    if (!b)
    {
      if (!bForceUseOfOBBTree) std::cerr << "failure. Exiting\n";
      delete[] gpoints;
      return false;
    }
    else
    {
      if (!bForceUseOfOBBTree) std::cerr << "success.\n";
    }
  }

  delete[] gpoints;
  return true;
}

/*
* Correspondences: x,y,z position in oriented bbox frame
* 0,0,0 ->bounds[0]
* 1,0,0 ->bounds[1]
* 1,1,0 ->bounds[2]
* 0,1,0 ->bounds[3]
* 0,0,1 ->bounds[4]
* 1,0,1 ->bounds[5]
* 1,1,1 ->bounds[6]
* 0,1,1 ->bounds[7]
*/
vtkPolyData* CreateOrientedBBoxPolyData(double bnds[8][3])
{
  vtkPolyData* PolyData = vtkPolyData::New();
  vtkPoints* pts = vtkPoints::New(VTK_FLOAT);
  pts->SetNumberOfPoints(8);
  PolyData->Allocate(12); // 2 triangles per face

  for (unsigned int i = 0; i < 8; i++)
  {
    pts->SetPoint(i, bnds[i][0], bnds[i][1], bnds[i][2]);
  } // copy points

  PolyData->SetPoints(pts);

  vtkTriangle* triangle = vtkTriangle::New();

  triangle->GetPointIds()->SetId(0, 0);
  triangle->GetPointIds()->SetId(1, 1);
  triangle->GetPointIds()->SetId(2, 5);
  PolyData->InsertNextCell(triangle->GetCellType(), triangle->GetPointIds()); // Insert triangle

  triangle->GetPointIds()->SetId(0, 0);
  triangle->GetPointIds()->SetId(1, 5);
  triangle->GetPointIds()->SetId(2, 4);
  PolyData->InsertNextCell(triangle->GetCellType(), triangle->GetPointIds()); // Insert triangle

  triangle->GetPointIds()->SetId(0, 1);
  triangle->GetPointIds()->SetId(1, 2);
  triangle->GetPointIds()->SetId(2, 6);
  PolyData->InsertNextCell(triangle->GetCellType(), triangle->GetPointIds()); // Insert triangle

  triangle->GetPointIds()->SetId(0, 1);
  triangle->GetPointIds()->SetId(1, 6);
  triangle->GetPointIds()->SetId(2, 5);
  PolyData->InsertNextCell(triangle->GetCellType(), triangle->GetPointIds()); // Insert triangle

  triangle->GetPointIds()->SetId(0, 4);
  triangle->GetPointIds()->SetId(1, 5);
  triangle->GetPointIds()->SetId(2, 6);
  PolyData->InsertNextCell(triangle->GetCellType(), triangle->GetPointIds()); // Insert triangle

  triangle->GetPointIds()->SetId(0, 4);
  triangle->GetPointIds()->SetId(1, 6);
  triangle->GetPointIds()->SetId(2, 7);
  PolyData->InsertNextCell(triangle->GetCellType(), triangle->GetPointIds()); // Insert triangle

  triangle->GetPointIds()->SetId(0, 0);
  triangle->GetPointIds()->SetId(1, 7);
  triangle->GetPointIds()->SetId(2, 3);
  PolyData->InsertNextCell(triangle->GetCellType(), triangle->GetPointIds()); // Insert triangle

  triangle->GetPointIds()->SetId(0, 0);
  triangle->GetPointIds()->SetId(1, 4);
  triangle->GetPointIds()->SetId(2, 7);
  PolyData->InsertNextCell(triangle->GetCellType(), triangle->GetPointIds()); // Insert triangle

  triangle->GetPointIds()->SetId(0, 2);
  triangle->GetPointIds()->SetId(1, 7);
  triangle->GetPointIds()->SetId(2, 3);
  PolyData->InsertNextCell(triangle->GetCellType(), triangle->GetPointIds()); // Insert triangle

  triangle->GetPointIds()->SetId(0, 2);
  triangle->GetPointIds()->SetId(1, 6);
  triangle->GetPointIds()->SetId(2, 7);
  PolyData->InsertNextCell(triangle->GetCellType(), triangle->GetPointIds()); // Insert triangle

  triangle->GetPointIds()->SetId(0, 0);
  triangle->GetPointIds()->SetId(1, 1);
  triangle->GetPointIds()->SetId(2, 2);
  PolyData->InsertNextCell(triangle->GetCellType(), triangle->GetPointIds()); // Insert triangle

  triangle->GetPointIds()->SetId(0, 0);
  triangle->GetPointIds()->SetId(1, 2);
  triangle->GetPointIds()->SetId(2, 3);
  PolyData->InsertNextCell(triangle->GetCellType(), triangle->GetPointIds()); // Insert triangle

  PolyData->Modified();
  SafeDelete(triangle);

  vtkPolyDataNormals *normals = vtkPolyDataNormals::New();
  normals->SetInputData(PolyData);
  normals->SplittingOff();
  normals->ConsistencyOn();
  normals->ComputePointNormalsOn();
  normals->Update();

  SafeDelete(PolyData);
  PolyData = vtkPolyData::New();
  PolyData->DeepCopy(normals->GetOutput());
  SafeDelete(normals);

  return PolyData;
}

// Get vtkpolydata of oriented bounding box in *model* space
vtkPolyData* CreateOrientedBBoxPolyDataInModelSpace(
  double minb[3], double maxb[3], // min and max of obbox in *local* space
  double obbox2model_M[16]) // local obbox to model space coordinates system transform 
{
  double bnds[8][3];

  double vi[3];
  double vo[3];

  size_t size_ = 3 * sizeof(double);

  // bnds[0]
  memcpy(vi, minb, size_);
  Transform< double >(vi, vo, obbox2model_M);
  memcpy(bnds[0], vo, size_);

  // bnds[1]
  vi[0] = minb[0] + maxb[0]; vi[1] = minb[1]; vi[2] = minb[2];
  Transform< double >(vi, vo, obbox2model_M);
  memcpy(bnds[1], vo, size_);

  // bnds[2]
  vi[0] = minb[0] + maxb[0]; vi[1] = minb[1] + maxb[1]; vi[2] = minb[2];
  Transform< double >(vi, vo, obbox2model_M);
  memcpy(bnds[2], vo, size_);

  // bnds[3]
  vi[0] = minb[0]; vi[1] = minb[1] + maxb[1]; vi[2] = minb[2];
  Transform< double >(vi, vo, obbox2model_M);
  memcpy(bnds[3], vo, size_);

  // bnds[4]
  vi[0] = minb[0]; vi[1] = minb[1]; vi[2] = minb[2] + maxb[2];
  Transform< double >(vi, vo, obbox2model_M);
  memcpy(bnds[4], vo, size_);

  // bnds[5]
  vi[0] = minb[0] + maxb[0]; vi[1] = minb[1]; vi[2] = minb[2] + maxb[2];
  Transform< double >(vi, vo, obbox2model_M);
  memcpy(bnds[5], vo, size_);

  // bnds[6]
  vi[0] = minb[0] + maxb[0]; vi[1] = minb[1] + maxb[1]; vi[2] = minb[2] + maxb[2];
  Transform< double >(vi, vo, obbox2model_M);
  memcpy(bnds[6], vo, size_);

  // bnds[7]
  vi[0] = minb[0]; vi[1] = minb[1] + maxb[1]; vi[2] = minb[2] + maxb[2];
  Transform< double >(vi, vo, obbox2model_M);
  memcpy(bnds[7], vo, size_);

  return CreateOrientedBBoxPolyData(bnds);
}

// Rasterize the surface of a poly data mesh
vtkStructuredPoints* fillStencil(double* bounds, double* spac, vtkPolyData* model)
{
  int i;
  double min[3], max[3];

  //double dblp[3];
  double spacing[3];
  double tol = 10;

  if (spac != NULL) memcpy(spacing, spac, 3 * sizeof(double));
  else { spacing[0] = 1; spacing[1] = 1; spacing[2] = 1; }

  if (bounds == NULL)
  {
    GetModelBBox(tol, model, min, max);
  }
  else
  {
    min[0] = bounds[0]; max[0] = bounds[1];
    min[1] = bounds[2]; max[1] = bounds[3];
    min[2] = bounds[4]; max[2] = bounds[5];
  }

  double minspac = spacing[0];
  minspac = spacing[1] < minspac ? spacing[1] : minspac;
  minspac = spacing[2] < minspac ? spacing[2] : minspac;

  vtkSmartPointer<vtkImageData> volume =
    vtkSmartPointer<vtkImageData>::New();

  volume->SetOrigin(min[0], min[1], min[2]);
  volume->SetDimensions((int)((max[0] - min[0]) / spacing[0] + 0.5), (int)((max[1] - min[1]) / spacing[1] + 0.5), (int)((max[2] - min[2]) / spacing[2] + 0.5));
  volume->SetSpacing(spacing);
  volume->AllocateScalars(VTK_UNSIGNED_CHAR, 1);

  unsigned char* ptr_uc = (unsigned char*)volume->GetScalarPointer(); for (i = 0; i < volume->GetNumberOfPoints(); i++) *(ptr_uc + i) = 255;


  // polygonal data --> image stencil:
  vtkSmartPointer<vtkPolyDataToImageStencil> pol2stenc =
    vtkSmartPointer<vtkPolyDataToImageStencil>::New();
#if VTK_MAJOR_VERSION <= 5
  pol2stenc->SetInput(model);
#else
  pol2stenc->SetInputData(model);
#endif
  pol2stenc->SetOutputOrigin(min);
  pol2stenc->SetOutputSpacing(spacing);
  pol2stenc->SetOutputWholeExtent(volume->GetExtent());
  pol2stenc->Update();

  // cut the corresponding white image and set the background:
  vtkSmartPointer<vtkImageStencil> imgstenc =
    vtkSmartPointer<vtkImageStencil>::New();
#if VTK_MAJOR_VERSION <= 5
  imgstenc->SetInput(volume);
  imgstenc->SetStencil(pol2stenc->GetOutput());
#else
  imgstenc->SetInputData(volume);
  imgstenc->SetStencilConnection(pol2stenc->GetOutputPort());
#endif
  imgstenc->ReverseStencilOff();
  imgstenc->SetBackgroundValue(0);
  imgstenc->Update();

  vtkSmartPointer<vtkImageToStructuredPoints> im2vpoints = vtkSmartPointer<vtkImageToStructuredPoints>::New();
  im2vpoints->SetInputConnection(imgstenc->GetOutputPort());
  im2vpoints->Update();

  vtkStructuredPoints* ov = vtkStructuredPoints::New();
  ov->DeepCopy(im2vpoints->GetOutput());
  return ov;
}


// Rasterize the surface of a poly data mesh
vtkStructuredPoints* voxelize(double* bounds, double* spac, vtkPolyData* model)
{

  int i;
  double min[3], max[3];

  //double dblp[3];
  double spacing[3];
  double tol = 10;

  if (spac != NULL) memcpy(spacing, spac, 3 * sizeof(double));
  else { spacing[0] = 1; spacing[1] = 1; spacing[2] = 1; }

  if (bounds == NULL)
  {
    GetModelBBox(tol, model, min, max);
  }
  else
  {
    min[0] = bounds[0]; max[0] = bounds[1];
    min[1] = bounds[2]; max[1] = bounds[3];
    min[2] = bounds[4]; max[2] = bounds[5];
  }

  double minspac = spacing[0];
  minspac = spacing[1] < minspac ? spacing[1] : minspac;
  minspac = spacing[2] < minspac ? spacing[2] : minspac;

  vtkStructuredPoints* volume = vtkStructuredPoints::New();
  volume->SetOrigin(min[0], min[1], min[2]);
  volume->SetDimensions((int)((max[0] - min[0]) / spacing[0] + 0.5), (int)((max[1] - min[1]) / spacing[1] + 0.5), (int)((max[2] - min[2]) / spacing[2] + 0.5));
  volume->SetSpacing(spacing);
  volume->AllocateScalars(VTK_UNSIGNED_CHAR, 1);

  unsigned char* ptr_uc = (unsigned char*)volume->GetScalarPointer(); for (i = 0; i < volume->GetNumberOfPoints(); i++) *(ptr_uc + i) = 0;
  int ids[3];
  double p1[3], p2[3], p3[3], p1p2[3], p1p3[3], p2p1[3], p2p3[3], p3p1[3], p3p2[3], ps[3], pt[3];


  //voxelize triangles using 3D bresenham
  for (i = 0; i < model->GetNumberOfCells(); i++)
  {
    ids[0] = model->GetCell(i)->GetPointIds()->GetId(0); ids[1] = model->GetCell(i)->GetPointIds()->GetId(1); ids[2] = model->GetCell(i)->GetPointIds()->GetId(2);
    p1[0] = model->GetPoint(ids[0])[0]; p1[1] = model->GetPoint(ids[0])[1]; p1[2] = model->GetPoint(ids[0])[2]; p2[0] = model->GetPoint(ids[1])[0]; p2[1] = model->GetPoint(ids[1])[1]; p2[2] = model->GetPoint(ids[1])[2]; p3[0] = model->GetPoint(ids[2])[0]; p3[1] = model->GetPoint(ids[2])[1]; p3[2] = model->GetPoint(ids[2])[2];
    p1p2[0] = p2[0] - p1[0]; p1p2[1] = p2[1] - p1[1]; p1p2[2] = p2[2] - p1[2]; p1p3[0] = p3[0] - p1[0]; p1p3[1] = p3[1] - p1[1]; p1p3[2] = p3[2] - p1[2];
    p2p1[0] = p1[0] - p2[0]; p2p1[1] = p1[1] - p2[1]; p2p1[2] = p1[2] - p2[2]; p2p3[0] = p3[0] - p2[0]; p2p3[1] = p3[1] - p2[1]; p2p3[2] = p3[2] - p2[2];
    p3p1[0] = p1[0] - p3[0]; p3p1[1] = p1[1] - p3[1]; p3p1[2] = p1[2] - p3[2]; p3p2[0] = p2[0] - p3[0]; p3p2[1] = p2[1] - p3[1]; p3p2[2] = p2[2] - p3[2];


#if 0

    BresenhamLine< 3 > Line1, Line2;
    BresenhamLine< 3 >::LType dir1, dir2;
    dir1[0] = p1p2[0]; dir1[1] = p1p2[1]; dir1[2] = p1p2[2];
    dir2[0] = p1p3[0]; dir2[1] = p1p3[1]; dir2[2] = p1p3[2];
    unsigned int length1 = dir1.GetNorm();
    unsigned int length2 = dir2.GetNorm();
    dir1 /= (float)length1;
    dir2 /= (float)length2;

    // let's create a line for each side of the triangle
    BresenhamLine< 3 >::OffsetArray offsets1 = Line1.BuildLine(dir1, length1);
    BresenhamLine< 3 >::OffsetArray offsets2 = Line2.BuildLine(dir2, length2);

    // traverse 
#else
    const double p1p2x = p1p2[0];
    const double p1p2y = p1p2[1];
    const double p1p2z = p1p2[2];

    const double p1p3x = p1p3[0];
    const double p1p3y = p1p3[1];
    const double p1p3z = p1p3[2];

    const double p2p1x = p2p1[0];
    const double p2p1y = p2p1[1];
    const double p2p1z = p2p1[2];

    const double p2p3x = p2p3[0];
    const double p2p3y = p2p3[1];
    const double p2p3z = p2p3[2];

    const double p3p1x = p3p1[0];
    const double p3p1y = p3p1[1];
    const double p3p1z = p3p1[2];

    const double p3p2x = p3p2[0];
    const double p3p2y = p3p2[1];
    const double p3p2z = p3p2[2];

    // compute minimal step
    // ensure that it is below min spacing direction to avoid holes
    double step = 0.1;
    while (step * p1p2x > minspac) step *= 0.5;
    while (step * p1p2y > minspac) step *= 0.5;
    while (step * p1p2z > minspac) step *= 0.5;
    while (step * p1p3x > minspac) step *= 0.5;
    while (step * p1p3y > minspac) step *= 0.5;
    while (step * p1p3z > minspac) step *= 0.5;
    while (step * p2p1x > minspac) step *= 0.5;
    while (step * p2p1y > minspac) step *= 0.5;
    while (step * p2p1z > minspac) step *= 0.5;
    while (step * p2p3x > minspac) step *= 0.5;
    while (step * p2p3y > minspac) step *= 0.5;
    while (step * p2p3z > minspac) step *= 0.5;
    while (step * p3p1x > minspac) step *= 0.5;
    while (step * p3p1y > minspac) step *= 0.5;
    while (step * p3p1z > minspac) step *= 0.5;
    while (step * p3p2x > minspac) step *= 0.5;
    while (step * p3p2y > minspac) step *= 0.5;
    while (step * p3p2z > minspac) step *= 0.5;

    for (double s = 0.0; s <= 1.0; s += step)
    {
      ps[0] = p1[0] + s*p1p2x; ps[1] = p1[1] + s*p1p2y; ps[2] = p1[2] + s*p1p2z;
      pt[0] = p1[0] + s*p1p3x; pt[1] = p1[1] + s*p1p3y; pt[2] = p1[2] + s*p1p3z;
      Line3D_UC(volume, ps[0], ps[1], ps[2], pt[0], pt[1], pt[2]);
      ps[0] = p2[0] + s*p2p1x; ps[1] = p2[1] + s*p2p1y; ps[2] = p2[2] + s*p2p1z;
      pt[0] = p2[0] + s*p2p3x; pt[1] = p2[1] + s*p2p3y; pt[2] = p2[2] + s*p2p3z;
      Line3D_UC(volume, ps[0], ps[1], ps[2], pt[0], pt[1], pt[2]);
      ps[0] = p3[0] + s*p3p1x; ps[1] = p3[1] + s*p3p1y; ps[2] = p3[2] + s*p3p1z;
      pt[0] = p3[0] + s*p3p2x; pt[1] = p3[1] + s*p3p2y; pt[2] = p3[2] + s*p3p2z;
      Line3D_UC(volume, ps[0], ps[1], ps[2], pt[0], pt[1], pt[2]);
    }
#endif
  }
  
  return volume;
}

// Fill Interior
vtkStructuredPoints* fillInterior(vtkStructuredPoints* input)
{
  int *dim = input->GetDimensions();
  unsigned char* ptr = (unsigned char*)input->GetScalarPointer();

  vtkStructuredPoints* volume = vtkStructuredPoints::New();
  volume->SetOrigin(input->GetOrigin());
  volume->SetDimensions(input->GetDimensions());
  volume->SetSpacing(input->GetSpacing());
  volume->AllocateScalars(VTK_UNSIGNED_CHAR, 1);

  // Fill the interior: fill the exterior with 0 and the rest (interior) with 255
  // surface has been rasterized and painted with 255
  vtkImageSeedConnectivity* connect = vtkImageSeedConnectivity::New();
  connect->SetInputData(input);
  connect->SetDimensionality(3);
  connect->SetInputConnectValue(0);
  connect->SetOutputConnectedValue(0);
  connect->SetOutputUnconnectedValue(255);

  // place seeds at four of the eight corners of image space
  connect->AddSeed(0, 0, dim[2] - 1);
  connect->AddSeed(0, dim[1] - 1, dim[2] - 1);
  connect->AddSeed(dim[0] - 1, 0, dim[2] - 1);
  connect->AddSeed(dim[0] - 1, dim[1] - 1, dim[2] - 1);
  connect->SetOutput(volume);
  connect->Update();

  SafeDelete(connect);
  return volume;
}

bool isInBounds(double *bounds, double* pt)
{
  if (pt[0] < bounds[0]) return false;
  if (pt[0] > bounds[1]) return false;
  if (pt[1] < bounds[2]) return false;
  if (pt[1] > bounds[3]) return false;
  if (pt[2] < bounds[4]) return false;
  if (pt[2] > bounds[5]) return false;
  return true;
}

void paintOutputImage(vtkStructuredPoints* fillImage,
  vtkStructuredPoints* outputImage,
  unsigned short label)
{
  // longer but more accurate: traverse output image and look in
  // in fill image for the filled areas

  unsigned int nOutput = outputImage->GetNumberOfPoints();
  int nFill = fillImage->GetNumberOfPoints();
  double pt[3];

  unsigned short* optr;
  if (!bUchar) optr = (unsigned short*)outputImage->GetScalarPointer();
  unsigned char* optr_uchar;
  if (bUchar) optr_uchar = (unsigned char*)outputImage->GetScalarPointer();
  unsigned char* fptr = (unsigned char*)fillImage->GetScalarPointer();
  double bounds[6];
  fillImage->GetBounds(bounds);

  for (unsigned int i = 0; i < nOutput; i++)
  {
    // Get output coordinates
    outputImage->GetPoint(i, pt);

    // fast test to test if it is in bounds
    if (!isInBounds(bounds, pt))
    {
      if (bUchar) optr_uchar++; else optr++;
      continue;
    }

    // find id of corresponding position in fill image (nearest neigh interpolation)
    vtkIdType id = fillImage->FindPoint(pt);

    // Get value
    if (id < 0 || id >= nFill)
    {
      if (bUchar) optr_uchar++; else optr++;
      continue;
    }

    unsigned char v = fptr[id];

    // if v is diff from zero then we paint
    if (v != 0)
    {
      if (bUchar) *optr_uchar = (unsigned char)label;
      else *optr = label;
    }
    if (bUchar) optr_uchar++; else optr++;
  }
}

bool writeImageToFile(vtkStructuredPoints* image, const std::string filename)
{
  // write image to file
  if (image->GetScalarType() == VTK_UNSIGNED_CHAR) {
    typedef itk::Image< unsigned char, 3 > ImageType;
    typedef itk::VTKImageToImageFilter< ImageType > VtkToItkFilterType;
    VtkToItkFilterType::Pointer converter = VtkToItkFilterType::New();
    converter->SetInput(image);
    converter->Update();

    typedef itk::ChangeInformationImageFilter< ImageType > FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    filter->SetInput(converter->GetOutput());
    typename ImageType::DirectionType newDirection = converter->GetOutput()->GetDirection();
    newDirection(0, 0) = directionx[0]; newDirection(0, 1) = directiony[0]; newDirection(0, 2) = directionz[0];
    newDirection(1, 0) = directionx[1]; newDirection(1, 1) = directiony[1]; newDirection(1, 2) = directionz[1];
    newDirection(2, 0) = directionx[2]; newDirection(2, 1) = directiony[2]; newDirection(2, 2) = directionz[2];

    filter->SetOutputDirection(newDirection);
    filter->ChangeDirectionOn();
    filter->Update();

    typedef itk::ImageFileWriter< ImageType > WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetInput(filter->GetOutput());
    writer->SetFileName(filename);
    writer->Update();
  }
  else if (image->GetScalarType() == VTK_UNSIGNED_SHORT)
  {
    typedef itk::Image< unsigned short, 3 > ImageType;
    typedef itk::VTKImageToImageFilter< ImageType > VtkToItkFilterType;
    VtkToItkFilterType::Pointer converter = VtkToItkFilterType::New();
    converter->SetInput(image);
    converter->Update();

    typedef itk::ChangeInformationImageFilter< ImageType > FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    filter->SetInput(converter->GetOutput());
    typename ImageType::DirectionType newDirection = converter->GetOutput()->GetDirection();
    newDirection(0, 0) = directionx[0]; newDirection(0, 1) = directiony[0]; newDirection(0, 2) = directionz[0];
    newDirection(1, 0) = directionx[1]; newDirection(1, 1) = directiony[1]; newDirection(1, 2) = directionz[1];
    newDirection(2, 0) = directionx[2]; newDirection(2, 1) = directiony[2]; newDirection(2, 2) = directionz[2];

    filter->SetOutputDirection(newDirection);
    filter->ChangeDirectionOn();
    filter->Update();

    typedef itk::ImageFileWriter< ImageType > WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetInput(filter->GetOutput());
    writer->SetFileName(filename);
    writer->Update();
  }

  return true;
}

bool GetPolyDataFromUGrid(vtkUnstructuredGrid* ugrid, vtkPolyData*& PolyData)
{
  if (ugrid == NULL) return false;
  int i;
  int j = 0;
  for (i = 0; i < ugrid->GetNumberOfCells(); i++)
  {
    j += ugrid->GetCell(i)->GetNumberOfPoints();
  }

  SafeDelete(PolyData);
  PolyData = vtkPolyData::New();
  vtkPoints* pts = vtkPoints::New();
  pts->SetNumberOfPoints(ugrid->GetNumberOfPoints() + ugrid->GetNumberOfCells());
  PolyData->SetPoints(pts);

  PolyData->Allocate(j, j);

  vtkTriangle* triangle = vtkTriangle::New();

  double p[3];
  for (i = 0; i < ugrid->GetNumberOfPoints(); i++)
  {
    ugrid->GetPoint(i, p);
    pts->SetPoint(i, p);
  } // copy points

  for (i = 0; i < ugrid->GetNumberOfCells(); i++)
  {
    triangle->GetPointIds()->SetId(0, ugrid->GetNumberOfPoints() + i);
    for (j = 0; j < ugrid->GetCell(i)->GetNumberOfPoints(); j++)
    {
      triangle->GetPointIds()->SetId(1, ugrid->GetCell(i)->GetPointId(j));
      if (j == ugrid->GetCell(i)->GetNumberOfPoints() - 1)
        triangle->GetPointIds()->SetId(2, ugrid->GetCell(i)->GetPointId(0));
      else
        triangle->GetPointIds()->SetId(2, ugrid->GetCell(i)->GetPointId(j + 1));

      PolyData->InsertNextCell(triangle->GetCellType(), triangle->GetPointIds()); // Insert triangle
    }

    // compute cell center
    double c[3] = { 0,0,0 };

    for (unsigned int k = 0; k < ugrid->GetCell(i)->GetNumberOfPoints(); ++k)
    {
      c[0] += ugrid->GetPoint(ugrid->GetCell(i)->GetPointId(k))[0];
      c[1] += ugrid->GetPoint(ugrid->GetCell(i)->GetPointId(k))[1];
      c[2] += ugrid->GetPoint(ugrid->GetCell(i)->GetPointId(k))[2];
    }

    c[0] /= (double)ugrid->GetCell(i)->GetNumberOfPoints();
    c[1] /= (double)ugrid->GetCell(i)->GetNumberOfPoints();
    c[2] /= (double)ugrid->GetCell(i)->GetNumberOfPoints();

    pts->SetPoint(ugrid->GetNumberOfPoints() + i, c); // insert cell center
  }

  PolyData->SetPoints(pts);
  SafeDelete(triangle);

  vtkPolyDataNormals *normals = vtkPolyDataNormals::New();
  normals->SetInputData(PolyData);
  normals->SplittingOff();
  normals->ConsistencyOn();
  normals->ComputePointNormalsOn();
  normals->Update();

  SafeDelete(PolyData);
  PolyData = vtkPolyData::New();
  PolyData->DeepCopy(normals->GetOutput());
  SafeDelete(normals);

  return true;
}

bool IsTriangulated(vtkPolyData* v)
{
  // traverse all cells. If at least one is not a triangle exit
  for (unsigned int i = 0; i < v->GetNumberOfCells(); ++i)
  {
    vtkCell* cell = v->GetCell(i);
    int np = cell->GetNumberOfPoints();

    if (np != 3)
      return false;
  }

  return true;
}

bool loadMesh(const std::string& filename, vtkSmartPointer< vtkPolyData >& poly, vtkSmartPointer< vtkUnstructuredGrid >& ugrid)
{
  // read mesh
  poly = vtkSmartPointer< vtkPolyData >::New();
  ugrid = vtkSmartPointer< vtkUnstructuredGrid >::New();

  bool bIsPolyData = LoadPolyData(filename, poly) == 0;
  bool bIsUGrid = false;
  if (!bIsPolyData) {
    // vtk unstructured grid?
    bIsUGrid = LoadUnstructuredGridFromVTK(filename, ugrid) == 0;
  }

  return bIsPolyData || bIsUGrid;
}


bool loadPolyData(const std::string& filename, vtkPolyData*& v)
{
  vtkSmartPointer< vtkPolyData > poly = vtkSmartPointer< vtkPolyData >::New();
  vtkSmartPointer< vtkUnstructuredGrid > ugrid = vtkSmartPointer< vtkUnstructuredGrid >::New();

  bool b = loadMesh(filename, poly, ugrid);

  if (b && ugrid.GetPointer()) {
    std::cout << "Filename " << filename << " loaded.\n";
    v->DeepCopy(poly);
  }
  else
  {
    std::cerr << "Cannot read triangular mesh " << filename << "\n";
    return false;
  }

  // test if mesh is triangulated otherwise conversion will fail...!
  bool bIsTri = IsTriangulated(v);

  if (!bIsTri)
  {
    std::cout << "Mesh is NOT triangular! Triangulate it before using this executable.\n";
    return false;
  }

  return true;
}

bool computeOrUpdatePolyDataFittedOrientedBboxFromMesh(
  vtkPolyData* v,
  const double* spacing,
  double* origin,
  int* size,
  double model2obbox_M[16],
  double obbox2model_M[16],
  double tol = 10, bool bUpdate = false)
{
  if (v->GetNumberOfPoints() == 0) return false;

  double min[3];
  double max[3];

  double epsilon = 0.00001;

  if (bUpdate)
  {
    // we need to transform the current image dimensions in model space
    // so that we can compute an updated obbox
    min[0] = origin[0];
    min[1] = origin[1];
    min[2] = origin[2];

    max[0] = size[0] * spacing[0] + min[0];
    max[1] = size[1] * spacing[1] + min[1];
    max[2] = size[2] * spacing[2] + min[2];

    vtkPolyData* imagev = CreateOrientedBBoxPolyDataInModelSpace(min, max, obbox2model_M);

    // Create temp pointer structure to compute obbox
    double* gpoints = new double[v->GetNumberOfPoints() * 3 + imagev->GetNumberOfPoints() * 3];
    unsigned int cter = 0;
    for (unsigned int k = 0; k < v->GetNumberOfPoints(); ++k)
    {
      double* p = v->GetPoint(k);
      gpoints[cter] = p[0]; ++cter;
      gpoints[cter] = p[1]; ++cter;
      gpoints[cter] = p[2]; ++cter;
    }

    for (unsigned int k = 0; k < imagev->GetNumberOfPoints(); ++k)
    {
      double* p = imagev->GetPoint(k);
      gpoints[cter] = p[0]; ++cter;
      gpoints[cter] = p[1]; ++cter;
      gpoints[cter] = p[2]; ++cter;
    }

    // extract oriented fitted bbox encompassing new polydata and image in model space
    bool b = false;

    if (!bForceUseOfOBBTree)
      GetModelFittedOrientedBBox(tol, epsilon, gpoints, v->GetNumberOfPoints() + imagev->GetNumberOfPoints(), min, max, model2obbox_M, obbox2model_M);

    if (!b)
    {
      if (!bForceUseOfOBBTree)	std::cout << "Unable to update oriented bounding box. Attempting to fix this by using an alternative method...";

      // gdiam failed... Let's switch to our backup option with vtkOBBTree
      // for that we need a unique polydata set which encompassses both v and imagev
      vtkSmartPointer< vtkAppendPolyData > append = vtkSmartPointer< vtkAppendPolyData >::New();
      append->AddInputData(v);
      append->AddInputData(imagev);
      append->Update();

      b = GetModelFittedOrientedBBoxWithVtkOBBTree(tol, append->GetOutput(), min, max, model2obbox_M, obbox2model_M);

      if (!bForceUseOfOBBTree)
      {
        if (b)
        {
          std::cout << "success.\n";
        }
        else
        {
          std::cout << "failure.\n";
        }
      }
    }

    SafeDelete(imagev);
    delete[] gpoints;

    if (!b)
    {
      return false;
    }
  }
  else
  {
    // only compute bbox for current model
    bool b = GetModelFittedOrientedBBox(tol, epsilon, v, min, max, model2obbox_M, obbox2model_M);

    if (!b) return false;
  }

  // update image characteristics
  origin[0] = min[0];
  origin[1] = min[1];
  origin[2] = min[2];

  size[0] = (max[0] - min[0]) / spacing[0];
  size[1] = (max[1] - min[1]) / spacing[1];
  size[2] = (max[2] - min[2]) / spacing[2];
  return true;
}


bool computeOrUpdatePolyDataFittedOrientedBboxFromFile(
  const std::string& filename,
  const double* spacing,
  double* origin,
  int* size,
  double model2obbox_M[16],
  double obbox2model_M[16],
  double tol = 10, bool bUpdate = false)
{
  vtkPolyData* v = vtkPolyData::New();
  bool bv = loadPolyData(filename, v);
  if (!bv && v)
  {
    SafeDelete(v);
    return false;
  }

  bool b = computeOrUpdatePolyDataFittedOrientedBboxFromMesh(v, spacing, origin, size, model2obbox_M, obbox2model_M, tol, bUpdate);
  SafeDelete(v);
  return b;
}

bool computeOrUpdatePolyDataBboxFromMesh(
  vtkPolyData* v,
  const double* spacing,
  double* origin,
  int* size,
  double tol = 10, bool bUpdate = false)
{
  if (v->GetNumberOfPoints() == 0) return false;

  int i;
  double min[3] = { std::numeric_limits< double >::max(),std::numeric_limits< double >::max(),std::numeric_limits< double >::max() };
  double max[3] = { -std::numeric_limits< double >::max(),-std::numeric_limits< double >::max(),-std::numeric_limits< double >::max() };

  if (bUpdate)
  {
    min[0] = origin[0];
    min[1] = origin[1];
    min[2] = origin[2];

    max[0] = size[0] * spacing[0] + min[0];
    max[1] = size[1] * spacing[1] + min[1];
    max[2] = size[2] * spacing[2] + min[2];
  }

  double dblp[3];
  for (i = 0; i < v->GetNumberOfPoints(); i++)
  {
    v->GetPoint(i, dblp);
    if (dblp[0] < min[0]) min[0] = dblp[0]; if (dblp[1] < min[1]) min[1] = dblp[1]; if (dblp[2] < min[2]) min[2] = dblp[2];
    if (dblp[0] > max[0]) max[0] = dblp[0]; if (dblp[1] > max[1]) max[1] = dblp[1]; if (dblp[2] > max[2]) max[2] = dblp[2];
  }

  min[0] -= tol; min[1] -= tol; min[2] -= tol;
  max[0] += tol; max[1] += tol; max[2] += tol;

  origin[0] = min[0];
  origin[1] = min[1];
  origin[2] = min[2];

  size[0] = (int)((max[0] - min[0]) / spacing[0] + 0.5);
  size[1] = (int)((max[1] - min[1]) / spacing[1] + 0.5);
  size[2] = (int)((max[2] - min[2]) / spacing[2] + 0.5);
  return true;
}

bool computeOrUpdatePolyDataBboxFromFile(
  const std::string& filename,
  const double* spacing,
  double* origin,
  int* size,
  double tol = 10, bool bUpdate = false)
{
  vtkPolyData* v = vtkPolyData::New();
  bool bv = loadPolyData(filename, v);
  if (!bv && v)
  {
    SafeDelete(v);
    return false;
  }

  bool b = computeOrUpdatePolyDataBboxFromMesh(v, spacing, origin, size, tol, bUpdate);
  SafeDelete(v);
  return b;
}

bool meshToLabel(const std::string& meshFilename, unsigned short label, double* M = NULL)
{
  // read mesh
  vtkPolyData* v = vtkPolyData::New();
  bool bv = loadPolyData(meshFilename, v);

  if (!bv) { SafeDelete(v); return false; }

  std::cout << "Mesh to obbox transform:\n";
  std::cout << M[0] << " " << M[1] << " " << M[2] << " " << M[3] << std::endl;
  std::cout << M[4] << " " << M[5] << " " << M[6] << " " << M[7] << std::endl;
  std::cout << M[8] << " " << M[9] << " " << M[10] << " " << M[11] << std::endl;
  std::cout << M[12] << " " << M[13] << " " << M[14] << " " << M[15] << std::endl;


  // Apply a possible transform to the mesh
  if (M)
  {
    vtkPolyData* vout = vtkPolyData::New();
    Transform(v, vout, M);
    SafeDelete(v);
    v = vout;
  }

  // before any voxelisation lets test
  // if mesh bbox intersects output image, 
  // this should fasten significantly the process
  double min[3];
  double max[3];
  GetModelBBox(10, v, min, max);
  double mc[3] = { 0.5*(min[0] + max[0]),0.5*(min[1] + max[1]) ,0.5*(min[2] + max[2]) };
  double mr[3] = { 0.5*(-min[0] + max[0]),0.5*(-min[1] + max[1]) ,0.5*(-min[2] + max[2]) };

  double mino[3];
  mino[0] = origin[0];
  mino[1] = origin[1];
  mino[2] = origin[2];

  double maxo[3];
  maxo[0] = origin[0] + size[0] * spacing[0];
  maxo[1] = origin[1] + size[1] * spacing[1];
  maxo[2] = origin[2] + size[2] * spacing[2];

  double ic[3] = { 0.5*(mino[0] + maxo[0]),0.5*(mino[1] + maxo[1]) ,0.5*(mino[2] + maxo[2]) };
  double ir[3] = { 0.5*(-mino[0] + maxo[0]),0.5*(-mino[1] + maxo[1]) ,0.5*(-mino[2] + maxo[2]) };

  /*
  bool b =
      min[0] > maxo[0] ||
      max[0] < mino[0] ||
      min[1] > maxo[1] ||
      max[1] < mino[1] ||
      min[2] > maxo[2] ||
      max[2] < mino[2];*/

      // https://studiofreya.com/3d-math-and-physics/simple-aabb-vs-aabb-collision-detection/
  bool bx = abs(mc[0] - ic[0]) <= (mr[0] + ir[0]);
  bool by = abs(mc[1] - ic[1]) <= (mr[1] + ir[1]);
  bool bz = abs(mc[2] - ic[2]) <= (mr[2] + ir[2]);
  bool b = !(bx && by && bz);


  if (b && bTestMeshBboxInImage)
  {
    std::cout << "Mesh " << meshFilename << " is not present in the output image (=no painting)\n";
    SafeDelete(v);
    return true;
  }

  // voxelize mesh surface
  // by passing NULL bounds, a bbox expanded with a tol is computed
  vtkStructuredPoints* surfaceImage;
  vtkStructuredPoints* fillImage;
  if (bUseVtkStencil)
  {
    std::cout << "---Voxelize and fill interior\n";
    fillImage = fillStencil(0, spacing, v);
  }
  else
  {
    std::cout << "---Voxelize\n";
    surfaceImage = voxelize(0, spacing, v);

    // temp
    //writeImageToFile(surfaceImage, "surface.vtk");

    // fill interior
    std::cout << "---fillinterior\n";
    fillImage = fillInterior(surfaceImage);
  }


  // temp
  //writeImageToFile( fillImage, "fill.vtk" );

  // Update out image with sub image of the filled image
  std::cout << "label of mesh " << meshFilename << " is: " << label << std::endl;
  std::cout << "---paintOutputImage\n";
  paintOutputImage(fillImage, oImage, label);

  if (!bUseVtkStencil) SafeDelete(surfaceImage);
  SafeDelete(fillImage);
  SafeDelete(v);
  return true;
}


bool analyseMeshesFile(const std::string& filename,
  bool bfit = false,
  bool boriented = false)
{
  //typedef std::map< unsigned short, std::string > LabelRefMapType;
  typedef std::vector< std::pair< unsigned short, std::string > > LabelRefMapType;
  LabelRefMapType lmap;

  //std::vector< unsigned int > order;

  // parse filename
  std::ifstream ifs(filename.c_str());
  if (!ifs.good() || !ifs.is_open())
  {
    std::cout << "Cannot open meshes file " << filename << std::endl;
    return false;
  }

  std::string line;

  int label;
  char fname[1024];

  while (std::getline(ifs, line))
  {
    // process line
    sscanf(line.c_str(), "%s %d", fname, &label);

    if (label == 0)
    {
      std::cout << fname << " mesh has a forbidden label: ";
      std::cout << "Label 0 is reserved for NO LABEL\n";
      continue;
    }

    lmap.push_back(std::pair< unsigned short, std::string >(label, fname));
  }

  ifs.close();

  // traverse map and perform conversion
  // respect order in file so user can play
  // with it
  if (bfit)
  {
    // image will be created from meshes bboxes
    if (!boriented)
    {
      std::cout << "Computing global bbox for all meshes...\n";
      bool bUpdate = false;
      for (LabelRefMapType::const_iterator ito = lmap.begin(); ito != lmap.end(); ++ito)
      {
        std::cout << "  bbox of mesh " << ito->second << " processed." << std::endl;
        bool b = computeOrUpdatePolyDataBboxFromFile(ito->second, spacing, origin, size, padding, bUpdate);
        if (!b) std::cout << "Error detected, skipping mesh\n";
        if (!bUpdate && b) bUpdate = true;
      }
      std::cout << "Done.\n";
    }
    else
    {
      std::cout << "Computing global oriented bbox for all meshes...\n";
      bool bUpdate = false;
      for (LabelRefMapType::const_iterator ito = lmap.begin(); ito != lmap.end(); ++ito)
      {
        std::cout << "  oriented bbox of mesh " << ito->second << " processed." << std::endl;
        bool b = computeOrUpdatePolyDataFittedOrientedBboxFromFile(ito->second, spacing, origin, size, glob_model2obbox_M, glob_obbox2model_M, padding, bUpdate);
        if (!b) std::cout << "Error detected, skipping mesh\n";
        if (!bUpdate && b) bUpdate = true;
      }
      std::cout << "Done.\n";
    }

    // create new image    
    oImage->SetOrigin(origin);
    oImage->SetDimensions(size);
    oImage->SetSpacing(spacing);
    if (bUchar)
      // char for less than 255 labels
      oImage->AllocateScalars(VTK_UNSIGNED_CHAR, 1);
    else
      // short for more than 255 labels
      oImage->AllocateScalars(VTK_UNSIGNED_SHORT, 1);

    oImage->GetOrigin(origin);
    oImage->GetSpacing(spacing);
    oImage->GetDimensions(size);

    // reset output image
    if (!bUchar)
    {
      unsigned short* optr;
      optr = (unsigned short*)oImage->GetScalarPointer();
      for (unsigned int k = 0; k < oImage->GetNumberOfPoints(); ++k, ++optr) *optr = 0;
    }
    else
    {
      unsigned char* optr = (unsigned char*)oImage->GetScalarPointer();
      for (unsigned int k = 0; k < oImage->GetNumberOfPoints(); ++k, ++optr) *optr = 0;
    }

    std::cout << "dimensions: " << size[0] << ", " << size[1] << ", " << size[2] << std::endl;
    std::cout << "origin    : " << origin[0] << ", " << origin[1] << ", " << origin[2] << std::endl;
    std::cout << "spacing   : " << spacing[0] << ", " << spacing[1] << ", " << spacing[2] << std::endl;
  }

  bool bOk = false;

  // Perform the labeling
  for (LabelRefMapType::const_iterator ito = lmap.begin(); ito != lmap.end(); ++ito)
  {
    bOk = meshToLabel(ito->second, ito->first, glob_model2obbox_M);
    std::cout << "Done.\n";
  }


  if (bOk)
  {
    return true;
  }
  return false;
}

int main(int argc, char** argv)
{

  // Wrap everything in a try block.  Do this every time, 
  // because exceptions will be thrown for problems.
  try {

    TCLAP::CmdLine cmd(
      "Triangular meshes to label images conversion. "
      "For multiple conversions, a config file "
      "has to be used. Each line of the file has the following format: "
      "MeshFilename labelNumber\n"
      "The fit mode (--fit option) will automatically create the output image, "
      "so that it fits the mesh(es) like a bounding box. "
      "The bounding box can be oriented (--ofit option), this implies that "
      "2 files will be output containing the (inverse) transformation from "
      "the model space to the oriented bounding box space in which "
      "the label image is defined. \n"
      "Output transform files will be named as <NAME>_m2b.txt and <NAME>_b2m.txt, "
      "where NAME is either the config file or the name of the mesh to label.\n"
      "In all cases an output image need spacing, origin and size information.\n"
      "Depending on the options, this image information can be provided by:\n"
      "-the user: --origin, --spacing and --size options\n"
      "-a reference image: --ref option, supporting oriented images\n"
      "-the bounding box fittting method which computes size and origin but NOT the spacing\n"
      "Image information can be specified in combined ways knowing that:\n"
      "-image update superseeds user, fit and ref info,\n"
      "-fit info superseeds user and ref info\n"
      "-user info superseeds ref info\n"
      "For example if you want to use the spacing of a reference image but fit a bounding box around the mesh "
      "use the --ref option with the -f option.\n"
      "In fit mode, you can specify a padding with --padding option, padding being specified in physical units.\n"
      "meshes are drawn according to their order in the config file.\n"
      "Meshes supported are vtk polydata, stl, obj and ply.\n"
      "WARNING: Make sure that your meshes are fully triangular and closed. The app "
      "will detect problematic triangular meshes but will not perform the triangularization for you.\n\n"
      "Any 3D itk format is supported for the output image.\n"
      "Author: Jerome Schmid - HEDS - 2012-2021"
      , ' ', "0.2.6");

    // Define a switch and add it to the command line.
    TCLAP::SwitchArg vtkNoStencilArg("", "noVtkStencil", "do not use vtk stencil", cmd, false);

    TCLAP::SwitchArg obbtreeArg("", "useOBBTree", "use vtkOBBTree to compute oriented bounding box", cmd, false);
    TCLAP::SwitchArg ucharArg("", "uchar", "use 255 max labels (unsigned char output image)", cmd, false);
    TCLAP::SwitchArg ofitArg("", "ofit", "fit label image in mesh *oriented* bbox", cmd, false);
    TCLAP::SwitchArg fitArg("f", "fit", "fit label image in mesh bbox", cmd, false);
    TCLAP::SwitchArg multiArg("m", "multi", "multiple conversion", cmd, false);
    TCLAP::SwitchArg updateArg("u", "update", "update output image", cmd, false);

    TCLAP::ValueArg<double> paddingArg(
      "",
      "padding",
      "padding",
      false,
      0.0,
      "padding in case of fitting with bbox",
      cmd
    );

    TCLAP::ValueArg<std::string> refImageArg(
      "r",
      "refImage",
      "refImage",
      false,
      "",
      "Reference image used to get spacing, origin and dimensions",
      cmd
    );

    TCLAP::ValueArg<unsigned int> labelArg(
      "l",
      "label",
      "label",
      false,
      0,
      "label value of the rasterized mesh",
      cmd
    );

    TCLAP::ValueArg<std::string> spacingArg(
      "",
      "spacing",
      "spacing",
      false,
      "",
      "spacing of output image",
      cmd
    );

    TCLAP::ValueArg<std::string> sizeArg(
      "",
      "size",
      "size",
      false,
      "",
      "size of output image",
      cmd
    );

    TCLAP::ValueArg<std::string> originArg(
      "",
      "origin",
      "origin",
      false,
      "",
      "origin of output image",
      cmd
    );

    TCLAP::ValueArg<std::string> outputFileArg(
      "o",
      "out",
      "out file",
      true,
      "",
      "Output label image file", cmd);

    TCLAP::ValueArg<std::string> inputFileArg(
      "i",
      "in",
      "in file",
      true,
      "",
      "Input triangular mesh file or mutiple meshes config file", cmd);

    // Parse the argv array.
    cmd.parse(argc, argv);

    // Get the value parsed by each arg. 
    std::string ifile = inputFileArg.getValue();
    std::string ofile = outputFileArg.getValue();
    bool bUpdate = updateArg.getValue();
    bool bMulti = multiArg.getValue();
    bool bfit = fitArg.getValue();
    bool boriented = ofitArg.getValue();
    if (boriented) bfit = true;
    bUchar = ucharArg.getValue();
    std::string rfile = refImageArg.getValue();
    bool bRef = !rfile.empty() && refImageArg.isSet();
    padding = paddingArg.getValue();
    bForceUseOfOBBTree = obbtreeArg.getValue();
    bUseVtkStencil = !vtkNoStencilArg.isSet();

    if (bForceUseOfOBBTree && !boriented)
    {
      std::cout << "Option --useOBBTree can only be used with oriented bbox option --ofit.\n";
      bForceUseOfOBBTree = false;
    }

    if (bfit && bUpdate)
    {
      std::cout << "Warning: fitting and update modes are exclusive. Fit mode chosen.\n";
      bUpdate = false;
    }

    if (bUpdate && !bRef)
    {
      std::cout << "Update mode: you must specify an image to update with the --ref option";
      exit(EXIT_FAILURE);
    }

    if (bRef && bUpdate)
    {
      bRef = false;
    }

    std::string originStr = originArg.getValue();
    std::string sizeStr = sizeArg.getValue();
    std::string spacingStr = spacingArg.getValue();
    unsigned short label = labelArg.getValue();

    bool bUserSpecifiedSpacing = spacingArg.isSet();
    bool bUserSpecifiedSize = sizeArg.isSet();
    bool bUserSpecifiedOrigin = originArg.isSet();

    bool bSpacingSet = false;
    bool bSizeSet = false;
    bool bOriginSet = false;

    bool bAllocateImage = false;

    directionx[0] = 1; directionx[1] = 0; directionx[2] = 0;
    directiony[0] = 0; directiony[1] = 1; directiony[2] = 0;
    directionz[0] = 0; directionz[1] = 0; directionz[2] = 1;

    Identity(glob_model2obbox_M);
    Identity(glob_obbox2model_M);

    if (bRef)
    {
      // Read image info without loading it completely in memory
      std::cout << "Using reference image: " << rfile << "\n";

      itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(rfile.c_str(), itk::ImageIOFactory::ReadMode);

      if (!imageIO)
      {
        std::string errMsg;
        errMsg = "ImageIOFactory cannot read header of file : ";
        errMsg.append(rfile);
        std::cerr << errMsg << std::endl;
        return 1;
      }

      imageIO->SetFileName(rfile);
      imageIO->ReadImageInformation();

      /** Get the component type, number of components, dimension and pixel type. */
      unsigned int inputDimension = imageIO->GetNumberOfDimensions();
      unsigned int numberOfComponents = imageIO->GetNumberOfComponents();
      std::string inputPixelComponentType = imageIO->GetComponentTypeAsString(
        imageIO->GetComponentType());
      std::string pixelType = imageIO->GetPixelTypeAsString(
        imageIO->GetPixelType());

      if (inputDimension != 3)
      {
        std::cerr << "Image dimension is not equal to 3! Only volumetric images are handled!\n";
        return 1;
      }

      for (unsigned int i = 0; i < inputDimension; ++i)
      {
        origin[i] = imageIO->GetOrigin(i);
        spacing[i] = imageIO->GetSpacing(i);
        size[i] = imageIO->GetDimensions(i);
        directionx[i] = imageIO->GetDirection(0)[i];
        directiony[i] = imageIO->GetDirection(1)[i];
        directionz[i] = imageIO->GetDirection(2)[i];
      }

      glob_obbox2model_M[0] = directionx[0]; glob_obbox2model_M[1] = directiony[0]; glob_obbox2model_M[2] = directionz[0]; glob_obbox2model_M[3] = origin[0];
      glob_obbox2model_M[4] = directionx[1]; glob_obbox2model_M[5] = directiony[1]; glob_obbox2model_M[6] = directionz[1]; glob_obbox2model_M[7] = origin[1];
      glob_obbox2model_M[8] = directionx[2]; glob_obbox2model_M[9] = directiony[2]; glob_obbox2model_M[10] = directionz[2]; glob_obbox2model_M[11] = origin[2];
      glob_obbox2model_M[12] = 0; glob_obbox2model_M[13] = 0; glob_obbox2model_M[14] = 0; glob_obbox2model_M[15] = 1;

      Invert_M<double>(glob_obbox2model_M, glob_model2obbox_M);
      glob_model2obbox_M[3] += origin[0];
      glob_model2obbox_M[7] += origin[1];
      glob_model2obbox_M[11] += origin[2];

      Invert_M<double>(glob_model2obbox_M, glob_obbox2model_M);

      bAllocateImage = true;
      bSpacingSet = true;
      bSizeSet = true;
      bOriginSet = true;

    } // bRef


    if (bUserSpecifiedSpacing || bUserSpecifiedSize || bUserSpecifiedOrigin)
    {
      bool bExit = false;

      double *ox = origin;
      double *oy = origin + 1;
      double *oz = origin + 2;

      double *sx = spacing;
      double *sy = spacing + 1;
      double *sz = spacing + 2;

      int *dx = size;
      int *dy = size + 1;
      int *dz = size + 2;

      // analyse str
      if (bUserSpecifiedOrigin)
      {
        int res = sscanf(originStr.c_str(), "%lf,%lf,%lf", ox, oy, oz);
        if (res != 3)
        {
          std::cerr << "New mode: origin specified has incorrect format.\n";
          bExit = true;
        }
        bOriginSet = true;
      }

      // analyse str
      if (bUserSpecifiedSpacing)
      {
        int res = sscanf(spacingStr.c_str(), "%lf,%lf,%lf", sx, sy, sz);
        if (res != 3)
        {
          std::cerr << "New mode: spacing specified has incorrect format.\n";
          bExit = true;
        }
        bSpacingSet = true;
      }

      // analyse str
      if (bUserSpecifiedSize)
      {
        int res = sscanf(sizeStr.c_str(), "%d,%d,%d", dx, dy, dz);
        if (res != 3)
        {
          std::cerr << "New mode: size specified has incorrect format.\n";
          bExit = true;
        }
        bSizeSet = true;
      }

      if (bExit) return 1;

      // create new image with specification given by cmd line
      bAllocateImage = true;
    }

    if (bUpdate && !bfit)
    {
      if (bUchar)
      {
        typedef itk::Image< unsigned char, 3 >   ImageType;
        typedef itk::ImageFileReader< ImageType >    ReaderType;
        typedef itk::ImageIOBase                     ImageIOBaseType;

        /** Create a reader. */
        ReaderType::Pointer reader = ReaderType::New();
        reader->SetFileName(rfile);

        /** Generate all information. */
        try
        {
          reader->Update();
        }
        catch (itk::ExceptionObject  &  err)
        {
          std::cerr << "Update mode: cannot read file" << std::endl;
          std::cerr << err << std::endl;
          return 1;
        }

        typedef itk::ImageToVTKImageFilter< ImageType > ConnectorType;
        ConnectorType::Pointer connector = ConnectorType::New();
        connector->SetInput(reader->GetOutput());
        try
        {
          connector->Update();
        }
        catch (itk::ExceptionObject  &  err)
        {
          std::cerr << "Update mode: cannot read file" << std::endl;
          std::cerr << err << std::endl;
          return 1;
        }
        oImage = vtkStructuredPoints::New();
        oImage->DeepCopy(connector->GetOutput());
      }
      else
      {
        typedef itk::Image< unsigned short, 3 >   ImageType;
        typedef itk::ImageFileReader< ImageType >    ReaderType;
        typedef itk::ImageIOBase                     ImageIOBaseType;

        /** Create a reader. */
        ReaderType::Pointer reader = ReaderType::New();
        reader->SetFileName(rfile);

        /** Generate all information. */
        try
        {
          reader->Update();
        }
        catch (itk::ExceptionObject  &  err)
        {
          std::cerr << "Update mode: cannot read file" << std::endl;
          std::cerr << err << std::endl;
          return 1;
        }

        typedef itk::ImageToVTKImageFilter< ImageType > ConnectorType;
        ConnectorType::Pointer connector = ConnectorType::New();
        connector->SetInput(reader->GetOutput());
        try
        {
          connector->Update();
        }
        catch (itk::ExceptionObject  &  err)
        {
          std::cerr << "Update mode: cannot read file" << std::endl;
          std::cerr << err << std::endl;
          return 1;
        }
        oImage = vtkStructuredPoints::New();
        oImage->DeepCopy(connector->GetOutput());
      }

      bAllocateImage = false;
      bSizeSet = true;
      bOriginSet = true;
      bSpacingSet = true;
    }
    
    if (paddingArg.isSet() && !bfit)
    {
      std::cout << "Padding specified, but it is only valid in fitting mode, option ignored\n";
    }
    else
    {
      std::cout << "Padding of " << padding << " (physical units) specified for fitting mode\n";
    }

    if (bfit && !bMulti)
    {
      double *sx = spacing;
      double *sy = spacing + 1;
      double *sz = spacing + 2;

      bool bExit = false;

      if (!bSpacingSet)
      {
        std::cerr << "Fit mode: no spacing specified.\n";
        bExit = true;
      }

      if (bExit) exit(EXIT_FAILURE);

      bool bTest = false;

      if (boriented)
        bTest = computeOrUpdatePolyDataFittedOrientedBboxFromFile(ifile, spacing, origin, size, glob_model2obbox_M, glob_obbox2model_M, padding, false);
      else
        bTest = computeOrUpdatePolyDataBboxFromFile(ifile, spacing, origin, size, padding, false);

      if (!bTest) exit(EXIT_FAILURE);

      bAllocateImage = true;
      bSizeSet = true;
      bOriginSet = true;
      bSpacingSet = true;
    }

    if (bfit && bMulti)
    {
      double *sx = spacing;
      double *sy = spacing + 1;
      double *sz = spacing + 2;

      bool bExit = false;

      if (!bSpacingSet)
      {
        std::cerr << "Fit mode: no spacing specified.\n";
        bExit = true;
      }

      if (bExit) exit(EXIT_FAILURE);

      oImage = vtkStructuredPoints::New();
      bAllocateImage = false;
      bSizeSet = true;
      bOriginSet = true;
      bSpacingSet = true;
    }

    if (bAllocateImage)
    {
      if (!bSizeSet) {
        std::cerr << "output image allocation: size has not been specified!\n";
        exit(EXIT_FAILURE);;
      }

      if (!bSpacingSet) {
        std::cerr << "output image allocation: spacing has not been specified!\n";
        exit(EXIT_FAILURE);;
      }

      if (!bOriginSet) {
        std::cerr << "output image allocation: origin has not been specified!\n";
        exit(EXIT_FAILURE);;
      }
      // create new image with specification given by cmd line
      oImage = vtkStructuredPoints::New();
      oImage->SetOrigin(origin);
      oImage->SetDimensions(size);
      oImage->SetSpacing(spacing);
      if (bUchar)
        oImage->AllocateScalars(VTK_UNSIGNED_CHAR, 1);
      else
        oImage->AllocateScalars(VTK_UNSIGNED_SHORT, 1);
    }

    if (oImage == NULL && !bfit && !bRef)
    {
      std::cerr << "Update mode: cannot read image to update.\n";
      exit(EXIT_FAILURE);
    }
    else if (oImage != NULL && !bfit)
    {
      oImage->GetOrigin(origin);
      oImage->GetSpacing(spacing);
      oImage->GetDimensions(size);

      std::cout << "out image: " << ofile << std::endl;
      std::cout << "dimensions: " << size[0] << ", " << size[1] << ", " << size[2] << std::endl;
      std::cout << "origin    : " << origin[0] << ", " << origin[1] << ", " << origin[2] << std::endl;
      std::cout << "spacing   : " << spacing[0] << ", " << spacing[1] << ", " << spacing[2] << std::endl;
      std::cout << "directionx: " << directionx[0] << ", " << directionx[1] << ", " << directionx[2] << std::endl;
      std::cout << "directiony: " << directiony[0] << ", " << directiony[1] << ", " << directiony[2] << std::endl;
      std::cout << "directionz: " << directionz[0] << ", " << directionz[1] << ", " << directionz[2] << std::endl;
      bSizeSet = true;
      bOriginSet = true;
      bSpacingSet = true;
    }

    std::cout << "mode is ";
    if (bUpdate) std::cout << "UPDATE: the output image is read from file and labels are updated.\n";
    else
    {
      // set 0 value to new image to be created
      if (!bUchar && !bMulti)
      {
        unsigned short* optr;
        optr = (unsigned short*)oImage->GetScalarPointer();
        for (unsigned int k = 0; k < oImage->GetNumberOfPoints(); ++k, ++optr) *optr = 0;
      }

      if (!bUchar && !bMulti)
      {
        unsigned char* optr = (unsigned char*)oImage->GetScalarPointer();
        for (unsigned int k = 0; k < oImage->GetNumberOfPoints(); ++k, ++optr) *optr = 0;
      }

      if (bfit && !boriented) std::cout << "NEW+FIT: the output image is created according to a bbox around the mesh.\n";
      else if (bfit && boriented) std::cout << "NEW+ORIENTED FIT: the output image is created according to an oriented bbox around the mesh.\n";
      else std::cout << "NEW: the output image is created according to the user image parameters.\n";
    }

    bool bOk = false;
    if (bMulti) // multi conversion with or without fit
    {
      std::cout << "Meshfile " << ifile << " is parsed for multiple conversion. \n\n";
      bOk = analyseMeshesFile(ifile, bfit, boriented);
    }
    else // single conversion with or without fit
    {
      // convert mesh to label
      bOk = meshToLabel(ifile, label, glob_model2obbox_M);
    }

    if (boriented)
    {
      // Saving transform to files for possible future usage
      std::string filename = ifile + std::string("_m2b.txt");
      SaveHomogeneousMatrixToFile(filename, glob_model2obbox_M);
      filename = ifile + std::string("_b2m.txt");
      SaveHomogeneousMatrixToFile(filename, glob_obbox2model_M);
    }

    if (bOk) writeImageToFile(oImage, ofile);
    else std::cout << "Label image has not been created.\n";

    if (oImage) SafeDelete(oImage);
  }
  catch (TCLAP::ArgException &e)  // catch any exceptions
  {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    exit(EXIT_FAILURE);
  }
  catch (...)
  {
    std::cerr << "vtk or system error\n";
    exit(EXIT_FAILURE);
  }

  return EXIT_SUCCESS;
}