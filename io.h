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

#ifndef __HEDS_IO_H__
#define __HEDS_IO_H__

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <string>


class ErrorObserver{
public:
    void ErrorEventObserver(
        vtkObject* caller,
        long unsigned int eventId,
        void* callData ){ throw 1; }
};

// Description:
// Check vtkUnstructuredGrid if has tetrahedra cells
extern bool UnstructuredGridHasTets(vtkSmartPointer< vtkUnstructuredGrid > ugrid);

// Description:
// Extract surface from vtkUnstructuredGrid
extern bool ExtractSurfaceFromUnstructuredGrid(vtkSmartPointer< vtkUnstructuredGrid > ugrid, vtkSmartPointer< vtkPolyData > poly, bool withCerr);


// Description:
// Load an unstructured grid from a specified file (.vtk).
extern int LoadUnstructuredGridFromVTK(  
    std::string VTKFileName,
    vtkSmartPointer< vtkUnstructuredGrid > ReadUnstructuredGrid,
        bool withCerr = false);

// Description:
// Load a polydata from a specified file (.vtk).
extern int LoadPolyDataFromVTK(  
    std::string VTKFileName,
    vtkSmartPointer< vtkPolyData > ReadPolyData,
        bool withCerr = false);

// Description:
// Load a vtkPolyData from an STL file.
extern int LoadPolyDataFromSTL(
    std::string STLFileName,
    vtkSmartPointer< vtkPolyData > readPolyData,
    bool withCerr = false);

// Description:
// Load a vtkPolyData from an OBJ file.
extern int LoadPolyDataFromOBJ(
    std::string OBJFileName,
    vtkSmartPointer< vtkPolyData > readPolyData,
    bool withCerr = false);

// Description:
// Load a vtkPolyData from a PLY file.
extern int LoadPolyDataFromPLY(
    std::string PLYFileName,
    vtkSmartPointer< vtkPolyData > readPolyData,
    bool withCerr = false);


// Description:
// Load a polydata from a file.
extern int LoadPolyData(  
    std::string fileName,
    vtkSmartPointer< vtkPolyData > ReadPolyData,
    bool bUseExtensionOnly = true );

// Description:
// Save an unstructured grid to a specified file (.vtk).
extern int SaveUnstrcturedGridToVTK( 
    vtkSmartPointer< vtkUnstructuredGrid > UGridToSave, 
    std::string VTKFileName,
    bool ascii = false,
    bool withCerr = false);

// Description:
// Save a polydata to a specified file (.vtk).
extern int SavePolyDataToVTK( 
    vtkSmartPointer< vtkPolyData > PolyDataToSave, 
    std::string VTKFileName,
    bool ascii = false,
    bool withCerr = false);

// Description:
// Save a vtkPolyData to an STL file.
extern int SavePolyDataToSTL(
    vtkSmartPointer< vtkPolyData > polyDataToSave,
    std::string STLFileName,
    bool ascii = false,
    bool withCerr = false);

// Description:
// Save a vtkPolyData to an OBJ file.
extern int SavePolyDataToOBJ(
    vtkSmartPointer< vtkPolyData > polyDataToSave,
    std::string OBJFileName,
    bool withCerr = false);

// Description:
// Save a vtkPolyData to an PLY file.
int SavePolyDataToPLY(
    vtkSmartPointer< vtkPolyData > polyDataToSave,
    std::string PLYFileName,
    bool ascii = false,
    bool withCerr = false);

// Description:
// Save a polydata to a specified file (.vtk | .stl).
extern int SavePolyData( 
    vtkSmartPointer< vtkPolyData > PolyDataToSave, 
    std::string FileName,
    bool ascii = false, 
    bool withCerr = false );


#endif