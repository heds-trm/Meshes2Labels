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

#include <vtkSmartPointer.h>
#include <vtkPolyDataWriter.h>
#include <vtkSTLWriter.h>
#include <vtkPLYReader.h>
#include <vtkPLYWriter.h>
#include <vtkPolyData.h>
#include <vtkSTLReader.h>
#include <vtkPolyDataReader.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkOBJReader.h>
#include "vtkOBJWriter.h"
#include <vtkCommand.h>
#include <vtkCellIterator.h>
#include <algorithm>
#include <vtkDataSetSurfaceFilter.h>
#include "io.h"


bool UnstructuredGridHasTets(vtkSmartPointer< vtkUnstructuredGrid > ugrid)
{
  vtkCellIterator *it = ugrid->NewCellIterator();
  for (it->InitTraversal(); !it->IsDoneWithTraversal(); it->GoToNextCell())
  {
    if (it->GetCellType() == VTK_TETRA)
    {
      return true; 
    }
  }

  return false;
}

bool ExtractSurfaceFromUnstructuredGrid(
  vtkSmartPointer< vtkUnstructuredGrid > ugrid, 
  vtkSmartPointer< vtkPolyData > poly, bool withCerr)
{
  auto surfaceFilter = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
  surfaceFilter->SetInputData(ugrid);

  ErrorObserver VTKReaderErrorObserver;
  surfaceFilter->AddObserver(
    vtkCommand::ErrorEvent,
    &VTKReaderErrorObserver,
    &ErrorObserver::ErrorEventObserver);
  try {
    surfaceFilter->Update();
  }
  catch (...) {
    if (withCerr)
      cerr << "ERROR in ExtractSurfaceFromUnstructuredGrid";
    return false; // Failure
  }

  poly->DeepCopy(surfaceFilter->GetOutput());

  return true;
}

// Description:
// Load an unstructured grid from a specified file (.vtk).
int LoadUnstructuredGridFromVTK(  
    std::string VTKFileName,
    vtkSmartPointer< vtkUnstructuredGrid > ReadUnstructuredGrid,
        bool withCerr ){

   vtkSmartPointer< vtkUnstructuredGridReader > UnstructuredGridReader = 
        vtkSmartPointer< vtkUnstructuredGridReader >::New();
    UnstructuredGridReader->SetFileName( VTKFileName.c_str() );

    ErrorObserver VTKReaderErrorObserver;
    UnstructuredGridReader->AddObserver( 
        vtkCommand::ErrorEvent, 
        &VTKReaderErrorObserver, 
        &ErrorObserver::ErrorEventObserver );
    try{
        UnstructuredGridReader->Update();
    }catch(...){
        if (withCerr)
            cerr << "ERROR in LoadUnstructuredGrid - "
                << " Failed to load UnstructuredGrid from VTK file:  " 
                << VTKFileName <<"\n";
        return 1; // Failure
    }

    ReadUnstructuredGrid->DeepCopy( UnstructuredGridReader->GetOutput() );

    return 0; // Success
}

// Description:
// Load a polydata from a specified file (.vtk).
int LoadPolyDataFromVTK(  
    std::string VTKFileName,
    vtkSmartPointer< vtkPolyData > ReadPolyData,
        bool withCerr )
{
    vtkSmartPointer< vtkPolyDataReader > PolyDataReader = 
        vtkSmartPointer< vtkPolyDataReader >::New();
    PolyDataReader->SetFileName( VTKFileName.c_str() );

    ErrorObserver VTKReaderErrorObserver;
    PolyDataReader->AddObserver( 
        vtkCommand::ErrorEvent, 
        &VTKReaderErrorObserver, 
        &ErrorObserver::ErrorEventObserver );
    try{
        PolyDataReader->Update();
    }catch(...){
        if (withCerr)
            cerr << "ERROR in LoadPolyData - "
                << " Failed to load polydata from VTK file:  " 
                << VTKFileName <<"\n";
        return 1; // Failure
    }

    ReadPolyData->DeepCopy( PolyDataReader->GetOutput() );

    return 0; // Success
}

// Description:
// Load a vtkPolyData from an STL file.
int LoadPolyDataFromSTL(
    std::string STLFileName,
    vtkSmartPointer< vtkPolyData > readPolyData,
    bool withCerr )
{
    vtkSmartPointer< vtkSTLReader > STLReader = 
        vtkSmartPointer< vtkSTLReader >::New();
    STLReader->SetFileName( STLFileName.c_str() );

    ErrorObserver STLReaderErrorObserver;
    STLReader->AddObserver( 
        vtkCommand::ErrorEvent, 
        &STLReaderErrorObserver, 
        &ErrorObserver::ErrorEventObserver );
    try {
        STLReader->Update();
    } catch (...)
    {
        if (withCerr)
            cerr << "ERROR in LoadPolyDataFromSTL - "
                 << "Failed to load the polydata from the STL file: " 
                 << STLFileName << ".\n";
        return 1; // Failure
    }

    readPolyData->DeepCopy( STLReader->GetOutput() );

    return 0; // Success
}

// Description:
// Load a vtkPolyData from an OBJ file.
int LoadPolyDataFromOBJ(
    std::string OBJFileName,
    vtkSmartPointer< vtkPolyData > readPolyData,
    bool withCerr )
{
    vtkSmartPointer< vtkOBJReader > OBJReader = 
        vtkSmartPointer< vtkOBJReader >::New();
    OBJReader->SetFileName( OBJFileName.c_str() );

    ErrorObserver OBJReaderErrorObserver;
    OBJReader->AddObserver( 
        vtkCommand::ErrorEvent, 
        &OBJReaderErrorObserver, 
        &ErrorObserver::ErrorEventObserver );
    try {
        OBJReader->Update();
    } catch (...)
    {
        if (withCerr)
            cerr << "ERROR in LoadPolyDataFromOBJ - "
                 << "Failed to load the polydata from the OBJ file: " 
                 << OBJFileName << ".\n";
        return 1; // Failure
    }

    readPolyData->DeepCopy( OBJReader->GetOutput() );

    return 0; // Success
}

// Description:
// Load a vtkPolyData from a PLY file.
int LoadPolyDataFromPLY(
    std::string PLYFileName,
    vtkSmartPointer< vtkPolyData > readPolyData,
    bool withCerr )
{
    vtkSmartPointer< vtkPLYReader > PLYReader = 
        vtkSmartPointer< vtkPLYReader >::New();
    PLYReader->SetFileName( PLYFileName.c_str() );

    ErrorObserver PLYReaderErrorObserver;
    PLYReader->AddObserver( 
        vtkCommand::ErrorEvent, 
        &PLYReaderErrorObserver, 
        &ErrorObserver::ErrorEventObserver );
    try {
        PLYReader->Update();
    } catch (...)
    {
        if (withCerr)
            cerr << "ERROR in LoadPolyDataFromPLY - "
                 << "Failed to load the polydata from the PLY file: " 
                 << PLYFileName << ".\n";
        return 1; // Failure
    }

    readPolyData->DeepCopy( PLYReader->GetOutput() );

    return 0; // Success
}


inline bool ends_with(std::string const & value, std::string const & ending)
{
    if (ending.size() > value.size()) return false;
    return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

// Description:
// Load a polydata from a file.
int LoadPolyData(  
    std::string fileName,
    vtkSmartPointer< vtkPolyData > ReadPolyData,
    bool bUseExtensionOnly )
{
    if( bUseExtensionOnly )
    {
        int res = 1;
        std::string lowerstr = fileName;
        std::transform(lowerstr.begin(), lowerstr.end(), lowerstr.begin(), ::tolower);

        if( ends_with(lowerstr, ".vtk") )
        {
            res = LoadPolyDataFromVTK( fileName, ReadPolyData, false );
        }
        else if( ends_with(lowerstr, ".stl") )
        {
            res = LoadPolyDataFromSTL( fileName, ReadPolyData, false );
        }
        else if( ends_with(lowerstr, ".obj") )
        {
            res = LoadPolyDataFromOBJ( fileName, ReadPolyData, false );
        }
        else if( ends_with(lowerstr, ".ply") )
        {
            res = LoadPolyDataFromPLY( fileName, ReadPolyData, false );
        }

        if (res == 1)
        {  
            cerr << "ERROR in LoadPolyData - cannot load polydata from file: " << fileName << endl;
            return 1;
        }
    }
    else
    {
        int res = 0;
        res = LoadPolyDataFromVTK( fileName, ReadPolyData, false );
        if (res == 1)
            res = LoadPolyDataFromSTL( fileName, ReadPolyData, false );
        if( res == 1 )
            res = LoadPolyDataFromOBJ( fileName, ReadPolyData, false );
        if( res == 1 )
            res = LoadPolyDataFromPLY( fileName, ReadPolyData, false );
        if (res == 1)
        {  
            cerr << "ERROR in LoadPolyData - cannot load polydata from file: " << fileName << endl;
            return 1;
        }
    }

    return 0;
}

int SaveUnstrcturedGridToVTK( 
    vtkSmartPointer< vtkUnstructuredGrid > UGridToSave, 
    std::string VTKFileName,
    bool ascii,
    bool withCerr )
{
    vtkSmartPointer< vtkUnstructuredGridWriter > UnstructuredGridWriter = 
        vtkSmartPointer< vtkUnstructuredGridWriter >::New();
    UnstructuredGridWriter->SetFileName( VTKFileName.c_str() );
    UnstructuredGridWriter->SetInputData( UGridToSave );
    if( ascii ) UnstructuredGridWriter->SetFileTypeToASCII();
    else UnstructuredGridWriter->SetFileTypeToBinary();

    ErrorObserver UnstructuredGridWriterErrorObserver;
    UnstructuredGridWriter->AddObserver( 
        vtkCommand::ErrorEvent, 
        &UnstructuredGridWriterErrorObserver, 
        &ErrorObserver::ErrorEventObserver );
    try{
        UnstructuredGridWriter->Update();
    }catch(...){
        if (withCerr) 
        {
            cerr << "ERROR in SaveUnstructuredGrid - " 
                << "Failed to save UnstructuredGrid to VTK file: "<< VTKFileName <<"\n";
            return 1; // Failure
        }
    }

    return 0; // Success
}

// Description:
// Save a polydata to a specified file (.vtk).
int SavePolyDataToVTK( 
    vtkSmartPointer< vtkPolyData > PolyDataToSave, 
    std::string VTKFileName,
    bool ascii,
    bool withCerr )
{
    vtkSmartPointer< vtkPolyDataWriter > PolyDataWriter = 
        vtkSmartPointer< vtkPolyDataWriter >::New();
    PolyDataWriter->SetFileName( VTKFileName.c_str() );
    PolyDataWriter->SetInputData( PolyDataToSave );
    if( ascii ) PolyDataWriter->SetFileTypeToASCII();
    else PolyDataWriter->SetFileTypeToBinary();

    ErrorObserver PolyDataWriterErrorObserver;
    PolyDataWriter->AddObserver( 
        vtkCommand::ErrorEvent, 
        &PolyDataWriterErrorObserver, 
        &ErrorObserver::ErrorEventObserver );
    try{
        PolyDataWriter->Update();
    }catch(...){
        if (withCerr) 
        {
            cerr << "ERROR in SavePolyData - " 
                << "Failed to save polydata to VTK file: "<< VTKFileName <<"\n";
            return 1; // Failure
        }
    }

    return 0; // Success
}

// Description:
// Save a polydata to a specified file (.OBJ).
int SavePolyDataToOBJ(
    vtkSmartPointer< vtkPolyData > polyDataToSave,
    std::string OBJFileName,
    bool withCerr )
{
    vtkSmartPointer< vtkOBJWriter > PolyDataWriter = 
        vtkSmartPointer< vtkOBJWriter >::New();
    PolyDataWriter->SetFileName( OBJFileName.c_str() );
    PolyDataWriter->SetInputData( polyDataToSave );

    ErrorObserver PolyDataWriterErrorObserver;
    PolyDataWriter->AddObserver( 
        vtkCommand::ErrorEvent, 
        &PolyDataWriterErrorObserver, 
        &ErrorObserver::ErrorEventObserver );
    try{
        PolyDataWriter->Update();
    }catch(...){
        if (withCerr) 
        {
            cerr << "ERROR in SavePolyData - " 
                << "Failed to save polydata to OBJ file: "<< OBJFileName <<"\n";
            return 1; // Failure
        }
    }

    return 0; // Success
}

// Description:
// Save a vtkPolyData to an STL file.
int SavePolyDataToSTL(
    vtkSmartPointer< vtkPolyData > polyDataToSave,
    std::string STLFileName,
    bool ascii,
    bool withCerr )
{
    vtkSmartPointer< vtkSTLWriter > STLWriter = 
        vtkSmartPointer< vtkSTLWriter >::New();
    STLWriter->SetFileName( STLFileName.c_str() );
    STLWriter->SetInputData( polyDataToSave );
    if( ascii ) STLWriter->SetFileTypeToASCII();
    else STLWriter->SetFileTypeToBinary();

    ErrorObserver STLWriterErrorObserver;
    STLWriter->AddObserver( 
        vtkCommand::ErrorEvent, 
        &STLWriterErrorObserver, 
        &ErrorObserver::ErrorEventObserver );
    try {
        STLWriter->Update();
    } catch (...)
    {
        if (withCerr) 
        {
            cerr << "ERROR in SavePolyDataToSTL - "
                 << "Failed to save the polydata to the STL file: " 
                 << STLFileName << ".\n";
            return 1; // Failure
        }
    }

    return 0; // Success
}

// Description:
// Save a vtkPolyData to an PLY file.
int SavePolyDataToPLY(
    vtkSmartPointer< vtkPolyData > polyDataToSave,
    std::string PLYFileName,
    bool ascii,
    bool withCerr )
{
    vtkSmartPointer< vtkPLYWriter > PLYWriter = 
        vtkSmartPointer< vtkPLYWriter >::New();
    PLYWriter->SetFileName( PLYFileName.c_str() );
    PLYWriter->SetInputData( polyDataToSave );
    if( ascii ) PLYWriter->SetFileTypeToASCII();
    else PLYWriter->SetFileTypeToBinary();

    ErrorObserver PLYWriterErrorObserver;
    PLYWriter->AddObserver( 
        vtkCommand::ErrorEvent, 
        &PLYWriterErrorObserver, 
        &ErrorObserver::ErrorEventObserver );
    try {
        PLYWriter->Update();
    } catch (...)
    {
        if (withCerr) 
        {
            cerr << "ERROR in SavePolyDataToPLY - "
                 << "Failed to save the polydata to the PLY file: " 
                 << PLYFileName << ".\n";
            return 1; // Failure
        }
    }

    return 0; // Success
}

// Description:
// Save a polydata to a specified file (.vtk | .stl).
int SavePolyData( 
    vtkSmartPointer< vtkPolyData > PolyDataToSave, 
    std::string FileName,
    bool ascii, 
    bool withCerr )
{    
    std::string ext = FileName.substr(FileName.size()-4, FileName.size());
    if ( strcmp( ext.c_str(), ".stl" ) == 0)
        SavePolyDataToSTL( PolyDataToSave, FileName, ascii, withCerr );
    else if (strcmp( ext.c_str(), ".vtk" ) == 0)
        SavePolyDataToVTK( PolyDataToSave, FileName, ascii, withCerr );
    else if(strcmp( ext.c_str(), ".obj" ) == 0)
        SavePolyDataToOBJ( PolyDataToSave, FileName, withCerr );
    else if(strcmp( ext.c_str(), ".ply" ) == 0)
        SavePolyDataToPLY( PolyDataToSave, FileName, ascii, withCerr );
    else
        SavePolyDataToVTK( PolyDataToSave, FileName+".vtk", ascii, withCerr );

    return 0;
}
