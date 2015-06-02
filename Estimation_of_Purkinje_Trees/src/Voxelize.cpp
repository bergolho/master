/*
* Copyright (c) 2009,
* Computational Image and Simulation Technologies in Biomedicine (CISTIB),
* Universitat Pompeu Fabra (UPF), Barcelona, Spain. All rights reserved.
* See license.txt file for details.
*/

#include <vtkXMLPolyDataReader.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyData.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkImageStencilData.h>
#include <vtkImageStencil.h>
#include <vtkImageData.h>
#include <vtkImageExport.h>
#include <vtkStructuredPointsWriter.h>
#include <string.h>
#include <itkVTKImageImport.h>
#include <itkImage.h>
#include <itkImageFileWriter.h>
#include "vtkITKUtility.h"

int main( int argc, char* argv[] )
{
	double bounds[6]; // Mesh bounding box
	double spacing[3]; // Output image voxel spacing
	double origin[3]; // Output image origin
	int extent[6]; // Output image size in voxels (from..to voxel index, hence the 6 elements)
        int imagesize[3];
	imagesize[0] = -1;
	imagesize[1] = -1;
	imagesize[2] = -1;
        float input_bounds[6];
        input_bounds[0] = -1;
	// Default spacing
	spacing[0] = 1.0;
	spacing[1] = 1.0;
	spacing[2] = 1.0;
      
    if (argc < 3) {
		printf("Author: Ruben Cardenes, Feb. 2011 \n");
		printf("Description: ");
		printf("Usage: %s [options] input_surface.vtk output_vol.[mhd/vtk]\n",argv[0]);
		printf("              -sp X (isotropic)\n");
		printf("              -spXYZ X Y Z\n");
		printf("              -imSize X Y Z \n");
		printf("              -bounds X1 X2 Y1 Y2 Z1 Z2 \n");
		return 1;
    } 

    std::string filename = argv[1];
    argc--;
    argv++;
    std::string filename_out = argv[1];
    argc--;
    argv++;

    int ok = 0;
    // Parse input arguments
    while ( argc > 1 ) {
	ok = false;

	if ( ( ok == false ) && ( strcmp( argv[1], "-sp" ) == 0 ) )
	{
		argc--;
		argv++;
		spacing[0] = atof( argv[1] );
		argc--;
		argv++;
                spacing[1] = spacing[0];
 		spacing[2] = spacing[0];
		ok = true;
	}
 	if ( ( ok == false ) && ( strcmp( argv[1], "-spXYZ" ) == 0 ) )
	{
		argc--;
		argv++;
		spacing[0] = atof( argv[1] );
		argc--;
		argv++;
		spacing[1] = atof( argv[1] );
		argc--;
		argv++;
		spacing[2] = atof( argv[1] );
		argc--;
		argv++;
		ok = true;
	}
        if ( ( ok == false ) && ( strcmp( argv[1], "-imSize" ) == 0 ) )
	{
		argc--;
		argv++;
		imagesize[0] = atoi( argv[1] );
		argc--;
		argv++;
		imagesize[1] = atoi( argv[1] );
		argc--;
		argv++;
		imagesize[2] = atoi( argv[1] );
		argc--;
		argv++;
		ok = true;
	}   
	if ( ( ok == false ) && ( strcmp( argv[1], "-bounds" ) == 0 ) )
	{
		argc--;
		argv++;
		input_bounds[0] = atof( argv[1] );
		argc--;
		argv++;
		input_bounds[1] = atof( argv[1] );
		argc--;
		argv++;
		input_bounds[2] = atof( argv[1] );
		argc--;
		argv++;
		input_bounds[3] = atof( argv[1] );
		argc--;
		argv++;
		input_bounds[4] = atof( argv[1] );
		argc--;
		argv++;
		input_bounds[5] = atof( argv[1] );
		argc--;
		argv++;
		ok = true;
	}   	
	if ( ok == false )
	{
		std::cerr << "Can not parse argument " << argv[1] << std::endl;
		return 1;
	}
   }

	//vtkXMLPolyDataReader* dataReader			= vtkXMLPolyDataReader::New();
	vtkPolyDataReader* dataReader			= vtkPolyDataReader::New();
	vtkPolyData* shape							= vtkPolyData::New();
	vtkPolyDataToImageStencil* stencilFilter	= vtkPolyDataToImageStencil::New();
	vtkImageStencilData* stencilData			= vtkImageStencilData::New();
	vtkImageStencil* imageStencil				= vtkImageStencil::New();
	vtkImageData* outputImage					= vtkImageData::New();
	vtkStructuredPointsWriter* spWriter			= vtkStructuredPointsWriter::New();

	//if( !dataReader->CanReadFile( filename.c_str() ) )
	//{
	//	std::cerr
	//		<< "Error: could not read data file:" << std::endl
	//		<< argv[1] << std::endl
	//		<< "Abnormal program termination." << std::endl;
	//	exit( 1 );
	//}
	dataReader->SetFileName( filename.c_str() );
	dataReader->Update();

	shape = dataReader->GetOutput();
	shape->GetBounds( bounds );

	// Set the output image origin and size in voxels; the size is
	// now by default set to the mesh bounding box plus two full
	// voxels on each side.
        if (imagesize[0] == -1) { 
         if (input_bounds[0] != -1) {
           bounds[0] = input_bounds[0]; bounds[1] = input_bounds[1]; bounds[2] = input_bounds[2]; 
	   bounds[3] = input_bounds[3]; bounds[4] = input_bounds[4]; bounds[5] = input_bounds[5];
         }
	 origin[0] = bounds[0] - 2 * spacing[0];
	 origin[1] = bounds[2] - 2 * spacing[1];
	 origin[2] = bounds[4] - 2 * spacing[2];
	 extent[0] = 0;
	 extent[1] = (bounds[1] - origin[0]) / spacing[0] + 2;
	 extent[2] = 0;
	 extent[3] = (bounds[3] - origin[1]) / spacing[1] + 2;
	 extent[4] = 0;
	 extent[5] = (bounds[5] - origin[2]) / spacing[2] + 2;
	 std::cout << " bounds " << bounds[0] << " " << bounds[1] << " " << bounds[2] << " " << bounds[3] << " " << bounds[4] << " " << bounds[5] << " " << std::endl;
	 std::cout << " extent " << extent[1] << " " << extent[3] << " " << extent[5] << " " << std::endl;
	 std::cout << " origin " << origin[0] << " " << origin[1] << " " << origin[2] << " " << std::endl;
        } else {
         /// If extent (size) is given instead of spacing:
         if (input_bounds[0] != -1) {
           bounds[0] = input_bounds[0]; bounds[1] = input_bounds[1]; bounds[2] = input_bounds[2]; 
	   bounds[3] = input_bounds[3]; bounds[4] = input_bounds[4]; bounds[5] = input_bounds[5];
         }
         spacing[0] = (bounds[1] - bounds[0]) / ( (float)imagesize[0] - 4.0);
         spacing[1] = (bounds[3] - bounds[2]) / ( (float)imagesize[1] - 4.0);
         spacing[2] = (bounds[5] - bounds[4]) / ( (float)imagesize[2] - 4.0);;
         origin[0] = bounds[0] - 2 * spacing[0];
	 origin[1] = bounds[2] - 2 * spacing[1];
	 origin[2] = bounds[4] - 2 * spacing[2];
         extent[0] = 0;
         extent[1] = imagesize[0];
         extent[2] = 0;
         extent[3] = imagesize[1];
         extent[4] = 0;
         extent[5] = imagesize[2];
	 std::cout << " imagesize " << imagesize[0] << " " << imagesize[1] << " " << imagesize[2] << " " << std::endl;
	 std::cout << " spacing " << spacing[0] << " " << spacing[1] << " " << spacing[2] << " " << std::endl;
	 std::cout << " origin " << origin[0] << " " << origin[1] << " " << origin[2] << " " << std::endl;
         ///////////////////////////////////////////77
        }
	stencilFilter->SetInput( shape );
	stencilFilter->SetOutputOrigin( origin );
	stencilFilter->SetOutputSpacing( spacing );
	stencilFilter->SetOutputWholeExtent( extent );
	stencilFilter->Update();

	outputImage->SetOrigin( origin );
	outputImage->SetSpacing( spacing );
	outputImage->SetExtent( extent );
	outputImage->SetScalarTypeToUnsignedChar();
	outputImage->SetNumberOfScalarComponents( 1 );
	outputImage->AllocateScalars();

	imageStencil->SetInput( outputImage );
	imageStencil->SetStencil( stencilFilter->GetOutput() );

	// A problem regarding a difference in behavior of the debug
	// and release modes arose here; the stencilFilter produces
	// 0 values for voxels outside the mesh in debug mode, while
	// it produces 0 values for voxels INSIDE the mesh in release
	// mode.
//	imageStencil->SetBackgroundValue( 0 ); // in debug mode, this produces white meshes on a black background
	imageStencil->SetBackgroundValue( 255 ); // in release mode, this produces black meshes on a white background
	imageStencil->Update();

	//outputImage = imageStencil->GetOutput();

	//spWriter->SetFileTypeToBinary();
	//spWriter->SetInput( outputImage );
	//spWriter->SetFileName( (filename_out).c_str() );
	//spWriter->Write();
	

	typedef itk::Image<unsigned char,3> TOutputImage; 

        typedef itk::VTKImageImport< TOutputImage > ImportType;
	ImportType::Pointer importer = ImportType::New();

	vtkImageExport* exporter = vtkImageExport::New();
	exporter->SetInput(imageStencil->GetOutput());
	ConnectPipelines(exporter,importer);
	

	typedef itk::ImageFileWriter< itk::Image<unsigned char,3> > writerType; 
	writerType::Pointer itkwriter = writerType::New();
	itkwriter->SetFileName( (filename_out).c_str() );
	itkwriter->SetInput(importer->GetOutput());
	itkwriter->Update();

	



}
