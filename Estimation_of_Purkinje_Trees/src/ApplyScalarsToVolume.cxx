/*=========================================================================
 
 Copyright Dec. 2011 Ruben Cardenes
 
 All rights reserved. 
 
 =========================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <vector>
#include <assert.h>

#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkPolyDataWriter.h"
#include "vtkPolyDataReader.h"
#include "vtkPointData.h"
#include "vtkFloatArray.h"
#include "vtkPointLocator.h"
#include "vtkIntArray.h"
#include "vtkIdList.h"
#include "vtkCell.h"

#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

double ComputeDistance(double *p, double *q) {
	return sqrt( (p[0]-q[0])*(p[0]-q[0]) + (p[1]-q[1])*(p[1]-q[1]) + (p[2]-q[2])*(p[2]-q[2]) );
}

int repeated(int id, std::vector<int> list_id) {
	int result = 0;
	for (unsigned int i=0;i<list_id.size();i++) {
		if (list_id[i] == id) {
			result = 1;
			break;
		}
	}
	return result;
}

int main( int argc, char * argv[] )
{
    char *input_file,*domain_file,*output_file,*array_name;
	float sum = 0;
    int debug = 0;
    if (argc < 4) {
		printf("Author:      Rubén Cárdenes, Dec. 2011\n");
		printf("Description: Este programa lee una superficie: (surface.vtk) y un volumen \n");
		printf("             (domain.vtk) y guarda los valores escalares de la superficie \n");
		printf("             asociados al valor 1 de domain.vtk en el volumen output.vtk\n");
		printf("Usage: %s surface.vtk domain.[mhd/vtk] output.[mhd/vtk] -array_name [name]\n",argv[0]);
		printf("          -d debug\n");
		printf("          -sum [0]\n");
		printf("Notes: -sum 100 will add 100 to the scalar values\n");
        return 1;
    } 
	
	input_file = argv[1];
	argc--;
	argv++;
	domain_file = argv[1];
	argc--;
	argv++;
	output_file = argv[1];
	argc--;
	argv++;
	int ok = 0;
	
	// Parse input arguments
	while ( argc > 1 ) {
		ok = false;               
		if ( ( ok == false ) && ( strcmp( argv[1], "-d" ) == 0 ) )
		{
			argc--;
			argv++;
			debug = 1;
			ok = true;
		}     
	    if ( ( ok == false ) && ( strcmp( argv[1], "-sum" ) == 0 ) )
	    {
			argc--;
			argv++;
			sum = atof(argv[1]);
			argc--;
			argv++;
			ok = true;
	    }      
	    if ( ( ok == false ) && ( strcmp( argv[1], "-array_name" ) == 0 ) )
	    {
			argc--;
			argv++;
			array_name = argv[1];
			argc--;
			argv++;
			ok = true;
	    }      
		
		
		if ( ok == false )
		{
			std::cout << "Can not parse argument " << argv[1] << std::endl;
			return 1;
		}
	}
	
	// Read surface
	vtkPolyDataReader* reader = vtkPolyDataReader::New();
	reader->SetFileName(input_file);
	reader->Update();
	
	typedef itk::Image< short,  3 >  ImageType;	
	typedef itk::ImageFileReader< ImageType >  ReaderType;
	typedef itk::ImageRegionIterator<ImageType> IteratorType;
	typedef itk::ImageRegionIteratorWithIndex<ImageType> IteratorWithIndexType;
	
	typedef itk::Image< float,  3 >  FloatImageType;	
	typedef itk::ImageFileWriter< FloatImageType >  WriterType;
	typedef itk::ImageRegionIterator<FloatImageType> FloatIteratorType;
	
	ReaderType::Pointer vol_reader = ReaderType::New();
	vol_reader->SetFileName( domain_file );
	try
    {
		vol_reader->Update();
    }
	catch( itk::ExceptionObject & err )
    {
		std::cout <<  "Reading volume: Exception caught: " 
		+ std::string(err.GetDescription()) ;
		return EXIT_FAILURE;
	}
	
	ImageType::Pointer inputImage = ImageType::New();
	inputImage = vol_reader->GetOutput();
	
	int max1 = inputImage->GetLargestPossibleRegion().GetSize()[1];
	int max2 = inputImage->GetLargestPossibleRegion().GetSize()[0];
	int max3 = inputImage->GetLargestPossibleRegion().GetSize()[2];
	//unsigned char* domain = (unsigned char*)malloc(sizeof(unsigned char)*max3*max2*max1);
	//IteratorType It( inputImage, inputImage->GetLargestPossibleRegion() );
	
	printf("input file: %s, dims: %d %d %d\n",input_file,max1,max2,max3);
	//int i=0;
	//for( It.GoToBegin(); !It.IsAtEnd(); ++It  ) {
	//	domain[i] = (unsigned short)It.Get();
	//	i++;
	//}
	
	//float spacing[3];
	//spacing[0] = inputImage->GetSpacing()[1];
	//spacing[1] = inputImage->GetSpacing()[0];
	//spacing[2] = inputImage->GetSpacing()[2];
	//float origin[3];
	//origin[0] = inputImage->GetOrigin()[1];
	//origin[1] = inputImage->GetOrigin()[0];
	//origin[2] = inputImage->GetOrigin()[2];
	
	// compute:
	vtkPolyData* surface = reader->GetOutput();
	
	vtkFloatArray* float_array = vtkFloatArray::New();
	float_array->Allocate(0);
	float_array = (vtkFloatArray*)surface->GetPointData()->GetArray(array_name);
	
	if (float_array->GetNumberOfTuples() == 0 ) {
		printf("No scalars found in array: %s \n",array_name);
		return 1;
	}
	printf("Number of scalars: %d\n",float_array->GetNumberOfTuples());

	vtkPointLocator* locator = vtkPointLocator::New();
	locator->SetDataSet(surface);
	locator->BuildLocator();
	
	FloatImageType::Pointer outputImage = FloatImageType::New();
	outputImage->SetSpacing(inputImage->GetSpacing());
	outputImage->SetOrigin(inputImage->GetOrigin());
	outputImage->SetRegions(inputImage->GetLargestPossibleRegion());
	outputImage->Allocate();
	outputImage->FillBuffer(0);
	
	IteratorWithIndexType Itw( inputImage, inputImage->GetLargestPossibleRegion() );
	FloatIteratorType outIt( outputImage, outputImage->GetLargestPossibleRegion() );
	int count = 0;
	
	vtkPolyData* polydata_out = vtkPolyData::New();
	polydata_out->Allocate();
	
	printf("Traverse the volume \n");
	std::vector<int> list_id;
	for( Itw.GoToBegin(), outIt.GoToBegin(); !Itw.IsAtEnd(); ++Itw, ++outIt ) {
		if (Itw.Get() == 1) {
			count++;
			ImageType::IndexType point_index = Itw.GetIndex();
			ImageType::PointType point;
			inputImage->TransformIndexToPhysicalPoint(point_index,point);
			double p[3];
			
			p[0] = point[0];p[1] = point[1];p[2] = point[2];
			int id = locator->FindClosestPoint(p);         
			float val = float_array->GetValue(id);
			outIt.Set(val + sum);
			//if ( repeated(id,list_id) == 0 ) {
			//	list_id.push_back(id);
			//	surface->GetPoint(id,p);
			//}
		}
	}
	
	printf("Count %d list_id.size %d \n",count,list_id.size());
	// Write 
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(output_file);
	writer->SetInput(outputImage);
	writer->Write();
	
	
	printf("Ok\n");
	
}
