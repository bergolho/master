/* Copyright (c) Ruben Cardenes Almeida 9/02/2010 */
/* Interfaz para llamar a SurfaceVoronoi */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <iostream>

#include <itkImage.h>
#include <itkImageRegionIterator.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

int EdgeDetect(unsigned char *domain,int max1,int max2,int max3, int in, int out, int label){
	int x,y,z,i;
	i = 0;
	for(z=0; z<max3; z++) {
		for(y=0; y<max1; y++) {
			for(x=0; x<max2; x++) {
				
				if ((z==0)||(x==0)||(y==0)||(z==max3-1)||(y==max1-1)||(x==max2-1)) {
					// nothing to do
				} else {
					if ( (domain[i]==in )&&
						((domain[i+1]== out)||
						 (domain[i-1]== out)||
						 
						 (domain[i+max2]== out)||
						 (domain[i-max2]== out)||
						 
						 (domain[i+max2*max1]== out)||
						 (domain[i-max2*max1]== out))) {	   	    
							
							domain[i]=label;    
						}
				}
				i++;
			}
		}
	}
	printf("Edge Detect Ok!\n");
	return 0;
	
}

int EdgeDetect2(unsigned char *domain,int max1,int max2,int max3, int out, int label){
	int x,y,z,i;
	i = 0;
	for(z=0; z<max3; z++) {
		for(y=0; y<max1; y++) {
			for(x=0; x<max2; x++) {
				
				if ((z==0)||(x==0)||(y==0)||(z==max3-1)||(y==max1-1)||(x==max2-1)) {
					// nothing to do
				} else {
					if ( (domain[i]!=out )&&
						((domain[i+1]== out)||
						 (domain[i-1]== out)||
						 
						 (domain[i+max2]== out)||
						 (domain[i-max2]== out)||
						 
						 (domain[i+max2*max1]== out)||
						 (domain[i-max2*max1]== out))) {	   	    
							
							domain[i]=label;    
						}
				}
				i++;
			}
		}
	}
	printf("Edge Detect Ok!\n");
	return 0;
	
}


int main(int argc, char* argv[]) {
	
	int in_value=0,out_value=255;
	int boundary_value=1;	
	int only_out_value=-1;	
	std::string filename;
	char *inputfile,*outputfile;
	int debug = 0;
	
    if (argc < 3) {
		printf("Author: Ruben Cardenes, Feb 2010 \n");
		printf("Usage: %s inputfile.mhd output.mhd [options]\n",argv[0]);
		printf("          [-in_value (=0)] \n");    
		printf("          [-out_value (=255)]\n");
		printf("          [-boundary_value(=1)]\n");
		printf("          [-only_out_value ]\n");
		printf(" Notes: inputfile.mhd representa el volumen de entrada. Es tipo unsigned_char \n");
		printf("        Se buscan los bordes entre in_value y out_value \n");
		printf("        (in_value) es 0 por defecto.  (out_value) es 255 por defecto. \n");
		printf("        El valor de los bordes en la salida (boundary_value) será 1 por defecto.\n");
		return 1;
    }
	
    inputfile = argv[1];
    argc--;
    argv++;
    outputfile = argv[1];
    argc--;
    argv++;
    int ok = 0;
    // Parse input arguments
    while ( argc > 1 ) {
		ok = false;
		
		if ( ( ok == false ) && ( strcmp( argv[1], "-in_value" ) == 0 ) )
		{
			argc--;
			argv++;
			in_value = atoi( argv[1] );
			argc--;
			argv++;
			ok = true;
		}
		
		if ( ( ok == false ) && ( strcmp( argv[1], "-out_value" ) == 0 ) )
		{
			argc--;
			argv++;
			out_value = atoi( argv[1] );
			argc--;
			argv++;
			ok = true;
		}
		if ( ( ok == false ) && ( strcmp( argv[1], "-boundary_value" ) == 0 ) )
		{
			argc--;
			argv++;
			boundary_value = atoi( argv[1] );
			argc--;
			argv++;
			ok = true;
		}
		if ( ( ok == false ) && ( strcmp( argv[1], "-only_out_value" ) == 0 ) )
		{
			argc--;
			argv++;
			only_out_value = atoi( argv[1] );
			argc--;
			argv++;
			ok = true;
		}
		if ( ( ok == false ) && ( strcmp( argv[1], "-debug" ) == 0 ) )
		{
			argc--;
			argv++;
			debug = 1;
			ok = true;
		}
		
		if ( ok == false )
		{
			std::cout << "Can not parse argument " << argv[1] << std::endl;
			return 1;
		}
	}
	
	typedef itk::Image< short,  3 >  ImageType;	
	typedef itk::ImageFileReader< ImageType >  ReaderType;
	typedef itk::ImageFileWriter< ImageType >  WriterType;
	typedef itk::ImageRegionIterator<ImageType> IteratorType;
	
	
	ReaderType::Pointer vol_reader = ReaderType::New();
	vol_reader->SetFileName( inputfile );
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
	
	unsigned char* domain = (unsigned char*)malloc(sizeof(unsigned char)*max3*max2*max1);
	IteratorType It( inputImage, inputImage->GetLargestPossibleRegion() );
	
	printf("input file: %s, dims: %d %d %d\n",inputfile,max1,max2,max3);
	int i=0;
	for( It.GoToBegin(); !It.IsAtEnd(); ++It  ) {
		domain[i] = (unsigned char)It.Get();
		i++; 
	}
	
	if (only_out_value >= 0) {
		EdgeDetect2(domain,max1,max2,max3, only_out_value, boundary_value);
	} else {
		EdgeDetect(domain,max1,max2,max3, in_value , out_value, boundary_value);
	}
	i=0;
	for( It.GoToBegin(); !It.IsAtEnd(); ++It  ) {
		It.Set(domain[i]);
		i++;        
	}
	// Write 
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(outputfile);
	writer->SetInput(inputImage);
	writer->Write();
	
	
	
}

