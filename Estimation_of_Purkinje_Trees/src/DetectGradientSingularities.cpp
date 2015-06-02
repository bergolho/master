/*
 *  DetectGradientSingularities.cpp
 *  
 *
 *  Created by Ruben Cardenes Almeida on 26/02/13.
 *  Copyright 2013 Universitat Pompeu Fabra. All rights reserved.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <vector>
#include <sys/time.h>
#include <itkImage.h>
#include <itkImageRegionIterator.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <iostream>
#include <fstream>

#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkFloatArray.h>
#include <vtkIdList.h>

#include "DToptimo3d.h"

int maptox(int mapindex,int max1, int max2) {
	return (mapindex % (max1*max2) ) % max2;
}

int maptoy(int mapindex,int max1 ,int max2) {
	int a = mapindex % (max1*max2);
	return ( a - a % max2 ) /max2;
}

int maptoz(int mapindex,int max1,int max2) {
	return (mapindex - (mapindex % (max1*max2))) / (max1*max2);
}

float dot_prod(float *a, float *b ) {
	float c = a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
	return c;
}


//int mapIndex3D(int r,int c,int z, int nr,int nc,int nz)
//{
//	if (c >= nc) return -1;
//	if (c < 0) return -1;
//	if (r >= nr) return -1;
//	if (r < 0) return -1;
//	if (z >= nz) return -1;
//	if (z < 0) return -1;
//	return c + r * nc + z * nr * nc;
//}

float compute_det(float mat[3][3]) {
	float value = mat[0][0]*mat[1][1]*mat[2][2] + mat[0][1]*mat[1][2]*mat[2][0] + mat[0][2]*mat[1][0]*mat[2][1] -
	mat[2][0]*mat[1][1]*mat[0][2] - mat[1][0]*mat[0][1]*mat[2][2] - mat[2][1]*mat[1][2]*mat[0][0];
	return value;
}

float comp_norm(float *a) {
	float b = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
	return b;
}

void normalize(float *g) {
	float r = sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]);
	g[0] = g[0]/r;g[1] = g[1]/r;g[2] = g[2]/r;
}

void compute_normal(unsigned char* domain, int *v_index, float* n, float* maps_SD, int max1,int max2,int max3) {
	int index = mapIndex3D(v_index[0],v_index[1],v_index[2],max1,max2,max3);
	float w,w1 = 0,w2 = 0,w0 = 0;
	n[0] = 0;n[1] = 0;n[2] = 0;
	
	//      n[1] += maps_SD[index] - maps_SD[index+1];
	//	n[1] += maps_SD[index-1] - maps_SD[index];
	//	
	//	n[0] += maps_SD[index] - maps_SD[index+max2];
	//	n[0] += maps_SD[index-max2] - maps_SD[index];
	//	
	//	n[2] += maps_SD[index] - maps_SD[index+max1*max2];
	//	n[2] += maps_SD[index-max1*max2] - maps_SD[index];
	
	for (int x=-1;x<=1;x++) {
		for (int y=-1;y<=1;y++) {
			for (int z=-1;z<=1;z++) {
				if (x == 0 && y==0 && z == 0) continue;
				int newmapindex = index + z*max1*max2 + y*max2 + x;
				if (maps_SD[newmapindex] != 0) {
					w = 1/sqrt(abs(x)+abs(y)+abs(z));
					if (abs(x)>0) {
						if (x > 0) n[1] += w*(maps_SD[index] - maps_SD[newmapindex]);
						if (x < 0) n[1] += w*(maps_SD[newmapindex] - maps_SD[index]);
						w1 += w;
					}
					if (abs(y)>0) {
						if (y > 0) n[0] += w*(maps_SD[index] - maps_SD[newmapindex]);
						if (y < 0) n[0] += w*(maps_SD[newmapindex] - maps_SD[index]);
						w0 += w;
					}
					if (abs(z)>0) {
						if (z > 0) n[2] += w*(maps_SD[index] - maps_SD[newmapindex]);
						if (z < 0) n[2] += w*(maps_SD[newmapindex] - maps_SD[index]);
						w2 += w;
					}
				}
			}
		}
	}
	n[0] = n[0]/w0;
	n[1] = n[1]/w1;
	n[2] = n[2]/w2;
	
	// Si queremos mas precision, usar los vecinos en esquinas
}

void compute_gradient2(float* fm_output, int *v_index, float* g, int max1,int max2,int max3, int invert_gradient) {
	int index = mapIndex3D(v_index[0],v_index[1],v_index[2],max1,max2,max3);
	float w,w1 = 0,w2 = 0,w0 = 0;
	g[0] = 0;g[1] = 0;g[2] = 0;
	for (int x=-3;x<=3;x++) {
		for (int y=-3;y<=3;y++) {
			for (int z=-3;z<=3;z++) {
				if (x == 0 && y==0 && z == 0) continue;
				int newmapindex = index + z*max1*max2 + y*max2 + x;
				if (newmapindex < 0 || newmapindex > max1*max2*max3) continue;
				if (fm_output[newmapindex] > 0) {
					w = 1/sqrt(abs(x)+abs(y)+abs(z));
					if (abs(x)>0) {
						if (x > 0) g[1] += w*(fm_output[index] - fm_output[newmapindex]);
						if (x < 0) g[1] += w*(fm_output[newmapindex] - fm_output[index]);
						w1 += w;
					}
					if (abs(y)>0) {
						if (y > 0) g[0] += w*(fm_output[index] - fm_output[newmapindex]);
						if (y < 0) g[0] += w*(fm_output[newmapindex] - fm_output[index]);
						w0 += w;
					}
					if (abs(z)>0) {
						if (z > 0) g[2] += w*(fm_output[index] - fm_output[newmapindex]);
						if (z < 0) g[2] += w*(fm_output[newmapindex] - fm_output[index]);
						w2 += w;
					}
				}
			}
		}
	}
	
	if (w0 != 0)  g[0] = g[0]/w0;
	if (w1 != 0)  g[1] = g[1]/w1;
	if (w2 != 0)  g[2] = g[2]/w2;
	
	if (invert_gradient) {
		g[0] = -g[0];
		g[1] = -g[1];
		g[2] = -g[2];
	}
	
}


void compute_second_derivative(unsigned char* domain, float** grad, int *v_index, float sg[3][3], int max1,int max2,int max3, int invert_gradient) {
	int index = mapIndex3D(v_index[0],v_index[1],v_index[2],max1,max2,max3);
	float w,w1 = 0,w2 = 0,w0 = 0;
	sg[0][0] = 0;sg[1][0] = 0;sg[2][0] = 0;
	sg[0][1] = 0;sg[1][1] = 0;sg[2][1] = 0;
	sg[0][2] = 0;sg[1][2] = 0;sg[2][2] = 0;

	for (int x=-3;x<=3;x++) {
		for (int y=-3;y<=3;y++) {
			for (int z=-3;z<=3;z++) {
				if (x == 0 && y==0 && z == 0) continue;
				int newmapindex = index + z*max1*max2 + y*max2 + x;
				if (newmapindex < 0 || newmapindex > max1*max2*max3) continue;
				if (domain[newmapindex] == 1) {
					w = 1/sqrt(abs(x)+abs(y)+abs(z));
					if (abs(x)>0) {
						if (x > 0) sg[1][1] += w*(grad[1][index] - grad[1][newmapindex]);
						if (x < 0) sg[1][1] += w*(grad[1][newmapindex] - grad[1][index]);
						
						if (x > 0) sg[1][0] += w*(grad[0][index] - grad[0][newmapindex]);
						if (x < 0) sg[1][0] += w*(grad[0][newmapindex] - grad[0][index]);
						
						if (x > 0) sg[1][2] += w*(grad[2][index] - grad[2][newmapindex]);
						if (x < 0) sg[1][2] += w*(grad[2][newmapindex] - grad[2][index]);
						w1 += w;
					}
					if (abs(y)>0) {
						if (y > 0) sg[0][0] += w*(grad[0][index] - grad[0][newmapindex]);
						if (y < 0) sg[0][0] += w*(grad[0][newmapindex] - grad[0][index]);
						
						if (y > 0) sg[0][1] += w*(grad[1][index] - grad[1][newmapindex]);
						if (y < 0) sg[0][1] += w*(grad[1][newmapindex] - grad[1][index]);
						
						if (y > 0) sg[0][2] += w*(grad[2][index] - grad[2][newmapindex]);
						if (y < 0) sg[0][2] += w*(grad[2][newmapindex] - grad[2][index]);
						w0 += w;
					}
					if (abs(z)>0) {
						if (z > 0) sg[2][2] += w*(grad[2][index] - grad[2][newmapindex]);
						if (z < 0) sg[2][2] += w*(grad[2][newmapindex] - grad[2][index]);
						
						if (z > 0) sg[2][1] += w*(grad[1][index] - grad[1][newmapindex]);
						if (z < 0) sg[2][1] += w*(grad[1][newmapindex] - grad[1][index]);
						
						if (z > 0) sg[2][0] += w*(grad[0][index] - grad[0][newmapindex]);
						if (z < 0) sg[2][0] += w*(grad[0][newmapindex] - grad[0][index]);
						w2 += w;
					}
					
					
				}
			}
		}
	}
	
	if (w0 != 0) {
		sg[0][0] = sg[0][0]/w0;
	    sg[0][1] = sg[0][1]/w0;
	    sg[0][2] = sg[0][2]/w0;
	}
	
	if (w1 != 0) {
		sg[1][0] = sg[1][0]/w1;
	    sg[1][1] = sg[1][1]/w1;
	    sg[1][2] = sg[1][2]/w1;
	}
	
	if (w2 != 0) {
		sg[2][0] = sg[2][0]/w2;
		sg[2][1] = sg[2][1]/w2;
		sg[2][2] = sg[2][2]/w2;
	}
}



void ComputeHessian(unsigned char* domain, float* fm, float* out_vol, vtkPoints* out_points, 
								 int max1, int max2, int max3, float* spacing, float* origin, 
								 vtkFloatArray* Scalars, int invert_gradient, float threshold,
								 int mode_sink, int mode_source,float* maps_SD_float) {
	float g[3],normal[3];
	float sg[3][3];
	float v_index[3];
	int v_index_int[3];
	int count = 0;
	
	float **grad;
	grad = (float**)malloc(sizeof(float*)*3);
	for (int i=0;i<3;i++) {
		grad[i] = (float*)malloc(max1*max2*max3*sizeof(float));
	}
	
	for (int i = 0;i<max1*max2*max3;i++) {	
		if (domain[i] == 1 ) { 
			v_index_int[0] = maptoy(i,max1,max2);
			v_index_int[1] = maptox(i,max1,max2);
			v_index_int[2] = maptoz(i,max1,max2);
			v_index[0] = (float)v_index_int[0];
			v_index[1] = (float)v_index_int[1];
			v_index[2] = (float)v_index_int[2];
			
			compute_gradient2(fm, v_index_int, g, max1, max2, max3, invert_gradient);
			// project normals to surface
			compute_normal(domain, v_index_int, normal, maps_SD_float, max1, max2, max3);
			normalize(normal);
			
			float scalar_prod = g[0]*normal[0] + g[1]*normal[1] + g[2]*normal[2];
			g[0] = g[0] - scalar_prod*normal[0];
			g[1] = g[1] - scalar_prod*normal[1];
			g[2] = g[2] - scalar_prod*normal[2];
			
			//normalize(g);
			grad[0][i] = g[0];
			grad[1][i] = g[1];
			grad[2][i] = g[2];
		}
	}
	
	for (int i = 0;i<max1*max2*max3;i++) {	
		if (domain[i] == 1 ) { 
			v_index_int[0] = maptoy(i,max1,max2);
			v_index_int[1] = maptox(i,max1,max2);
			v_index_int[2] = maptoz(i,max1,max2);
			v_index[0] = (float)v_index_int[0];
			v_index[1] = (float)v_index_int[1];
			v_index[2] = (float)v_index_int[2];
						
			compute_second_derivative(domain, grad, v_index_int, sg, max1, max2, max3, invert_gradient);
			float det = compute_det(sg);
			out_vol[i] = det;
			
			//if (v_index_int[0] == 58 && v_index_int[1] == 148 && v_index_int[2] == 66) {
			//if (v_index_int[0] == 52 && v_index_int[1] == 141 && v_index_int[2] == 61) {
			//	printf("sg: [%f %f %f], [%f %f %f], [%f %f %f] det: %f\n",sg[0][0],sg[0][1],sg[0][2],
			//		   sg[1][0],sg[1][1],sg[1][2],sg[2][0],sg[2][1],sg[2][2],det);
			//}
			if (mode_sink) {
				if (det > threshold) {
					out_points->InsertNextPoint(v_index[1]*spacing[0]+origin[0],v_index[0]*spacing[1]+origin[1],v_index[2]*spacing[2]+origin[2]);
					Scalars->InsertNextValue(det+5);
				}
			}
			if (mode_source) {
				if (det < threshold) {
					out_points->InsertNextPoint(v_index[1]*spacing[0]+origin[0],v_index[0]*spacing[1]+origin[1],v_index[2]*spacing[2]+origin[2]);
					Scalars->InsertNextValue(det+5);
				}
			}
			
			count++;
		}
	}
	
	printf("Num sing points: %d\n",count);
	
}

float compute_std(std::vector<float> angle_vector) {
    float m = 0,v = 0;
	if (angle_vector.size() == 0) return 0;
	for (unsigned int i=0;i<angle_vector.size();i++) {
		m += angle_vector[i];
	}
	m = m/(float)angle_vector.size();
	for (unsigned int i=0;i<angle_vector.size();i++) {
	    v += (angle_vector[i]-m)*(angle_vector[i]-m);
	}
    v = sqrt(v/(float)angle_vector.size());
	return v;
}


void DetectGradientSingularities(unsigned char* domain, float* fm, vtkPoints* out_points, 
								 int max1, int max2, int max3, float* spacing, float* origin, 
								 vtkFloatArray* Scalars, int invert_gradient, float threshold,
								 int mode_sink, int mode_source, float* maps_SD_float) {
	float g[3],gn[3],d[3],normal[3];
	float v_index[3];
	int v_index_int[3],v_index_new[3];
	int count = 0;
	std::vector<float> angle_vector;
	for (int i = 0;i<max1*max2*max3;i++) {	
		if (domain[i] == 1 ) { 
			v_index_int[0] = maptoy(i,max1,max2);
			v_index_int[1] = maptox(i,max1,max2);
			v_index_int[2] = maptoz(i,max1,max2);
			v_index[0] = (float)v_index_int[0];
			v_index[1] = (float)v_index_int[1];
			v_index[2] = (float)v_index_int[2];
			
			compute_gradient2(fm, v_index_int, g, max1, max2, max3, invert_gradient);
			normalize(g);
			
			// project normals to surface
			compute_normal(domain, v_index_int, normal, maps_SD_float, max1, max2, max3);
			normalize(normal);
			float scalar_prod = g[0]*normal[0] + g[1]*normal[1] + g[2]*normal[2];
			g[0] = g[0] - scalar_prod*normal[0];
			g[1] = g[1] - scalar_prod*normal[1];
			g[2] = g[2] - scalar_prod*normal[2];
			/////////////
			angle_vector.clear();
			for (int x=-1;x<=1;x++) {
				for (int y=-1;y<=1;y++) {
					for (int z=-1;z<=1;z++) {
						if (x == 0 && y==0 && z == 0) continue;
						int newmapindex = i + z*max1*max2 + y*max2 + x;
						if (newmapindex < 0 || newmapindex > max1*max2*max3) continue;
						if (domain[newmapindex] == 1) {
							v_index_new[0] = maptoy(newmapindex,max1,max2);
							v_index_new[1] = maptox(newmapindex,max1,max2);
							v_index_new[2] = maptoz(newmapindex,max1,max2);
							compute_gradient2(fm, v_index_new, gn, max1, max2, max3, invert_gradient);
							normalize(gn);
							
							// project normals to surface
							compute_normal(domain, v_index_int, normal, maps_SD_float, max1, max2, max3);
							normalize(normal);
							scalar_prod = gn[0]*normal[0] + gn[1]*normal[1] + gn[2]*normal[2];
							gn[0] = gn[0] - scalar_prod*normal[0];
							gn[1] = gn[1] - scalar_prod*normal[1];
							gn[2] = gn[2] - scalar_prod*normal[2];
							/////////////
							
							double product = (double)dot_prod(g,gn);
							double angle = acos(product);
							//printf("--> angle %f dot_product %f g %f %f %f, gn: %f %f %f\n",angle, product, g[0],g[1],g[2],gn[0],gn[1],gn[2]);
							if (angle > 3.1415927*threshold) {
								d[0] = (float)v_index_new[0] - (float)v_index_int[0];
								d[1] = (float)v_index_new[1] - (float)v_index_int[1];
								d[2] = (float)v_index_new[2] - (float)v_index_int[2];
								normalize(d);
								if ( dot_prod(d,gn) > 0 && mode_sink == 1) {
									// Sink points 
									angle_vector.push_back(angle);
								}
								if ( dot_prod(d,gn) < 0 && mode_source == 1) {
									// Source points
									angle_vector.push_back(angle);
								}
								//printf("angle %f\n",angle);
							}
						}
					}
				}
			}
			
			if (mode_sink == 1 && angle_vector.size() > 2) {
				float std_dev = compute_std(angle_vector);
				if (std_dev > 0.1) { 
			      //Scalars->InsertNextValue(0.75);
				  Scalars->InsertNextValue(fm[i]);
			      out_points->InsertNextPoint(v_index[1]*spacing[0]+origin[0],v_index[0]*spacing[1]+origin[1],v_index[2]*spacing[2]+origin[2]);
			      count++;
				}
			}
			if (mode_source == 1 && angle_vector.size() > 2) {
				float std_dev = compute_std(angle_vector);
				if (std_dev > 0.1) { 
			      //Scalars->InsertNextValue(1);
				  Scalars->InsertNextValue(fm[i]);
			      out_points->InsertNextPoint(v_index[1]*spacing[0]+origin[0],v_index[0]*spacing[1]+origin[1],v_index[2]*spacing[2]+origin[2]);
			      count++;
				}
			}
			
						
		}
	}
	printf("Num sing points: %d\n",count);
	
}

int main(int argc, char* argv[]) {
	
	char *domainfile,*fmfile,*outputfile;
	int debug = 0;
	int invert_gradient = 0;
	float threshold = 0.6;
	int mode_sink = 0, mode_source = 0; 
	if (argc < 4) {
		printf("Author: Ruben Cárdenes, Feb. 2010\n");
		printf("Description: program that computes the singularities in a surface gradient field\n");
		printf("Usage: %s domain_file.[vtk/mhd] fm.[vtk/mhd] output_points.vtk [options]\n",argv[0]);
		printf("                         -debug\n"); 
		printf("                         -threshold\n"); 
		printf("                         -invert_gradient\n");
		printf("                         -mode_sink\n");
		printf("                         -mode_source\n");
		printf("Notes: domain_file.mhd volume representing the surface, 1 is foreground,.\n");
		return 1;
	}
	
	
    // domainfile = edge_open
	domainfile = argv[1];
	argc--;
	argv++;
    
    // fmfile = LAT_volume
	fmfile = argv[1];
	argc--;
	argv++;
    
    // outputfile = singular_pointsLAT
	outputfile = argv[1];
	argc--;
	argv++;
	int ok = 0;
    
	// Parse input arguments
	while ( argc > 1 ) {
		ok = false;
		
		if ( ( ok == false ) && ( strcmp( argv[1], "-debug" ) == 0 ) )
		{
			argc--;
			argv++;
			debug = 1;
			ok = true;
		}
		if ( ( ok == false ) && ( strcmp( argv[1], "-threshold" ) == 0 ) )
		{
			argc--;
			argv++;
			threshold = atof(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		
		if ( ( ok == false ) && ( strcmp( argv[1], "-invert_gradient" ) == 0 ) )
		{
			argc--;
			argv++;
			invert_gradient = 1;
			ok = true;
		}
		if ( ( ok == false ) && ( strcmp( argv[1], "-mode_sink" ) == 0 ) )
		{
			argc--;
			argv++;
			mode_sink = 1;
			ok = true;
		}
		if ( ( ok == false ) && ( strcmp( argv[1], "-mode_source" ) == 0 ) )
		{
			argc--;
			argv++;
			mode_source = 1;
			ok = true;
		}
		
		if ( ok == false )
		{
			cerr << "Can not parse argument " << argv[1] << endl;
			return 1;
		}
    }
	
	if (mode_sink == 0 && mode_source == 0) mode_sink = 1;
	///////////////////////////////////////////////////////////////
	// Read 3D image
	
	typedef itk::Image< short,  3 >  ImageType;	
	typedef itk::ImageFileReader< ImageType >  ReaderType;
	typedef itk::ImageRegionIterator<ImageType> IteratorType;
	
	typedef itk::Image< float,  3 >  FloatImageType;	
	typedef itk::ImageFileReader< FloatImageType >  FloatReaderType;
	typedef itk::ImageFileWriter< FloatImageType >  FloatWriterType;
	typedef itk::ImageRegionIterator<FloatImageType> FloatIteratorType;
	
    // Load the readers for the images
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( domainfile );
	try
    {
		reader->Update();
    }
	catch( itk::ExceptionObject & err )
    {
		std::cout <<  "Reading volume: Exception caught: " 
		+ std::string(err.GetDescription()) ;
		return EXIT_FAILURE;
    }
	
    // Read input image containig the surface+LAT
    ImageType::Pointer inputImage = ImageType::New();
	inputImage = reader->GetOutput();

    // Obtaining dimension of each axis to build a solid cube with the largest size edge
	int max1 = inputImage->GetLargestPossibleRegion().GetSize()[1];
	int max2 = inputImage->GetLargestPossibleRegion().GetSize()[0];
	int max3 = inputImage->GetLargestPossibleRegion().GetSize()[2];
	unsigned char* domain = (unsigned char*)malloc(sizeof(unsigned char)*max3*max2*max1);
    
	IteratorType It( inputImage, inputImage->GetLargestPossibleRegion() );
	int i=0;
	for( It.GoToBegin(); !It.IsAtEnd(); ++It  ) {
		domain[i] = (unsigned short)It.Get();		
		i++;
	}
	float spacing[3];
	spacing[0] = inputImage->GetSpacing()[0];
	spacing[1] = inputImage->GetSpacing()[1];
	spacing[2] = inputImage->GetSpacing()[2];
	float origin[3];
	origin[0] = inputImage->GetOrigin()[0];
	origin[1] = inputImage->GetOrigin()[1];
	origin[2] = inputImage->GetOrigin()[2];
	
	/////// Read FM
	FloatReaderType::Pointer float_reader = FloatReaderType::New();
	float_reader->SetFileName( fmfile );
	try
    {
		float_reader->Update();
    }
	catch( itk::ExceptionObject & err )
    {
		std::cout <<  "Reading volume: Exception caught: " 
		+ std::string(err.GetDescription()) ;
		return EXIT_FAILURE;
    }
	FloatImageType::Pointer float_inputImage = FloatImageType::New();
	float_inputImage = float_reader->GetOutput();
	
	float* fm = (float*)malloc(sizeof(float)*max3*max2*max1);
	FloatIteratorType FIt( float_inputImage, float_inputImage->GetLargestPossibleRegion() );
	
	i=0;
	for( FIt.GoToBegin(); !FIt.IsAtEnd(); ++FIt  ) {
		fm[i] = (float)FIt.Get();
		i++;
	}
    /////////// 
		
	vtkPolyData* out_data = vtkPolyData::New();
	out_data->Allocate();
	
	vtkPoints* out_points = vtkPoints::New();
	out_points->Allocate(0);
	
	vtkFloatArray* Scalars = vtkFloatArray::New();
	Scalars->Allocate(0);
	
	// Compute maps_SD_float
	char * maps_aux = (char*)malloc(sizeof(char)*max1*max2*max3);
	float * maps_SD_float = (float*)malloc(sizeof(float)*max1*max2*max3);
	unsigned char* domain_SD = (unsigned char*)malloc(sizeof(unsigned char)*max1*max2*max3);
	for (int i = 0;i<max1*max2*max3;i++) {
		if (domain[i] == 1) {
			domain_SD[i] = domain[i];
		} else {
			domain_SD[i] = 0;
		}
	}
    // returns maps_SD_float
	DToptimo3d(domain_SD, max1, max2, max3, 2, maps_aux, maps_SD_float, 5);
	free(maps_aux);
	free(domain_SD);
	
	//Le doy la vuelta al signo de los de dentro:
	for (int i = 0;i<max1*max2*max3;i++) {
		if (domain[i] == 2 || domain[i] == 0) {
			maps_SD_float[i] = -maps_SD_float[i];
		}
	}
	/////////
	
    // returns out_points
	DetectGradientSingularities(domain, fm, out_points, max1, max2, max3, spacing, 
								origin, Scalars, invert_gradient, threshold, 
								mode_sink, mode_source,maps_SD_float);
		
	
	//float* out_vol = (float*)calloc(max3*max2*max1,sizeof(float));
//	for (int i=0;i<max1*max2*max3;i++) out_vol[i]=-1;
//	ComputeHessian(domain, fm, out_vol, out_points, max1, max2, max3, spacing, 
//								origin, Scalars, invert_gradient, threshold, 
//								mode_sink, mode_source,maps_SD_float);
//	
//	FloatImageType::Pointer outputImage = FloatImageType::New();
//	outputImage->SetSpacing(inputImage->GetSpacing());
//	outputImage->SetOrigin(inputImage->GetOrigin());
//	outputImage->SetRegions(inputImage->GetLargestPossibleRegion());
//	outputImage->Allocate();
//	outputImage->FillBuffer(0);
//	FloatIteratorType outIt( outputImage, outputImage->GetLargestPossibleRegion() );
//	
//	i=0;
//	for( outIt.GoToBegin(); !outIt.IsAtEnd(); ++outIt ) {
//		if (out_vol[i] == -1) {
//			outIt.Set(0);
//		} else {
//			outIt.Set(out_vol[i]+5);
//		}
//		i++;
//	}
//	
//	FloatWriterType::Pointer imwriter = FloatWriterType::New();
//	imwriter->SetFileName(outputfile);
//	imwriter->SetInput(outputImage);
//	imwriter->Write();
	
	out_data->SetPoints(out_points);
	out_data->GetPointData()->SetScalars(Scalars);
	
	vtkPolyDataWriter* writer = vtkPolyDataWriter::New();
	writer->SetInput(out_data);
	writer->SetFileName(outputfile);
	writer->Update();
	
	printf("Ok\n");

}