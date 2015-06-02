/*=========================================================================
 
 Copyright Dec. 2011 Ignacio Larrabide/Ruben Cardenes
 
 All rights reserved. 
 
 =========================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <vector>

#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkPolyDataWriter.h"
#include "vtkPolyDataReader.h"
#include "vtkPointData.h"
#include "vtkFloatArray.h"
#include "vtkIdList.h"
#include "vtkType.h"

double Distance(double* p1,double* p2) {
	return sqrt( (p1[0]-p2[0])*(p1[0]-p2[0]) + (p1[1]-p2[1])*(p1[1]-p2[1]) + (p1[2]-p2[2])*(p1[2]-p2[2]) );
}

double dotprod(double* v1,double* v2) {
	return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}

void vectprod(double* u,double* v, double* out) {
	out[0] = u[1]*v[2] - v[1]*u[2];
	out[1] = -u[0]*v[2] + v[0]*u[2];
	out[2] = u[0]*v[1] - v[0]*u[1];
}

int CheckIDInList(int id,std::vector<int> lista_elim) {
	if (lista_elim.size() == 0) return 0; 
	for (unsigned int i=0;i<lista_elim.size();i++) {
		if (id == lista_elim[i]) return 1;
	}
	return 0;
}

void GroupPoints(vtkPolyData* input, vtkPoints* output, float distance, vtkFloatArray* scalars_in, vtkFloatArray* scalars_out) {
	int npts = input->GetPoints()->GetNumberOfPoints();
	double pt0[3],pt1[3];
	std::vector<int> lista_elim;
	
	for (int id0=0;id0<npts-1;id0++) {
		if (CheckIDInList(id0,lista_elim) == 1 ) continue; 
        input->GetPoint(id0,pt0);
		for (int id1=id0+1;id1<npts;id1++) {
			input->GetPoint(id1,pt1);
			float aux = Distance(pt0,pt1);
			if (aux <= distance) {
			    // Eliminamos ese punto (id1) de la lista 
				lista_elim.push_back(id1);
			}
        } 
    }
	
	printf("lista_elim.size %d npts %d\n",(int)lista_elim.size(),npts);
	for (int id0=0;id0<npts;id0++) {
	   	if (CheckIDInList(id0,lista_elim) == 1 ) continue;
		input->GetPoint(id0,pt1);
        output->InsertNextPoint(pt1);
		float value = scalars_in->GetValue(id0);
		scalars_out->InsertNextValue(value);
	}
	
}

int main( int argc, char * argv[] )
{
    char *input_file,*output_file;
    int debug = 0;
	float distance = 10;
    if (argc < 3) {
		printf("Uso: %s input_points.vtk output_points.vtk \n",argv[0]);
		printf("        -d debug\n");
		printf("        -dist distance\n");
        return 1;
    } 
	
	input_file = argv[1];
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
		if ( ( ok == false ) && ( strcmp( argv[1], "-dist" ) == 0 ) )
		{
			argc--;
			argv++;
			distance = atof(argv[1]);
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
	
	vtkPolyDataReader* reader = vtkPolyDataReader::New();
	reader->SetFileName(input_file);
	reader->Update();
	
	if (reader->GetOutput()->GetPoints()->GetNumberOfPoints() == 0) {
		printf("file %s, has 0 points",input_file);
		return -1;
	}
	
	vtkPoints* points_out = vtkPoints::New();
	points_out->Allocate(0);
	vtkFloatArray* scalars_in = vtkFloatArray::New();
	scalars_in->Allocate(0);
	scalars_in = (vtkFloatArray*)reader->GetOutput()->GetPointData()->GetScalars();
	vtkFloatArray* scalars_out = vtkFloatArray::New();
	scalars_out->Allocate(0);
	
	//////////////////////////////////////////////////////////////
	printf("Number of input points %d \n",reader->GetOutput()->GetPoints()->GetNumberOfPoints());
	
	GroupPoints(reader->GetOutput(), points_out, distance, scalars_in, scalars_out);
	
	printf("Number of output points %d \n",points_out->GetNumberOfPoints());
	
	
	vtkPolyData* polydata_out=vtkPolyData::New();
	polydata_out->Allocate();
	polydata_out->SetPoints(points_out);
	polydata_out->GetPointData()->SetScalars(scalars_out);
	
	vtkPolyDataWriter* writer = vtkPolyDataWriter::New();
	writer->SetFileName(output_file);
	writer->SetInput(polydata_out);
	writer->Update();
	printf("Ok\n");
	
}
