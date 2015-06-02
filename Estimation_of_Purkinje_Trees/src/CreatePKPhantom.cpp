/*
 *  Backtrace.cpp
 *  
 *
 *  Created by Ruben Cardenes Almeida on 20/11/12.
 *  Copyright 2012 Universitat Pompeu Fabra. All rights reserved.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <vector>
#include <list>
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
#include <vtkFloatArray.h>
#include <vtkPolyDataWriter.h>

struct punto_struct_type {
	float x;
	float y;
	float z;
};

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


int mapIndex3D(int r,int c,int z, int nr,int nc,int nz)
{
	if (c >= nc) return -1;
	if (c < 0) return -1;
	if (r >= nr) return -1;
	if (r < 0) return -1;
	if (z >= nz) return -1;
	if (z < 0) return -1;
	return c + r * nc + z * nr * nc;
}

float dot_prod(float *a, float *b ) {
	float c = a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
	return c;
}

double distance(double* p, double* q) {
	return sqrt( (q[0]-p[0])*(q[0]-p[0]) + (q[1]-p[1])*(q[1]-p[1]) + (q[2]-p[2])*(q[2]-p[2])  );
}



float comp_norm(float *a) {
	float b = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
	return b;
}

void normalize(float *g) {
	float r = sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]);
	g[0] = g[0]/r;g[1] = g[1]/r;g[2] = g[2]/r;
}

int check_neighborhood(int index_a,int index_b,int max1,int max2) {
	// Devuelve 0 si no son vecinos, 1 si lo son:
	int x,y,z,newmapindex;
	for (x=-1;x<=1;x++) {
		for (y=-1;y<=1;y++) {
			for (z=-1;z<=1;z++) {
				if (x == 0 && y == 0 && z == 0) continue;
				newmapindex = index_a + z*max1*max2 + y*max2 + x;
				if (newmapindex == index_b) {
					return 1;
				}
			}
		}
	}
	return 0;
}

int check_num_connected_neighbors(std::vector<int> neigh,int max1,int max2) {
	std::vector<int> aux,prop;
	int index;
	
	if (neigh.size() == 0) return 0;
	for (unsigned int i=0;i<neigh.size();i++) {
		aux.push_back(-1);
	}
	
	int num_connected_neighbors = 0;
	for (unsigned int j=0;j<neigh.size();j++) {
		if (aux[j] == -1) {
			prop.push_back(neigh[j]);
			num_connected_neighbors++;
			aux[j] = num_connected_neighbors;
			while (prop.size() > 0) {
				int end = prop.size()-1;
				index = prop[end];
				prop.pop_back();
				for (unsigned int i=0;i<neigh.size();i++) {
					if (index == neigh[i] || aux[i] > 0) continue;
					if (check_neighborhood(index,neigh[i],max1,max2) == 1) {
						prop.push_back(neigh[i]);
						aux[i] = num_connected_neighbors;
					}
				}
			}
		}
	}
	
	return num_connected_neighbors;
	
}  

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

// funcion que a partir de un pto devuelve:
// 0 (ningun vecino):          pto aislado
// 1 (un vecino):              pto extremo
// 2 (dos vecinos conexos)   : pto superfluo 
// 3 (dos vecinos disconexos): pto intermedio
// 4 (tres vecinos disconexos): pto triple 
// 5 (tres vecinos disconexos): pto intermedio 
// 6 (tres vecinos disconexos): pto superfluo 
// 7 (dos vecinos disconexos, que son bifurcacion, con otro pto en commun): pto superfluo 
// 8 (mas de tres vecinos): pto superfluo 
int CheckPointState(unsigned short *vol_esqueleto,int mapindex,int max1, int max2, int max3) {
	int state = -1;
	int i,j,k,x,y,z,newindex,vecino_indep = 0,flag_aux;
	int n=0;
	int count = 0;
	std::vector<int> neigh,neighN0,neighN1;
	std::vector<int> index_aux;
	
	k = maptoz(mapindex,max1,max2);
	i = maptoy(mapindex,max1,max2);
	j = maptox(mapindex,max1,max2);
	neigh.clear();
	for (x=-1;x<=1;x++) {
		for (y=-1;y<=1;y++) {
			for (z=-1;z<=1;z++) {
				if (x == 0 && y == 0 && z == 0) continue;
				if (k+z < 0 || i+y < 0 || j+x < 0 || k+z > max3-1 || i+y > max1-1 || j+x > max2-1) continue; // bound checking
				newindex = mapIndex3D(i+y,j+x,k+z,max1,max2,max3);
				count++;
				if (vol_esqueleto[newindex] > 0 ) {
					n++;
					neigh.push_back(newindex);
				}
			}// end z  
		} // end y  
	} // end x
	
	
	if (n==0) { // pto aislado
        state = 0;
		return state;
	}
	if (n==1) { // pto extremo
        state = 1;
		return state;
	}
	if (n==2) { // 2 vecinos
		if (check_neighborhood(neigh[0],neigh[1],max1,max2) != 0) { 
			state = 2; 
			return state; // pto superfluo 
		} else { // Both neihbors are separated 
			if (vol_esqueleto[neigh[0]] == 5 && vol_esqueleto[neigh[1]] == 5) { // both are separated bifurcations 
				// miramos los vecinos de neigh[0] y los vecinos de neihg[1]
				neighN1.clear();
				k = maptoz(neigh[1],max1,max2);
				i = maptoy(neigh[1],max1,max2);
				j = maptox(neigh[1],max1,max2);
				for (x=-1;x<=1;x++) {
					for (y=-1;y<=1;y++) {
						for (z=-1;z<=1;z++) {
							if (x == 0 && y == 0 && z == 0) continue;
							if (k+z < 0 || i+y < 0 || j+x < 0 || k+z > max3-1 || i+y > max1-1 || j+x > max2-1) continue; // bound checking
							newindex = mapIndex3D(i+y,j+x,k+z,max1,max2,max3);
							if (vol_esqueleto[newindex] > 0 && newindex != mapindex) {
								neighN1.push_back(newindex);
							}
						}// end z  
					} // end y  
                } // end x
				neighN0.clear();
				k = maptoz(neigh[0],max1,max2);
				i = maptoy(neigh[0],max1,max2);
				j = maptox(neigh[0],max1,max2);
				for (x=-1;x<=1;x++) {
					for (y=-1;y<=1;y++) {
						for (z=-1;z<=1;z++) {
							if (x == 0 && y == 0 && z == 0) continue;
							if (k+z < 0 || i+y < 0 || j+x < 0 || k+z > max3-1 || i+y > max1-1 || j+x > max2-1) continue; // bound checking
							newindex = mapIndex3D(i+y,j+x,k+z,max1,max2,max3);
							if (vol_esqueleto[newindex] > 0 && newindex != mapindex) {
								neighN0.push_back(newindex);
							}
						}// end z  
					} // end y  
                } // end x
                
                int neighbor_in_common = 0;
                for (unsigned int u=0;u<neighN0.size();u++) {
                    for (unsigned int v=0;v<neighN1.size();v++) {
						if (neighN0[u] == neighN1[v]) neighbor_in_common = 1;
                    } 
                }
                if (neighbor_in_common == 1) state = 7;
			} else {
				state = 3;
				return state; // pto intermedio
			}
        }
	}
	
	if (n>=5) {
        state = 8;
     	return state;
	}
	
	if (n>=3) { // puede ser pto doble, triple o eliminable 
		// Sera pto triple si tres de sus vecinos no se tocan entre si (son independientes)
		// O si uno es independiente y los otros dos estan separados por una distancia > 1
		int num_independent_neighbors =0;
		for (i=0;i<n;i++) {
			flag_aux = 0;
	        for (j=0;j<n;j++) {
				if (j==i) continue;
				if (check_neighborhood(neigh[i],neigh[j],max1,max2) != 0) {
					flag_aux = 1;
					break;
				}
			}
			if (flag_aux == 0) num_independent_neighbors++;
		}
		int num_connected_neighbors = check_num_connected_neighbors(neigh, max1, max2);
		//std::cout << "num_independent_neighbors " << num_independent_neighbors << std::endl; 
		if (n==3 && num_independent_neighbors == 1 ) { // pto triple 
			// Busco el que es independiente:
			for (i=0;i<n;i++) {
				flag_aux = 0;
				for (j=0;j<n;j++) {
					if (j==i) continue;
					if (check_neighborhood(neigh[i],neigh[j],max1,max2) != 0) {
						flag_aux = 1;
						break;
					}
				}
				if (flag_aux == 0) vecino_indep = i;
			}
			// Comparamos los no independientes:
			i =(vecino_indep+1)%3;
			j =(vecino_indep+2)%3;	
			if (distance2(neigh[i],neigh[j],max1,max2) > 1 ) {
				state = 4;  // pto triple
				return state;
			} else {
				state = 5;  // pto intermedio
				return state;
			}
		} 
		
		if (num_independent_neighbors >= 3 ) { // pto triple 
			state = 4;
			return state;
		} 
		if (num_independent_neighbors == 0 &&  num_connected_neighbors == 1) { // pto eliminable: los tres vecinos estan juntos
			state = 6;
			return state;
		} 	
    } // end if n>=3
	
    return state;
	
} 


// Label continuously a tree phantom starting from a seed using a time step  
void ComputePKPhantom2(unsigned char* domain,float *fm,int start_index[3],int max1,int max2,int max3, float step,float spacing) {
	
	std::list<int> front;
	for (int i=0;i<max1*max2*max3;i++) {
		if (fm[i] > 0 && domain[i] == 1 && maptox(i,max1,max2) == start_index[1] 
			&& maptoy(i,max1,max2) == start_index[0] && maptoz(i,max1,max2) == start_index[2]) {
			fm[i] = 3;
			front.push_back(i);
			printf("found start index \n");
			break;
		} 
	}
	float dtot = 0;
	int count = 0, count_ant = 0;
	while (front.size() != 0) {
		int index= front.front();
		front.pop_front();
		
		count_ant = count; 
		for (int x=-1;x<=1;x++) {
			for (int y=-1;y<=1;y++) {
				for (int z=-1;z<=1;z++) {
					if (x == 0 && y==0 && z == 0) continue;
					int newmapindex = index + z*max1*max2 + y*max2 + x;
					if (newmapindex < 0 || newmapindex > max1*max2*max3) continue;
					if (fm[newmapindex] == 1) {
						front.push_back(newmapindex);
						fm[newmapindex] = fm[index] + step; 
						dtot += sqrt(abs(x)+abs(z)+abs(y))*spacing; 
						count++;
					}
				}
			}
		}
		
		if (count_ant == count) {
			for (int x=-2;x<=2;x++) {
				for (int y=-2;y<=2;y++) {
					for (int z=-2;z<=2;z++) {
						if (x == 0 && y==0 && z == 0) continue;
						int newmapindex = index + z*max1*max2 + y*max2 + x;
						if (newmapindex < 0 || newmapindex > max1*max2*max3) continue;
						if (fm[newmapindex] == 1) {
							front.push_back(newmapindex);
							fm[newmapindex] = fm[index] + step; 
							dtot += sqrt(abs(x)+abs(z)+abs(y))*spacing; 
							count++;
						}
					}
				}
			}
		}
		
	}
	printf("Estimated vel: %f (dtot %f step %f spacing %f\n",dtot/(step*(float)count),dtot,step,spacing);
	printf("count %d\n",count);
	
}


void ComputePKPhantom(unsigned char* domain,float *fm,int start_index[3],int max1,int max2,int max3) {
	
	for (int i=0;i<max1*max2*max3;i++) {
		if (fm[i] > 0 && domain[i] != 1) {
			// Buscamos el pto mas en domain
			int d[27],candidate[27];
			int k = 0;
			for (int x=-1;x<=1;x++) {
				for (int y=-1;y<=1;y++) {
					for (int z=-1;z<=1;z++) {
						if (x == 0 && y==0 && z == 0) continue;
						int newmapindex = i + z*max1*max2 + y*max2 + x;
						if (newmapindex < 0 || newmapindex > max1*max2*max3) continue;
						if (domain[newmapindex] == 1 && fm[newmapindex] == 0) {
							d[k] = abs(x) + abs(y) + abs(z);
							candidate[k] = newmapindex;
							k++;
						}
					}
				}
			}
			int min_dist = 100;
			int final_index = 1;
			for (int j=0;j<k;j++) {
				if (d[j] < min_dist) {
					final_index = candidate[j];
					min_dist = d[j];
				}
			}
			fm[final_index] = 1;
			fm[i]=0;
			if (maptox(i,max1,max2) == 15 && maptoy(i,max1,max2) == 108 && maptoz(i,max1,max2) == 120 ) {
			    printf("final_index %d %d %d \n",maptox(final_index,max1,max2),maptoy(final_index,max1,max2),maptoz(final_index,max1,max2));
			}
		} // end if
	}
	
}

void IncludeMidpoints(unsigned short* outcenterline_short, float* fm,
					  std::vector<int> vec_index_esq, std::vector<int> &vec_index_esq_selected,
					  int max1,int max2,int max3, float dist_min_entre_puntos) {
	
	for (unsigned int k=0;k<vec_index_esq.size();k++) {
		int index = vec_index_esq[k];
		int state = CheckPointState(outcenterline_short, index, max1, max2, max3);
		if (state == 3 || state == 5) {
			float min_dist = 10000;
			for (unsigned int j=0;j<vec_index_esq_selected.size();j++) {
				int index_selected = vec_index_esq_selected[j];
				float dist = distance2(index,index_selected,max1,max2);
				if (dist < min_dist) {
					min_dist = dist;
				}
			}
			if (min_dist > dist_min_entre_puntos && fm[index] > 22) {
				vec_index_esq_selected.push_back(index);
			}
		}
	}
	
}

void ReduceToBifAndExtrema(unsigned char* domain,float* fm,int max1,int max2,int max3, float dist_min_entre_puntos,
						   vtkPoints* points,vtkFloatArray* scalars, float* spacing,float* origin) {
	std::vector<int> vec_index_esq;
	unsigned short* outcenterline_short = (unsigned short*)malloc(sizeof(unsigned short)*max1*max2*max3);
	for (int i=0;i<max1*max2*max3;i++) {
		if (fm[i] > 0) {
			outcenterline_short[i] = 1;
			vec_index_esq.push_back(i);
		} else {
			outcenterline_short[i] = 0;
		}
	}
	
	int flag = 1,iteracion =0,removed = 0;
	while (flag) { // While any superfluous exists
		printf("vec_index_esq.size() %d it %d removed %d\n",(int)vec_index_esq.size(),iteracion,removed);
		flag = 0;
		for (unsigned k=0;k<vec_index_esq.size();k++) {
			int index = vec_index_esq[k];
			//printf("index, %d outcenterline_short[index] %d\n",index,outcenterline_short[index]);getchar();
			if (outcenterline_short[index] > 0) { 
				int state = CheckPointState(outcenterline_short, index, max1, max2, max3);
				//printf("index, %d state %d\n",index,state);getchar();
				if (state == 0) {outcenterline_short[index] = 0;fm[index] = 0;flag=1;removed++;}
				if (state == 1) outcenterline_short[index] = 4;
				if (state == 2 || state == 6) {outcenterline_short[index] = 0;fm[index] = 0;flag=1;removed++;}
				if (state == 3 || state == 5) outcenterline_short[index] = 2;
				if (state == 4) outcenterline_short[index] = 5;
				if (state == 7) {outcenterline_short[index] = 0;fm[index] = 0;flag=1;removed++;}
				if (state == 8) {outcenterline_short[index] = 0;fm[index] = 0;flag=1;removed++;}
			}
		}
		iteracion++;
	}	
	
	// Nos quedamos solo con los extremos y las bifurcaciones
	std::vector<int> vec_index_esq_selected;
	for (unsigned int k=0;k<vec_index_esq.size();k++) {
		int index = vec_index_esq[k];
		if (outcenterline_short[index] > 0) { 
			int state = CheckPointState(outcenterline_short, index, max1, max2, max3);
			if ((state == 1 || state == 4) && fm[index] > 30) {
				vec_index_esq_selected.push_back(index);
			}
		}
	}
	
	printf("vec_index_esq_selected.size %d\n",(int)vec_index_esq_selected.size());
	IncludeMidpoints(outcenterline_short, fm, vec_index_esq, vec_index_esq_selected,
					 max1, max2, max3, dist_min_entre_puntos);
	
	printf("vec_index_esq_selected.size %d\n",(int)vec_index_esq_selected.size());
	
	// For some reason some points of the tree are outsite the edge. We put them inside with this piece of code	
	int count1 = 0,count2 = 0;
	int replace_index = 0;
	for (unsigned int k=0;k<vec_index_esq_selected.size();k++) {
		int index = vec_index_esq_selected[k];
		if (fm[index] > 0 && domain[index] != 1) {
			count1++;
			int continue_flag = 1;
			for (int x=-1;x<1;x++) {
				for (int y=-1;y<1;y++) {
					for (int z=-1;z<1;z++) {
						if (x == 0 && y == 0 && z == 0) continue;
						int newmapindex = index + z*max1*max2 + y*max2 + x;
						if (newmapindex < 0 || newmapindex > max1*max2*max3) continue;
						if (domain[newmapindex] == 1 && continue_flag == 1) {
							count2++;
							fm[newmapindex] = fm[index];
							fm[index] = 0;
							continue_flag = 0;
							replace_index = newmapindex;
						}
					}
				}
			}
			if (continue_flag == 0)	vec_index_esq_selected[k] = replace_index;
		}
	}
	printf("count1 %d, count2 %d\n",count1,count2);
	count1 = 0,count2 = 0;
	for (unsigned int k=0;k<vec_index_esq_selected.size();k++) {
		int index = vec_index_esq_selected[k];
		if (fm[index] > 0 && domain[index] != 1) {
			count1++;
			int continue_flag = 1;
			for (int x=-2;x<2;x++) {
				for (int y=-2;y<2;y++) {
					for (int z=-2;z<2;z++) {
						if (x == 0 && y == 0 && z == 0) continue;
						int newmapindex = index + z*max1*max2 + y*max2 + x;
						if (newmapindex < 0 || newmapindex > max1*max2*max3) continue;
						if (domain[newmapindex] == 1 && continue_flag == 1) {
							count2++;
							fm[newmapindex] = fm[index];
							fm[index] = 0;
							continue_flag = 0;
							replace_index = newmapindex;
						}
					}
				}
			}
			if (continue_flag == 0)	vec_index_esq_selected[k] = replace_index;
		}
	}
	printf("count1 %d, count2 %d\n",count1,count2);
	
	printf("vec_index_esq_selected.size %d\n",(int)vec_index_esq_selected.size());
	
	// Put the selected values into fm 
	std::vector<float> vec_value_esq_selected;
	for (unsigned int k=0;k<vec_index_esq_selected.size();k++) {
		int index = vec_index_esq_selected[k];
		vec_value_esq_selected.push_back(fm[index]);
	}
	
	float p[3];
	for (int i=0;i<max1*max2*max3;i++) {fm[i] = 0;}
	for (unsigned int k=0;k<vec_index_esq_selected.size();k++) {
		int index = vec_index_esq_selected[k];
		fm[index] = vec_value_esq_selected[k];
		
		int y = maptoy(index,max1,max2);
		int x = maptox(index,max1,max2);
		int z = maptoz(index,max1,max2);
		p[0] = x*spacing[0]+origin[0];
		p[1] = y*spacing[1]+origin[1];
		p[2] = z*spacing[2]+origin[2];
		points->InsertNextPoint(p);
		scalars->InsertNextValue(fm[index]);
	}
	
	printf("vec_index_esq_selected.size %d\n",(int)vec_index_esq_selected.size());
	
}

int main(int argc, char* argv[]) {
	
	char *domainfile,*PKfile,*outputfile,*vtkpoints_outputfile;
	int debug = 0,vtkpoints_storing = 0;
	float dist_min_entre_puntos = 10, baseline_val = 0;
	//int start_index[3] = {165,105,103};	
	int start_index[3] = {28,103,167};	
    float step = 0.1;
	if (argc < 4) {
		printf("Author: Ruben Cárdenes, Feb. 2010\n");
		printf("Description: program that creates simple tree phantom for Purkinje simple\n");
		printf("             experiment the input is the domain file (edge.vtk) defining \n");
		printf("             the edges of the LV  and the volume defining a tree (PKtree_vol).\n");
		printf("             The output will be same PK tree volume with scalar values on  \n");
		printf("             it defining the arrival times. In this version the generated tree\n");
		printf("             is discontinuous, it generates points every point_interdist \n");
		printf("Usage: %s edge.[mhd/vtk] PKtree_vol.[mhd/vtk] output_vol.vtk [options]\n",argv[0]);
		printf("                         -debug\n"); 
		printf("                         -start_index [x y z]\n"); 
		printf("                         -step [0.1]\n"); 
		printf("                         -point_interdist [10]\n"); 
		printf("                         -vtkpoints_outputfile [filename]\n"); 
		printf("Notes: edge.[mhd/vtk] volume representing the LV surface, with labels:  \n");
		printf("            1 for the surface, 0 for the interior, 255 for the exterior.\n");
		printf("       point_interdist: minima distance entre end terminals (default = 10)  \n");
		printf("       step: time step used to label the tree from the HIS \n");
		printf("       Optionally the output can be stored as vtkpoints in a file using\n"); 
		printf("            the -vtkpoints_outputfile [filename] option\n"); 
		return 1;
	}
	
	
	domainfile = argv[1];
	argc--;
	argv++;
	PKfile = argv[1];
	argc--;
	argv++;
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
		
		
		if ( ( ok == false ) && ( strcmp( argv[1], "-start_index" ) == 0 ) )
		{
			argc--;
			argv++;
			start_index[0] = atoi( argv[1] );
			argc--;
			argv++;
			start_index[1] = atoi( argv[1] );
			argc--;
			argv++;
			start_index[2] = atoi( argv[1] );
			argc--;
			argv++;
			ok = true;
		}
		
		if ( ( ok == false ) && ( strcmp( argv[1], "-step" ) == 0 ) )
		{
			argc--;
			argv++;
			step = atof(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ( ( ok == false ) && ( strcmp( argv[1], "-point_interdist" ) == 0 ) )
		{
			argc--;
			argv++;
			dist_min_entre_puntos = atof(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ( ( ok == false ) && ( strcmp( argv[1], "-vtkpoints_outputfile" ) == 0 ) )
		{
			argc--;
			argv++;
			vtkpoints_outputfile = argv[1];
			vtkpoints_storing = 1;
			argc--;
			argv++;
			ok = true;
		}
		if ( ( ok == false ) && ( strcmp( argv[1], "-baseline_val" ) == 0 ) )
		{
			argc--;
			argv++;
			baseline_val = atof(argv[1]);
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
	
	///////////////////////////////////////////////////////////////
	// Read 3D image
	typedef itk::Image< short,  3 >  ImageType;	
	typedef itk::ImageFileReader< ImageType >  ReaderType;
	typedef itk::ImageRegionIterator<ImageType> IteratorType;
	
	typedef itk::Image< float,  3 >  FloatImageType;	
	typedef itk::ImageFileReader< FloatImageType >  FloatReaderType;
	typedef itk::ImageFileWriter< FloatImageType >  FloatWriterType;
	typedef itk::ImageRegionIterator<FloatImageType> FloatIteratorType;
	
	// Read domain
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
	ImageType::Pointer inputImage = ImageType::New();
	inputImage = reader->GetOutput();
	
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
	
	//printf("domain file: %s, dims: %d %d %d\n",domainfile,max1,max2,max3);
	std::vector<int> seeds_list;
	int mapindex = mapIndex3D(start_index[0],start_index[1],start_index[2],max1,max2,max3);
	seeds_list.push_back(mapindex);
	printf("index %d %d %d, domain[mapindex] %d\n",start_index[0],start_index[1],start_index[2],domain[mapindex]);
	
	/////// Read PK 
	FloatReaderType::Pointer float_reader = FloatReaderType::New();
	float_reader->SetFileName( PKfile );
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
	
	float* fm = (float*)malloc(sizeof(float)*max3*max2*max1);
	FloatIteratorType FIt( float_reader->GetOutput(), float_reader->GetOutput()->GetLargestPossibleRegion() );
	
	i=0;
	int count_aux = 0;
	for( FIt.GoToBegin(); !FIt.IsAtEnd(); ++FIt  ) {
		fm[i] = (float)FIt.Get();
		if (fm[i] > 0) {count_aux++;}
		i++;
	}
	
	printf("Pkfile %s PK[mapindex] %f, count_aux %d, i %d\n",PKfile,fm[mapindex],count_aux,i);

	/////// Create the phantom tree with increasing arrival time values 
	ComputePKPhantom2(domain,fm,start_index,max1,max2,max3,step,spacing[0]);
	
	/////// Reduce the phantom tree to the extrema and bifurcations
	vtkPoints* points = vtkPoints::New();
	points->Allocate(0);
	vtkFloatArray* scalars = vtkFloatArray::New();
	scalars->Allocate(0);
	
	float p[3];
	for (int k=0;k<max1*max2*max3;k++) {
		if (fm[k] > 0) {
		  fm[k] = fm[k]+baseline_val;
  		  int y = maptoy(k,max1,max2);
		  int x = maptox(k,max1,max2);
		  int z = maptoz(k,max1,max2);
		  p[0] = x*spacing[0]+origin[0];
		  p[1] = y*spacing[1]+origin[1];
		  p[2] = z*spacing[2]+origin[2];
		  points->InsertNextPoint(p);
		  scalars->InsertNextValue(fm[k]);
		}
	}

	//ReduceToBifAndExtrema(domain,fm,max1,max2,max3,dist_min_entre_puntos,points,scalars,spacing,origin);
	
	FloatImageType::Pointer outputImage = FloatImageType::New();
	outputImage->SetSpacing(float_reader->GetOutput()->GetSpacing());
	outputImage->SetOrigin(float_reader->GetOutput()->GetOrigin());
	outputImage->SetRegions(float_reader->GetOutput()->GetLargestPossibleRegion());
	outputImage->Allocate();
	outputImage->FillBuffer(0);
	
	FloatIteratorType outIt( outputImage, outputImage->GetLargestPossibleRegion() );
	
	i = 0;
	int num_sing_points = 0;
	float max_value = -1;
	for( outIt.GoToBegin(); !outIt.IsAtEnd(); ++outIt  ) {
		outIt.Set(fm[i]);
		if (fm[i] > 0) num_sing_points++;
		if (max_value < fm[i]) {
			max_value = fm[i];
		}
		i++;
	}
	printf("max_value %f\n",max_value);
	printf("num sing points %d\n",num_sing_points);
	
	// Write result
	printf("output vol file %s\n",outputfile);
	FloatWriterType::Pointer writer = FloatWriterType::New();
	writer->SetFileName(outputfile);
	writer->SetInput(outputImage);
	writer->Write(); 
	
	if (vtkpoints_storing == 1) {
		printf("output vtk_points file %s\n",vtkpoints_outputfile);
		vtkPolyData* polydata = vtkPolyData::New();
		polydata->Allocate();
		polydata->SetPoints(points);
		polydata->GetPointData()->SetScalars(scalars);
		
		vtkPolyDataWriter* writer = vtkPolyDataWriter::New();
		writer->SetInput(polydata);
		writer->SetFileName(vtkpoints_outputfile);
		writer->Update(); 
	}
	printf("Ok\n");
}
