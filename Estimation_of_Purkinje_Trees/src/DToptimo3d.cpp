/* Copyright (c) Ruben Cardenes Almeida 22/03/2002 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <assert.h>

#define MAX_ELEM_IN_BUCKET 1000000
#define NUM_BUCKETS 400
#define MAXPATTERNS (16384*4)
#define MAXCLASSNUMBER MAXPATTERNS
#define MAXDIM 20
#ifndef UCHAR
#define UCHAR(c) ((unsigned char)(c))
#endif

#define round(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))

struct element {
  int x;
  int y;
  int z;
  char icur;
  short dcur;
  char idcur;
  float df;
};

struct bucket {
  int num_elem;
  int *index_elem;
  int *index_l;
};

struct nodeDataNew {
  int id;
  int pclass;
  float d[MAXDIM];
  int row;
  int col;
  int slice;
};

struct nodeDataNew prototypeNodeData[MAXPATTERNS];

int numelembucket[NUM_BUCKETS];
int numrechazos = 0;
int numasignaciones = 0;
int asignacionesraras = 0;
int numPrototypes;
int highestIndexClass;
int actualDimension; 
int numPrototypesInClass[MAXCLASSNUMBER];
char buffer[2048];
int pdim;

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

int mapIndex2D(int r,int c, int nr,int nc)
{
  if (c >= nc) return -1;
  if (c < 0) return -1;
  if (r >= nr) return -1;
  if (r < 0) return -1;
  return c + r * nc;
}

//float distance(int x1,int y1,int x2,int y2) {
//  return sqrt ((x1-x2) * (x1-x2) + (y1-y2) * (y1-y2));
//}

float distance3d(int x1,int y1,int z1, int x2,int y2,int z2) {
  return sqrt ((x1-x2) * (x1-x2) + (y1-y2) * (y1-y2) + (z1-z2) * (z1-z2));
}

void assign(struct element *Element_p, int l,int mapindex,int d,char *maps,float* maps_float, struct element *proto) {
  //maps[(*Element_p).icur][mapindex] = proto[l].icur;
  maps[mapindex] = proto[l].icur;
  /* si queremos mapas de distancias */
  if ( maps_float[mapindex]  != 0 ) maps_float[mapindex] = (*Element_p).df;
  //maps_float[mapindex] = (*Element_p).df;
	
  (*Element_p).icur++;
  if ((*Element_p).dcur != d) {
    (*Element_p).idcur = (*Element_p).icur;
  }
  (*Element_p).dcur = d;
  numasignaciones++;
}

int propagate3d(struct element Element_p, int l,int max1,int max2, int max3, struct element *proto,struct bucket *Bucket, char* maps, int dcur,struct element *Element) {
  int new_mapindex,dist,x,y,z,siguiente;
  int icuraux;
	float distreal;

  for (x=-1;x<2;x=x+2) {
    y = 0; z = 0;
     if (Element_p.y + y >=0 && Element_p.x + x >=0 && Element_p.z + z >=0
	&& Element_p.y + y < max2 && Element_p.x + x < max1 && Element_p.z + z < max3) {
      new_mapindex = mapIndex3D(Element_p.x + x,Element_p.y + y,Element_p.z + z,max1,max2,max3);
      if (abs(Element_p.x - proto[l].x) < abs(Element_p.x + x - proto[l].x) || 
	  abs(Element_p.y - proto[l].y) < abs(Element_p.y + y - proto[l].y) || 
	  abs(Element_p.z - proto[l].z) < abs(Element_p.z + z - proto[l].z) ) {
	icuraux = Element[new_mapindex].icur>0?Element[new_mapindex].icur-1:0;
	if (Element[new_mapindex].icur < 1 && maps[new_mapindex] != l) { 
		distreal = distance3d(Element_p.x+x,Element_p.y+y,Element_p.z+z,proto[l].x,proto[l].y,proto[l].z);
		dist = round(distreal);	
		//dist = round(distance3d(Element_p.x+x,Element_p.y+y,Element_p.z+z,proto[l].x,proto[l].y,proto[l].z));
	  if (dist > dcur+1 || dist < dcur) continue;  
	  /* printf("%d Element_p.x+x %d, Element_p.y+y %d \n",dist,Element_p.x+x,Element_p.y+y);*/
	  /* put p,l in Bucket d=dist */
	  siguiente = Bucket[dist].num_elem;
	  if (siguiente >= MAX_ELEM_IN_BUCKET) {
	    printf("Excedidos num elem in bucket %d\n",MAX_ELEM_IN_BUCKET);
	    return 1;
	  }
	  /* control 
	  if (dist == 2) {
	    printf("y %d, x %d, d %d icuraux %d, icur %d l %d maps[new_mapindex] %d\n",Element_p.y+y,Element_p.x+x,dcur,icuraux,Element[new_mapindex].icur,l,maps[new_mapindex]);
	  }
	  */
	  Bucket[dist].index_elem[siguiente] = new_mapindex ;
	  Bucket[dist].index_l[siguiente] = l;
	  Bucket[dist].num_elem++;	
	  numelembucket[dist]++;
	  if ( Element[new_mapindex].df > distreal) {
		Element[new_mapindex].df = distreal;
	  }
	  }
      }
    } 
  }

  for (y=-1;y<2;y=y+2) {
    x = 0; z = 0;
    if (Element_p.y + y >=0 && Element_p.x + x >=0 && Element_p.z + z >=0
	&& Element_p.y + y < max2 && Element_p.x + x < max1 && Element_p.z + z < max3) {
      new_mapindex = mapIndex3D(Element_p.x + x,Element_p.y + y,Element_p.z + z,max1,max2,max3);
      if (abs(Element_p.x - proto[l].x) < abs(Element_p.x + x - proto[l].x) || 
	  abs(Element_p.y - proto[l].y) < abs(Element_p.y + y - proto[l].y) ||
	  abs(Element_p.z - proto[l].z) < abs(Element_p.z + z - proto[l].z) ) {
        icuraux = Element[new_mapindex].icur>0?Element[new_mapindex].icur-1:0;
	if (Element[new_mapindex].icur < 1 && maps[new_mapindex] != l) { 
		distreal = distance3d(Element_p.x+x,Element_p.y+y,Element_p.z+z,proto[l].x,proto[l].y,proto[l].z);
		dist = round(distreal);	
		//dist = round(distance3d(Element_p.x+x,Element_p.y+y,Element_p.z+z,proto[l].x,proto[l].y,proto[l].z));
	  if (dist > dcur+1 || dist < dcur) continue;  
	  /*printf("%d Element_p.x+x %d, Element_p.y+y %d \n",dist,Element_p.x+x,Element_p.y+y);*/
	  /* put p,l in Bucket d=dist */
	  siguiente = Bucket[dist].num_elem;
	  if (siguiente >= MAX_ELEM_IN_BUCKET) {
	    printf("Excedidos num elem in bucket: %d\n",MAX_ELEM_IN_BUCKET);
	    return 1;
	  }
	  /* control 
	  if (dist == 2) {
	     printf("y %d, x %d, d %d icuraux %d, icur %d l %d maps[new_mapindex] %d\n",Element_p.y+y,Element_p.x+x,dcur,icuraux,Element[new_mapindex].icur,l,maps[new_mapindex]);
	  }
	  */ 
	  Bucket[dist].index_elem[siguiente] = new_mapindex ;
	  Bucket[dist].index_l[siguiente] = l;
	  Bucket[dist].num_elem++;
	  numelembucket[dist]++;
		if ( Element[new_mapindex].df > distreal) {
			Element[new_mapindex].df = distreal;
		}
	  
	}
      }
    } 
  }

  for (z=-1;z<2;z=z+2) {
    x = 0; y=0;
    if (Element_p.y + y >=0 && Element_p.x + x >=0 && Element_p.z + z >=0
	&& Element_p.y + y < max2 && Element_p.x + x < max1 && Element_p.z + z < max3) {
      new_mapindex = mapIndex3D(Element_p.x + x,Element_p.y + y,Element_p.z + z,max1,max2,max3);
      if (abs(Element_p.x - proto[l].x) < abs(Element_p.x + x - proto[l].x) || 
	  abs(Element_p.y - proto[l].y) < abs(Element_p.y + y - proto[l].y) ||
	  abs(Element_p.z - proto[l].z) < abs(Element_p.z + z - proto[l].z) ) {
        icuraux = Element[new_mapindex].icur>0?Element[new_mapindex].icur-1:0;
	if (Element[new_mapindex].icur < 1 && maps[new_mapindex] != l) { 
		distreal = distance3d(Element_p.x+x,Element_p.y+y,Element_p.z+z,proto[l].x,proto[l].y,proto[l].z);
		dist = round(distreal);	
	  //  dist = round(distance3d(Element_p.x+x,Element_p.y+y,Element_p.z+z,proto[l].x,proto[l].y,proto[l].z));
	  if (dist > dcur+1 || dist < dcur) continue;  
	  /*printf("%d Element_p.x+x %d, Element_p.y+y %d \n",dist,Element_p.x+x,Element_p.y+y);*/
	  /* put p,l in Bucket d=dist */
	  siguiente = Bucket[dist].num_elem;
	  if (siguiente >= MAX_ELEM_IN_BUCKET) {
	    printf("Excedidos num elem in bucket: %d\n",MAX_ELEM_IN_BUCKET);
	    return 1;
	  }
	  /* control 
	  if (dist == 2) {
	     printf("y %d, x %d, d %d icuraux %d, icur %d l %d maps[new_mapindex] %d\n",Element_p.y+y,Element_p.x+x,dcur,icuraux,Element[new_mapindex].icur,l,maps[new_mapindex]);
	  }
	  */ 
	  Bucket[dist].index_elem[siguiente] = new_mapindex ;
	  Bucket[dist].index_l[siguiente] = l;
	  Bucket[dist].num_elem++;
	  numelembucket[dist]++;
		if ( Element[new_mapindex].df > distreal) {
			Element[new_mapindex].df = distreal;
		}
	  
	}
      }
    } 
  }
  
  return 0;
}

int DToptimo3d(unsigned char* prototypes,int max1, int max2, int max3, int K,char* maps,float* maps_float, int d_max ) {
  int i,j,x,y,z,l,r,count;
  int siguiente,indice_actual;
  int d,mapindex,buckets_empty;
  struct element *Element;
  struct element *Proto;
  struct bucket Bucket[NUM_BUCKETS];
  float distancia_from_l,distancia_ultima;

  /* inicializamos los mapas de etiquetas a -1 */
	for (j=0;j<max1*max2*max3;j++) {
      maps[j] = -1;
	  maps_float[j] = -1;
	  
    }

  /* reservamos memoria para los prototipos, para los Elementos 
     y para el bucket inicial (bucket 0) */
  Proto=(struct element*)malloc(sizeof(struct element)*max1*max2*max3);
  Element=(struct element*)malloc(sizeof(struct element)*max1*max2*max3);
  Bucket[0].index_elem = (int*)malloc(sizeof(int)*MAX_ELEM_IN_BUCKET);
  Bucket[0].index_l = (int*)malloc(sizeof(int)*MAX_ELEM_IN_BUCKET);  
  memset(Bucket[0].index_elem,0,sizeof(int)*MAX_ELEM_IN_BUCKET);
  memset(Bucket[0].index_l,0,sizeof(int)*MAX_ELEM_IN_BUCKET);
  Bucket[0].num_elem = 0;

  if (Element == (struct element*)NULL || Proto == (struct element*)NULL) {
    printf("Out of memory \nClose aplications \n and try again\n");
    return 1;
  }

  /* inicializamos los elementos, los prototipos, y el bucket inicial */
  i=0;
  for (z=0;z<max3;z++) {
    for (x=0;x<max1;x++) {
      for (y=0;y<max2;y++) {
		  if (prototypes[j] == 0 || prototypes[j] == 255) {
			  Element[i].df = 999999;
		  } else {
			  Element[i].df = 0;
		  }
		  if (prototypes[i] == 255) { 
	 	    Element[i].icur = 1;
		  } else {
	        Element[i].icur = 0;
 	        Element[i].dcur = -1;
	        Element[i].idcur = 0;
	        Element[i].x = x;
	        Element[i].y = y;
	        Element[i].z = z;
		  }
		  i++;
      }
    }
  }
  l = 0;
  count = 0;
  for (z=0;z<max3;z++) {
    for (x=0;x<max1;x++) {
       for (y=0;y<max2;y++) {
		 if (prototypes[count] > 0 && prototypes[count] < 255) {
			 Proto[l].x= x;
			 Proto[l].y= y;
			 Proto[l].z= z;
			 Proto[l].icur = prototypes[count];
			 //Proto[l].icur = 1;
			 Bucket[0].index_elem[l] = count;
			 Bucket[0].index_l[l] = l;
			 //printf("prototypes[count] %d\n",prototypes[count]);
			 //Bucket[0].index_l[l] = prototypes[count];
			 Bucket[0].num_elem++;
			 maps_float[count] = 0;
			 l++;
		 }
		 count++;
      }
    }
  }

  printf("hay %d prototipos \n",l);
  for (i=1;i<NUM_BUCKETS;i++) { 
    Bucket[i].num_elem = 0;
    numelembucket[i] = 0;
  }
  numelembucket[0] = Bucket[0].num_elem;
  
  /* Fin de inicializacion */
  d=0;
  while (1) {
    /*printf("Distancia actual: %d\n",d);*/
    Bucket[d+1].index_elem = (int*)malloc(sizeof(int)*MAX_ELEM_IN_BUCKET);
    Bucket[d+1].index_l = (int*)malloc(sizeof(int)*MAX_ELEM_IN_BUCKET);
    /* printf("reservada memoria para bucket d+1\n" ); */
    while (Bucket[d].num_elem != 0) {      
      /* Get (p,l) from Bucket d */
      /* printf("Obteniendo de bucket: %d\n",d); */
      siguiente = Bucket[d].num_elem-1;
      /* printf("num_elem = %d \n",Bucket[d].num_elem); */
      indice_actual = Bucket[d].index_elem[siguiente];
      l=Bucket[d].index_l[siguiente];      
      Bucket[d].index_elem[siguiente] = -1;      
      Bucket[d].num_elem--;     
      mapindex =mapIndex3D(Element[indice_actual].x,Element[indice_actual].y,Element[indice_actual].z,max1,max2,max3);
      if (mapindex != indice_actual) {
	printf("Error catastrofico: mapindex %d != indice_actual %d \n",mapindex,indice_actual);
	return 1;
      }
      /* printf("mapindex %d\n",mapindex);  */
      /*         *****         */  

      if (Element[indice_actual].icur < 1) {
	if (Element[indice_actual].dcur < d) {	  
	  assign(&Element[indice_actual],l,mapindex,d,maps,maps_float,Proto);
	  /* printf("asignado pixel\n");*/
	  if (propagate3d(Element[indice_actual],l,max1,max2,max3,Proto,Bucket,maps,d,Element) != 0) return 1;
	} else {
	  for (j=Element[indice_actual].idcur;j<=Element[indice_actual].icur;j++) {
	    if (maps[mapindex] == (char)l) { numrechazos++; break;}		
	    if (j == Element[indice_actual].icur) {  
	      asignacionesraras++;	
	      distancia_from_l = distance3d(Element[indice_actual].x,Element[indice_actual].y,Element[indice_actual].z,Proto[l].x,Proto[l].y,Proto[l].z);
	      r = maps[mapindex];
	      distancia_ultima = distance3d(Element[indice_actual].x,Element[indice_actual].y,Element[indice_actual].z,Proto[r].x,Proto[r].y,Proto[r].z);
	      if (distancia_from_l < distancia_ultima ) {		
		//maps[mapindex] = (char)l;
			  maps[mapindex] = Proto[l].icur;
			  
		assign(&Element[indice_actual],r,mapindex,d,maps,maps_float,Proto);
	      } else {	          
		assign(&Element[indice_actual],l,mapindex,d,maps,maps_float,Proto);
	      }
	      if (propagate3d(Element[indice_actual],l,max1,max2,max3,Proto,Bucket,maps,d,Element) != 0) return 1;
	      break;
	    }
	  }/* end for */
	} /*end else*/
      } /* end if (Element_p.icur < 1) */
    } /*end while */

    if (Bucket[d].num_elem == 0) {
      /* printf("num_elem en bucket %d = %d\n",d,numelembucket[d]);*/
      free(Bucket[d].index_elem);
      free(Bucket[d].index_l);	
    }
    d++;
    if ( d >= NUM_BUCKETS ) {
      printf("excedido el maximo numero de buckets: %d\n",NUM_BUCKETS);
      return 1;
    }
    buckets_empty = 1;
    for (i=0;i<=d;i++) {
      if (Bucket[d].num_elem != 0) {
	buckets_empty = 0;
	break;
      }
    }
    /* if all buckets are empty break the main bucle and finalize */
    if (buckets_empty == 1) break;
    if (d >= d_max) break;
  }


  printf("Distancia maxima alcanzada: %d\n",d); 
  
  free(Proto);
  free(Bucket[d].index_elem);
  free(Bucket[d].index_l);
  free(Element);
  printf("DT optimo3d OK\n");
  return d;
}
