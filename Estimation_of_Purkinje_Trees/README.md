
###**(c) Ruben Cardenes Almeida, 2014  Code to estimate Purkinje trees on LV surfaces** 

##Abstract:

The electrical activation of the heart is a complex physiological process that is essential for the understanding of several cardiac dysfunctions, such as ventricular tachycardia (VT). Nowadays, patient-specific activation times on ventricular chambers can be estimated from electro-anatomical maps, providing crucial information to clinicians for guiding cardiac radio-frequency ablation treatment. However, some relevant electrical pathways such as those of the Purkinje system are very difficult to interpret from these maps due to sparsity of data and the limited spatial resolution of the system. We present here a novel method to estimate these fast electrical pathways from the local activations maps (LATs) obtained from electro-anatomical maps. The location of Purkinje-myocardial junctions (PMJs) is estimated considering them as critical points of a distance map defined by the activation maps, and then minimal cost geodesic paths are computed on the ventricular surface between the detected junctions. Experiments to validate the proposed method have been carried out in simplified and realistic simulated data, showing good performance on recovering the main characteristics of simulated Purkinje networks (e.g. PMJs). A feasibility study with real cases of fascicular VT was also performed, showing promising results.

##First step: Singular points (PMJs) detection:  

Usually we start from a single vtk file, a surface with the LAT values on it (surface_LAT.vtk):  

Voxelize surface_LAT.vtk domain.vtk -sp 0.05
-sp 0.05 is the spacing we want for our volume. If the resulting volume is more than aprox 200x200x200 voxels

Increase the spacing to reduce the volume size, otherwise, the subsequent computations will require too much memory and time. 

###Now create an edge version of the volume as this:   

edgedetect3D domain.vtk edge.vtk -in_value 0 -out_value 255 


###Now manually remove the open part of the LV to create edge_open.vtk (this I should do it automatically but I was too lazy for that). Then, create the volumetric version of the LAT: LAT_volume.vtk   

ApplyScalarsToVolume surface_LAT.vtk edge.vtk LAT_volume.vtk -array_name Activation

We will use edge_open.vtk and LAT_volume.vtk in the following, so domain.vtk can be erased. 

### Now, we compute the singular points 
**INPUT:** 
edge_open.vtk  ----> a vtk volume (STRUCTURED_POINTS) with values = 1 on the LV surface = 0 inside, = 255 outside    
LAT_volume.vtk ----> a vtk volume (STRUCTURED_POINTS) with the LAT values on the LV surface  
**OUTPUT:** 
singular_points.vtk  a point set (POLYDATA)  with the detected points
**PARAMETERS:** 
-threshold (threshold used. In the paper we used Espsilon=threshold*180) 
-mode_source (to detect only source points, not sink points)  
#### Example: 

DetectGradientSingularities edge_open.vtk LAT_volume.vtk singular_pointsLAT.vtk -threshold 0.35 -mode_source

### After that, we group the points which were detected too close. The paremeter dist is related to the spacing of the data. When the spacing is 0.5, then I used dist =1  

`enter code here`GroupPoints singular_pointsLAT.vtk singular_pointsLAT_grouped.vtk -dist 1

## Second step: Streamlie construction:  

**INPUT:** 
edge_open.vtk  ----> a vtk volume (STRUCTURED_POINTS) with values = 1 on the LV surface = 0 inside, = 255 outside    
LAT_volume.vtk ----> a vtk volume (STRUCTURED_POINTS) with the LAT values on the LV surface  
singular_pointsLAT_grouped.vtk ----> a vtk points (POLYDATA) with singular points obtained in step 1  
-his_point HISpoint.vtk ----> HISpoints a vtk points (POLYDATA) with the HIS point manually selected (so, this file has to be created using for instance ParaView by selecting the point where the HIS is located).   
**OUTPUT:** 
streamlines.vtk  a point set (POLYDATA)  with the detected streamlines
**PARAMETERS:** 
-second_seed 0 (it says the algorithm to start joining from the first lowest value singular point detected)    
-alpha (parameter used in the method)  
-verbose (print information on screen) 
#### Example: (this one takes a bit more time) 

ComputePKTree edge_open.vtk LAT_volume.vtk singular_pointsLAT_grouped.vtk streamlines.vtk -verbose -second_seed 0 -alpha 15 -his_point HISpoint.vtk

### I use to smooth a little bit the resulting streamlines: 

SmoothPolyLine streamlines.vtk streamlines.vtk

### And sometimes, I reproject the points to the initial surface (surface_LAT): 

ProjectPointsToSurface streamlines.vtk surface_LAT.vtk streamlines_final.vtk

#### FIles Included: 
*Voxelize.cpp*
*edgedetect3D.cpp* 
*ApplyScalarsToVolume.cxx*
*DetectGradientSingularities.cpp*
*DToptimo3d.cpp*
*DToptimo3d.h*
*GroupPoints.cxx*
*ComputePKTree.cpp*
*SmoothPolyLine.cpp*
*ProjectPointsToSurface.cxx*
*vtkITKUtility.h*
> Written with [StackEdit](https://stackedit.io/).
