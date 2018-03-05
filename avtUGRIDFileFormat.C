/*****************************************************************************
*
* Copyright (c) 2000 - 2015, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-442911
* All rights reserved.
*
* This file is  part of VisIt. For  details, see https://visit.llnl.gov/.  The
* full copyright notice is contained in the file COPYRIGHT located at the root
* of the VisIt distribution or at http://www.llnl.gov/visit/copyright.html.
*
* Redistribution  and  use  in  source  and  binary  forms,  with  or  without
* modification, are permitted provided that the following conditions are met:
*
*  - Redistributions of  source code must  retain the above  copyright notice,
*    this list of conditions and the disclaimer below.
*  - Redistributions in binary form must reproduce the above copyright notice,
*    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
*    documentation and/or other materials provided with the distribution.
*  - Neither the name of  the LLNS/LLNL nor the names of  its contributors may
*    be used to endorse or promote products derived from this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
* ARE  DISCLAIMED. IN  NO EVENT  SHALL LAWRENCE  LIVERMORE NATIONAL  SECURITY,
* LLC, THE  U.S.  DEPARTMENT OF  ENERGY  OR  CONTRIBUTORS BE  LIABLE  FOR  ANY
* DIRECT,  INDIRECT,   INCIDENTAL,   SPECIAL,   EXEMPLARY,  OR   CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
* SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
* CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
* LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
* OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
* DAMAGE.
*
*****************************************************************************/

// ************************************************************************* //
//                            avtUGRIDFileFormat.C                           //
// ************************************************************************* //

#include <avtUGRIDFileFormat.h>

#include <FileFunctions.h>

#include <netcdf.h>

#include <DebugStream.h>

#include <algorithm> // sort, transform
#include <vector>

#include <vtkFloatArray.h>
#include <vtkRectilinearGrid.h>
#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkCellType.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#include <avtDatabaseMetaData.h>

#include <DBOptionsAttributes.h>
#include <Expression.h>

#include <InvalidVariableException.h>

#include <vtkCellDataToPointData.h>

#if defined(_WIN32)
 #include <win32-regex.h>
#else
 #include <regex.h>
#endif

#include <string>


// for some reason, there are two ways to read string attributes
// and it's difficult to predict which is correct.
std::string get_att_as_string(int ncid,int varid,const char *att_name) {
  int retval;
  char *att_data;
  std::string result("");

  if ( !(retval=nc_get_att_string(ncid, varid, att_name,&att_data)) ) {
    result=att_data; // is that legal?
    return result;
  }

  // debug5 << "nc_get_att_string failed, will trying get_att_text" << endl;

  // is it possible that we have to test for text or string??
  size_t att_len;
  if ( nc_inq_attlen(ncid,varid,att_name,&att_len) ) {
    debug5 << "nc_get_att_string _and_ attlen failed for " << att_name << endl;
    return result;
  }

  char *buff=new char[att_len+1];
  if ( nc_get_att_text(ncid,varid,att_name,buff ) ) {
    delete[] buff;
    return result;
  }

  buff[att_len]='\0';
  result=buff;
  delete[] buff;
  return result;
}


std::vector<std::string> split(const std::string src, char c = ' ')
{
  const char *str=src.c_str();

  std::vector<std::string> result;

  do {
    const char *begin = str;

    while(*str != c && *str)
      str++;
    
    result.push_back(std::string(begin, str));
  } while (0 != *str++);

  return result;
}

//////////---- MeshInfo ----//////////

MeshInfo::MeshInfo(int _ncid, int mesh_var,int z_var)
{
  char var_name[NC_MAX_NAME];
  char zvar_name[NC_MAX_NAME];
  ncid=_ncid;
  varid=mesh_var;
  layer_z_var=z_var;
  active_timestate=-1;

  nc_inq_varname(ncid,varid,var_name);
  if ( layer_z_var < 0 ) {
    name=var_name;
    layer_dim=-1;
  } else {
    nc_inq_varname(ncid,layer_z_var,zvar_name);
    name=std::string(var_name) + "." + std::string(zvar_name);
  }

  debug1 << "Creating mesh info " << name << endl;

  std::string cell_dim_name=get_att_as_string(ncid,varid,"face_dimension");
  if ( nc_inq_dimid(ncid,cell_dim_name.c_str(),&cell_dim) ) {
    debug1 << "Failed to read id of element dimension" << endl;
    return;
  }

  if ( nc_inq_dimlen(ncid,cell_dim,&n_cells2d) ) {
    debug1 << "Failed to read number of 2D cells" << endl;
    return;
  }

  debug1 << "Read n_cells="<< n_cells2d << endl;

  if(layer_z_var<0) {
    n_cells3d=n_cells2d;
    debug1 << "MeshInfo constructor: grid is 2D, so set n_cells3d=n_cells2d" << endl;
  }

  std::string node_coords=get_att_as_string(ncid,varid,"node_coordinates");
  if( node_coords=="" ) return;

  debug1 << "Node coordinate values " << node_coords << endl;

  std::string node_x=node_coords.substr(0,node_coords.find_first_of(" "));
  // debug1 << "node x variable '" << node_x << "'" << endl;
  std::string node_y=node_coords.substr(node_coords.find_last_of(" ")+1,
                                   node_coords.length());
  // debug1 << "node y variable '" << node_y << "'" << endl;

  if( nc_inq_varid(ncid,node_x.c_str(),&node_x_var) ) {
    debug1 << "Failed to find node_x variable" << endl;
    return;
  }
  if( nc_inq_varid(ncid,node_y.c_str(),&node_y_var) ) {
    debug1 << "Failed to find node_y variable" << endl;
    return;
  }

  // _assume_ that they have exactly one dimension.
  nc_inq_vardimid(ncid,node_x_var,&node_dim);
  nc_inq_dimlen(ncid,node_dim,&n_nodes);
  
  debug1 << "Found " << n_nodes << " nodes" << endl;

}

vtkPoints *MeshInfo::GetNodes(void) {
  float *xcoords=new float[n_nodes];
  float *ycoords=new float[n_nodes];

  if ( nc_get_var_float(ncid, node_x_var, xcoords) ) {
    debug1 << "Failed to read x coordinate" << endl;
    return NULL;
  }
  if ( nc_get_var_float(ncid, node_y_var, ycoords) ) {
    debug1 << "Failed to read y coordinate" << endl;
    return NULL;
  }
  
  debug1 << "Read coordinates" << endl;

  vtkPoints *points = vtkPoints::New();
  points->SetDataTypeToFloat();
  points->SetNumberOfPoints(n_nodes);
  
  float pnt[3];
  for(int n=0;n<n_nodes;n++){
    pnt[0]=xcoords[n];
    pnt[1]=ycoords[n];
    pnt[2]=0.0; // unused 
    points->SetPoint(n,pnt);
  }
  
  delete[] xcoords;
  delete[] ycoords;

  return points;
}

vtkDataSet *
MeshInfo::GetMesh(int timestate) 
{
  vtkUnstructuredGrid *surface=GetMesh2D(timestate);

  debug1 << "GetMesh..." << endl;
  if( layer_z_var < 0) {
    return surface;
  } else {
    debug1 << "...ExtrudeTo3D" << endl;
  
    vtkUnstructuredGrid *full=ExtrudeTo3D(timestate,surface);
    // without GetMesh2D calling mesh->Register(), this
    // crashes.  But that was because it was being Deleted
    // in multiple places.
    surface->Delete(); // is this okay? seems to make it crash.
    return full;
  }
}

vtkUnstructuredGrid *
MeshInfo::GetMesh2D(int timestate) 
{
  // the others are built on top of the node mesh
  vtkPoints *points = GetNodes();
  
  vtkUnstructuredGrid *mesh=vtkUnstructuredGrid::New();
  mesh->SetPoints(points);
  points->Delete() ; // pretty sure this is correct...
  mesh->Allocate();
  // maybe I shouldn't do this?
  // mesh->Register(NULL); // returning a pointer which owns itself

  // read the node_face info, build up triangles/quads.
  debug1 << "GetMesh: ugrid mesh name " << name << endl;

  nc_type xtype;
  
  std::string face_node=get_att_as_string(ncid,varid,"face_node_connectivity");

  debug1 << "Face-Node connectivity " << face_node << endl;

  int face_node_var;

  if( nc_inq_varid(ncid,face_node.c_str(),&face_node_var) ) {
    debug1 << "Failed to find face_node_connectivity variable" << endl;
    return NULL;
  }

  int f_n_dims[2]; // _assume_ face_node_connectivity has two dimensions, [Nfaces,MaxNodePerFace]
  nc_inq_vardimid(ncid,face_node_var,f_n_dims);

  int face_node_start=0;
  if ( nc_get_att_int(ncid,face_node_var,"start_index", &face_node_start) ) {
    debug1 << "Failed to find start_index for face_node_connectivity - will assume " 
           << face_node_start << endl;
  }

  int face_node_fill=-1; 
  if ( nc_get_att_int(ncid,face_node_var,"_FillValue", &face_node_fill) ) {
    debug1 << "Failed to find fill value for face_node_connectivity - will assume " 
           << face_node_fill << endl;
  }

  size_t max_node_per_face;
  // assumes that face_node_var is [cells,nodes]
  nc_inq_dimlen(ncid,f_n_dims[1],&max_node_per_face);
  
  int *faces=new int[n_cells2d*max_node_per_face];

  if ( nc_get_var_int(ncid, face_node_var, faces) ) {
    debug1 << "Failed to read face_node_var" << endl;
    return NULL;
  }

  debug1 << "UGRID: max node per face: " << max_node_per_face << endl;

  // Load as 2D, then optionally extrude to 3D
  vtkIdType *vertices=new vtkIdType[max_node_per_face];
    
  for(int f=0;f<n_cells2d;f++) {
    int n;
    for(n=0;n<max_node_per_face;n++) {
      if (faces[f*max_node_per_face+n] == face_node_fill ) 
        break;
      // last dimension varies fastest
      // convert to zero-based index
      // have to see if this is the right order
      vertices[n] = faces[f*max_node_per_face + n] - face_node_start; 
    }
    // This was failing on ERROR: bad cell type -- what fixed that?
    switch ( n ) {
    case 3:
      // if ( f==0 ) {
      //   debug1 << "UGRID:GetMesh: inserting triangle cell, vertices=" 
      //          << vertices[0] << " " << vertices[1] << " " << vertices[2] << endl;
      // }
      mesh->InsertNextCell(VTK_TRIANGLE,n,vertices);
      break;
    case 4: 
      mesh->InsertNextCell(VTK_QUAD,n,vertices);
      // debug1 << "UGRID: Trying to insert VTK_QUAD" << endl;
      break;
    default:
      // debug1 << "UGRID: Trying to insert VTK_POLYGON" << endl;
      mesh->InsertNextCell(VTK_POLYGON,n,vertices);
      break;
    }
  }
  delete[] vertices;
  delete[] faces;

  debug1 << "Returning 2D mesh" << endl;
  return mesh;
}

/**
   Dirty work of creating a prism with more than 6 nodes on top and bottom
   Ordering of the points making up faces is in the sense of a positive 
   outward facing normal.  i.e. a right-hand rule curling in the order of
   nodes gives the outward normal.  Equivalently, when viewed from outside
   the volume, nodes are ordered CCW.

   Assume that internally the nodes of a cell are stored in CCW order.  So
   the top face of the prism is already in the proper order.
 **/ 
void
insertNPrism(vtkUnstructuredGrid *full_mesh,
             int npoints2d,
             vtkIdType *point_ids) {
  // following http://www.vtk.org/Wiki/VTK/Examples/Cxx/GeometricObjects/Polyhedron
  vtkSmartPointer<vtkCellArray> faces =
    vtkSmartPointer<vtkCellArray>::New();

  vtkIdType *face=new vtkIdType[npoints2d];

  // create the top:
  for(int i=0;i<npoints2d;i++) {
    face[i]=point_ids[i]; // already proper order, I think.
  }
  faces->InsertNextCell(npoints2d, face);
  for(int i=0;i<npoints2d;i++) {
    // and the bottom - reverse the order
    face[i]=point_ids[2*npoints2d - 1 - i]; // old: npoints2d+i
  }
  faces->InsertNextCell(npoints2d, face);

  // and each side facet:
  for(int facet=0;facet<npoints2d;facet++) {
    face[3]=point_ids[facet];
    face[2]=point_ids[(facet+1)%npoints2d];
    face[1]=point_ids[npoints2d+(facet+1)%npoints2d];
    face[0]=point_ids[npoints2d+facet];
    faces->InsertNextCell(4,face);
  }

  full_mesh->InsertNextCell(VTK_POLYHEDRON, 2*npoints2d, point_ids,
                            2+npoints2d, faces->GetPointer());
  delete[] face;
}

void
MeshInfo::activateTimestate(int timestate) {
  // get some 3D mesh info for a given time step

  if ( (active_timestate == timestate) || 
       (layer_dim<0) ) 
    return;

  if ( cell_dim < 0 ) {
    debug1 << "cell_dim not set - don't know how to set timestate for 3D, nodal field" << endl;
    return;
  }

  n_cells3d=0; // to be incremented below
  
  if(cell_kmin.size() != n_cells2d )
    cell_kmin.resize(n_cells2d,-1);

  if(cell_kmax.size() != n_cells2d ) 
    cell_kmax.resize(n_cells2d,-1);

  // for transcribed delwaq output, typically have
  // volumes, which can be used to filter out dry segments.
  // for DFM output, it's probably sigma level and we can assume
  // all layers are active?

  // older DELWAQ code, need to figure out how to support both.

  // float *volumes=read_cell_z_full("volume",timestate);
  // if ( !volumes ) return;
  // for(int cell2d=0;cell2d<n_cells2d;cell2d++)  {
  //   // cell_kmin from first non-zero volume
  //   int k;
  //   for(k=0;
  //       (k<n_layers) && (volumes[cell2d*n_layers+k] <= 0);
  //       k++) ;
  // 
  //   cell_kmin[cell2d]=k;
  //   for(; (k<n_layers) && ( volumes[cell2d*n_layers+k]==volumes[cell2d*n_layers+k] );
  //       k++) ;
  //   cell_kmax[cell2d]=k;
  //   ncells_3d += (cell_kmax[cell2d] - cell_kmin[cell2d]);
  // }

  // new DFM code - sigma layers, so they're all full:
  for(int i=0;i<n_cells2d;i++){
    cell_kmax[i]=n_layers;
    cell_kmin[i]=0;
    n_cells3d += (cell_kmax[i] - cell_kmin[i]);
  }
  debug1 << "activateTimestate timestate=" << timestate 
         << " cell_kmax[0]=" << cell_kmax[0] 
         << " n_cells3d=" << n_cells3d << endl;

  active_timestate=timestate;
}

vtkDataArray *MeshInfo::ZoneToNode2D(vtkDataArray *zonal,vtkUnstructuredGrid *ds)
{
  // This is leaking
  // Taken from avtExpressionFilter.C, just the useful bits without
  // mamby pamby error checking.
  ds->GetCellData()->SetScalars(zonal);

  vtkCellDataToPointData *cd2pd = vtkCellDataToPointData::New();
  cd2pd->SetInputData(ds);
  // I think Update() is where the leak occurs, possibly leaking b/c
  // I was calling outv->Register(NULL) below?
  cd2pd->Update();
  vtkDataSet *ds3 = cd2pd->GetOutput();
  vtkDataArray *outv = ds3->GetPointData()->GetScalars();
  
  // This appears to be correct, but note that the caller must 
  // Delete() the returned data.
  outv->Register(NULL);

  cd2pd->Delete();

  return outv;
}


bool MeshInfo::set_bounds_from_variable(int bounds_var,
                                        std::vector<float> &layer_bottom,
                                        std::vector<float> &layer_top)
{
  float *layer_bounds=new float[2*n_layers];
  
  if ( nc_get_var_float(ncid,bounds_var,layer_bounds) ) {
    debug1 << "Failed to read layer bounds" << endl;
    bounds_var=-1; // fall through to next section
  } else {
    for(int i=0;i<n_layers;i++) {
      layer_bottom[i]=layer_bounds[2*i];
      layer_top[i]=layer_bounds[2*i+1];
    }
  }
  delete[] layer_bounds;

  return bounds_var>=0;
}

bool MeshInfo::set_bounds_from_dfm_flowelem_zw(std::vector<float> &layer_bottom,
                                               std::vector<float> &layer_top)
{
  int flowelem_zw_var;
  int k;

  if( nc_inq_varid(ncid,"FlowElem_zw",&flowelem_zw_var) ) {
    debug1 << "No FlowElem_zw, will not try DFM workaround" << endl;
    return false;
  }

  // in the netcdf map output FlowElem_zw has dimensions time, element, wdim
  size_t startp[3];
  size_t countp[3];

  startp[0]=0;
  countp[0]=1;
  startp[1]=0;
  countp[1]=n_cells2d;
  startp[2]=0;
  countp[2]=n_layers+1;
  
  float *elem_zw=new float[countp[0]*countp[1]*countp[2]];

  // scan full field, but only at first time step
  if ( nc_get_vara_float(ncid,flowelem_zw_var,startp,countp,elem_zw) ) {
    debug1 << "Failed to read FlowElem_zw" << endl;
    delete[] elem_zw;
    return false;
  }

  // looking for the deepest one
  float z;
  float z_min=1e10;
  for(int c=0;c<n_cells2d;c++) {
    z=elem_zw[c*(n_layers+1)];
    if (z>1e10) // missing value
      continue; 

    if( z >= z_min )
      continue;
    
    z_min=z; // this is the deepest so far

    // copy this one.  Go from the surface down, because DFM writes these
    // values padded on the surface end, whereas other variables and the mesh
    // are padded on the bed end of the array.
    int k_real=n_layers; // used for indexing layer_bottom and top
    // the loop k is indexing FlowElem_zw
    
    for(k=n_layers;k>=0;k--) {
      z=elem_zw[c*(n_layers+1)+k];
      if (z>1e10) { // watercolumn is missing surface cells. keep going deeper
        continue;
      }
      
      if(k_real<n_layers)
        layer_bottom[k_real]=z;
      if(k_real>0)
        layer_top[k_real-1]=z;
      k_real--;
    }
    
    // fill in bogus data below there:
    for(;k_real>=0;k_real--) {
      // fabricate unit-thickness layers from there down.
      z=layer_bottom[k_real+1] - 1;
      if(k_real<n_layers)
        layer_bottom[k_real]=z;
      if(k_real>0)
        layer_top[k_real-1]=z;
    }
  }

  for(k=0;k<n_layers;k++) {
    debug1 << " layer k="<<k<<" bottom="<< layer_bottom[k] << " top="<<layer_top[k]<<endl;
  }

  delete[] elem_zw;

  debug1 << "Set layer bounds based on FlowElem_zw" << endl;
  return true;
}


/* set_layer_bounds
 * based on UGRID layer coordinates and optional layer bounds,
 * populate elevations for top/bottom coordinates.  does not 
 * apply sigma transform.  Does not handle coordinates which differ
 * in space.
 */
void
MeshInfo::set_layer_bounds(std::string z_std_name,
                           std::vector<float> &layer_bottom,
                           std::vector<float> &layer_top) 
{
  // and the bounds variable? no guarantees that this exists.
  int bounds_var;
  std::string bounds_name=get_att_as_string(ncid,layer_z_var,"bounds");
  if ( bounds_name == "" ) {
    debug1 << "No bounds attribute" << endl;
    bounds_var=-1;
  } else {
    if( nc_inq_varid(ncid,bounds_name.c_str(),&bounds_var) ) {
      debug1 << "Failed to find bounds variable" << endl;
      bounds_var=-1;
    }
    debug5 << "Found vertical bounds variable" << endl;
  }

  // for now, assume one set of z coordinates for the whole grid
  // (for z layers, we're ignoring partial layers)
  // this is already implicitly assumed by looking only for z coordinates with
  // a single dimension.
  layer_bottom.resize(n_layers);
  layer_top.resize(n_layers);

  if ( bounds_var>= 0) {
    if ( set_bounds_from_variable(bounds_var,layer_bottom,layer_top) )
      return;
  }

  // DFM workaround, if FullGridOutput=1, and maybe just z-layers even then.
  if ( set_bounds_from_dfm_flowelem_zw(layer_bottom,layer_top) )
    return;
  
  //  bounds_var<0, or failed to use bounds_var
  // read layer centers, extrapolate to layer bounds
  float *layer_centers=new float[n_layers];

  if ( nc_get_var_float(ncid,layer_z_var,layer_centers) ) {
    debug1 << "Failed to read layer_z_var" << endl;
    debug1 << "Fabricating!" << endl;
    for(int i=0;i<n_layers;i++) {
      layer_centers[i]=(float)(i+0.5) / (float)n_layers;
    }
  }

  // Special DFlowFM workaround - layer coordinates are
  // incorrect after first output, at least in r50237
  if( fabs(layer_centers[0]) > 100000 ) {
    debug1 << "layer coordinate looks corrupted.  Will fabricate." << endl;
    for(int i=0;i<n_layers;i++) {
      layer_centers[i]=(float)(i+0.5) / (float)n_layers;
    }
  }

  if ( n_layers==1 ) {
    if ( z_std_name=="ocean_zlevel_coordinate" ) {
      // total punt - could get smarter and look for 
      // eta and bedlevel.
      layer_bottom[0]=layer_centers[0]-0.5;
      layer_top[0]=layer_centers[0]+0.5;
    } else {
      // reasonable guess for sigma coordinates
      layer_bottom[0]=0;
      layer_top[0]=1;
    }
  } else {
    for(int k=0;k<n_layers;k++) {
      debug1 << "k = " << k << " of n_layers = " << n_layers << endl;
      // Set the lower bound:
      if (k==0) {
        if( (z_std_name=="ocean_zlevel_coordinate") ) { 
          layer_bottom[k]=layer_centers[0] - 0.5*(layer_centers[1]-layer_centers[0]);
        } else {
          // sigma coordinate may be increasing or decreasing, independent of
          // the positive: up/down attribute.  This comes up with DWAQ output where
          // layers start at the surface and increasing layer dimensions goes towards
          // the bed.
          if( layer_centers[1] > layer_centers[0] ) {
            layer_bottom[k]=0; 
          } else {
            layer_bottom[k]=1;
          }
        }
      } else {
        // copy from previous layer
        layer_bottom[k]=layer_top[k-1];
      }

      // and the top of the layer
      if (k==n_layers-1) {
        if( (z_std_name=="ocean_zlevel_coordinate") ) { 
          layer_top[k]=layer_bottom[k] + (layer_centers[k]-layer_centers[k-1]);
        } else {
          // whatever the increasing/decreasing state of the sigma coordinate,
          // top of the top should be opposite bottom of the bottom.
          layer_top[k]=1-layer_bottom[0]; 
        }
      } else {
        layer_top[k]=0.5*(layer_centers[k]+layer_centers[k+1]);
      }
      debug1 << "Layer bounds: " << layer_bottom[k] << " to " << layer_top[k] << endl;
    }
  }
  delete[] layer_centers;
}


vtkUnstructuredGrid *
MeshInfo::ExtrudeTo3D(int timestate,
                      vtkUnstructuredGrid *surface)
{
  vtkDataArray *bedlevel=NULL;
  vtkDataArray *eta=NULL;

  ////// Now basically extrude the surface mesh to prisms
  // bounded by the cell_divisions elevations 
  activateTimestate(timestate);

  // the new points 
  // copy the z=0 points from the flat surface, and then repeat
  // for each z-value 
  int n_surf_points = surface->GetNumberOfPoints();
  int offset;

  // for z-grid:
  //   have a vertical dimension - in the current file, nFlowMesh_layers
  //   hopefully that has a bounds attribute, which points to a variable
  //   with dimensions nFlowMesh_layers,2, giving top and bottom elevation
  //   of the z-layers.
  
  // so far we only support layer_z_var having a single dimension.

  int layer_var_ndim=-1;
  nc_inq_varndims(ncid,layer_z_var,&layer_var_ndim);
  if ( layer_var_ndim != 1 ) {
    debug1 << "Whoa - layer_var_ndim must be 1, but it was " << layer_var_ndim << endl;
    return NULL;
  }  

  // z-coordinate or sigma terrain-following?
  std::string z_std_name=get_att_as_string(ncid,layer_z_var,"standard_name");
  if ( z_std_name == "" ) {
    debug1 << "No standard name on z coordinate variable.  Can't tell sigma vs. z" << endl;

    std::string z_long_name=get_att_as_string(ncid,layer_z_var,"long_name");
    if ( z_long_name.find("z layer")==0 ) {
      debug1 << "Long name suggests z coordinate" << endl;
      z_std_name="ocean_zlevel_coordinate";
    } else if( z_long_name.find("sigma layer")==0 ) {
      debug1 << "Long name suggests sigma coordinate" << endl;
      z_std_name="ocean_sigma_coordinate";
    } else {
      return NULL;
    }
  }
  
  if ( z_std_name=="ocean_sigma_coordinate" ) {
    debug1<< "Vertical coordinate indicates sigma coordinates" << endl;

    // which variables define the layers?
    std::string formula=get_att_as_string(ncid,layer_z_var,"formula_terms");
    if ( formula=="" ) {
      debug1 << "Sigma coordinate, but no formula_terms" << endl;
      sigma_eta=sigma_bedlevel="";
    } else {
      // something like:
      // LayCoord_cc:formula_terms = "sigma: LayCoord_cc eta: s1 bedlevel: flowelem_bl" ;
      std::vector<std::string> tokens = split(formula);
      for (int i=0;i<tokens.size();i++) {
        if( tokens[i] == "eta:" ) {
          i++;
          sigma_eta=tokens[i];
        } else if ( tokens[i]=="bedlevel:" ) {
          i++;
          sigma_bedlevel=tokens[i];
        } else if ( tokens[i]=="sigma:" ) {
          i++;
          sigma_sigma=tokens[i];
        } else {
          debug1 << "sigma formula terms not understood: " << tokens[i] << endl;
        }
      }
      debug1 << "Sigma formula terms, sigma is " << sigma_sigma 
             << ", eta is " << sigma_eta 
             << ", and bedlevel is " << sigma_bedlevel << endl;
      if ( (sigma_eta!="") && (sigma_bedlevel!="") ) {
        bedlevel=parent->GetVar(timestate,sigma_bedlevel.c_str());
        eta=parent->GetVar(timestate,sigma_eta.c_str());
        // sigma_sigma should aleady be loaded as layer_centers, based
        // on layer_z_var
        debug1 << "Fetched values for sigma formula terms" << endl;

        // going out on a limb here --
        vtkDataArray *eta_node=ZoneToNode2D(eta,surface);
        vtkDataArray *bed_node=ZoneToNode2D(bedlevel,surface);
        eta->Delete();
        bedlevel->Delete();
        eta=eta_node;
        bedlevel=bed_node;
      }
    }
  } else if (z_std_name=="ocean_zlevel_coordinate" ) {
    debug1<< "Vertical coordinate indicates z-level coordinates" << endl;
  }

  //-----  Extract vertical coordinate bounds -----
  std::vector<float> layer_bottom;
  std::vector<float> layer_top;
  set_layer_bounds(z_std_name,layer_bottom,layer_top);


  //----------- Create points ---------
  vtkPoints *all_points = vtkPoints::New();
  // Here, add +1 to n_layers because 10 layers have 11 unique
  // node elevations
  all_points->SetNumberOfPoints( (n_layers+1)*n_surf_points );
  float *ap_data = (float*)all_points->GetVoidPointer(0);
  double *surf_point;
  float *a_point;

  for(int surf_point_id=0;surf_point_id<n_surf_points;surf_point_id++) {
    surf_point = surface->GetPoint(surf_point_id);
    for(int k=0;k<n_layers+1;k++){
      // pointer to the point being defined
      a_point = ap_data + 3*(surf_point_id+k*n_surf_points);
      
      a_point[0] = surf_point[0];
      a_point[1] = surf_point[1];

      // assumes that bounds are continuous - bottom of one cell same
      // as top of the next
      if(k==0) {
        a_point[2] = layer_bottom[0];
      } else {
        a_point[2] = layer_top[k-1];
      }

      // apply the sigma transformation
      if ( (eta != NULL) && (bedlevel != NULL) ) {
        double z_bed=bedlevel->GetComponent(surf_point_id,0);
        double z_surf=eta->GetComponent(surf_point_id,0);

        a_point[2]= a_point[2]*(z_surf - z_bed) + z_bed;
      }
    }
  }

  debug1 << "Done creating 3D points" << endl;

  vtkUnstructuredGrid *full_mesh = vtkUnstructuredGrid::New();
   
  // Copy the cell structure, extruding through depth
  vtkIdType surf_cell_npoints;
  vtkIdType *surf_cell_point_ids;
  
  // 2*MAX_SIDES: top and bottom points, up to a hexagonal prism, or larger
  vtkIdType new_cell_point_ids[2*MAX_SIDES];
  
  // full_cell2valid provides a mapping of expected cell index to
  //  actual cell index;
  
  // settings['stairstep']=0 ->  The bottom of the bottom valid cell is defined
  //   by the voronoi center depth
  //                      =1 ->  Round to next z-level (deeper)
  // It's going to be a real pain to do partial depths, because it
  //  upsets the mapping of points to different depths (bottom points
  //   can't be reused, and the numbering gets all funky)
  
  // Fetch the depth values in order to evaluate the bottom cells
  // this is a little different for the ugrid stuff.
  //  probably have node-centered depth
  //  what's the CF convention say about this?
  //  could create the full mesh?
  //  probably it is supposed to be a [time,cell,layer] variable
  //  in the netcdf.
  //  there's the further consideration that really we want bounds, too,
  //  so it would be [time,cell,layer,d2]
  // regardless, this needs to be time-dependent.
  
  // the most immediate question is how to figure out how many 
  // cells are in each water column.
  
  // vtkFloatArray *cell_depths =(vtkFloatArray *)GetVarBathymetry();
  // float *bath_data = cell_depths->GetPointer(0);
  
  full_mesh->Allocate();
   
  int expected_cell_id=0;
  int real_cell_id=0;
     
  for(int surf_cell_id=0 ;
      surf_cell_id<n_cells2d ;
      surf_cell_id++ ) {
    surface->GetCellPoints(surf_cell_id,
                           surf_cell_npoints,surf_cell_point_ids);
  
    if( surf_cell_npoints > MAX_SIDES ) {
      debug1 << "Cell has too many points - truncating! " << surf_cell_npoints << " to " << MAX_SIDES << endl;
      surf_cell_npoints=MAX_SIDES;
    }
   
    int k; // index into vertical cells

    if( surf_cell_id == 0 ) {
      debug1 << " surf_cell_id=" << surf_cell_id
             << " cell_kmin=" << cell_kmin[surf_cell_id] 
             << " cell_kmax=" << cell_kmax[surf_cell_id] 
             << endl;
    }
    for(k=cell_kmin[surf_cell_id];
        k<cell_kmax[surf_cell_id];
        k++) {
      
      /// The Visit manual shows this as 0-1-2 being one triangular
      // face with outward-facing normal, and 3-4-5 being the other
      // end with inward-facing normal (where the normal is defined
      // as outward on CCW ordered vertices)
      // if this really angers Visit, may have to test ordering here
      // and reverse top/bottom if vertices are in wrong order.
  
      // but the general polyhedron wants entirely outward facing
      // normals.
      // and the illustration of thw VTK_HEXAHEDRON 
      // shows 0,1,2,3 having an inward facing normal.
      
      // at some point we might add new bottom points to allow
      // for partial cell depths, but for now just do stair-stepping
  
      // tricky loop modification to order the nodes in the right way
      // for the different cell types.
      if( surf_cell_npoints==3 || surf_cell_npoints>6 )
        offset=0;
      else
        offset=surf_cell_npoints;
      
      // the upper level:
      for (int vertex=0;vertex<surf_cell_npoints;vertex++) {
        new_cell_point_ids[vertex+offset] = k*n_surf_points + surf_cell_point_ids[vertex];
      }
  
      offset=surf_cell_npoints-offset;
      
      // the lower level:
      for (int vertex=0;vertex<surf_cell_npoints;vertex++) {
        new_cell_point_ids[vertex+offset] = (k+1)*n_surf_points + 
          surf_cell_point_ids[vertex];
      }
        
      // Insert the cell into the mesh.
      if ( surf_cell_npoints==3 ) {
        full_mesh->InsertNextCell(VTK_WEDGE, 6, new_cell_point_ids);
      } else if (surf_cell_npoints==4 ) {
        full_mesh->InsertNextCell(VTK_HEXAHEDRON, 8, new_cell_point_ids);
      } else if (surf_cell_npoints==5 ) {
        full_mesh->InsertNextCell(VTK_PENTAGONAL_PRISM, 10, new_cell_point_ids);
      } else if (surf_cell_npoints==6 ) {
        full_mesh->InsertNextCell(VTK_HEXAGONAL_PRISM, 12, new_cell_point_ids);
      } else {
        // debug1 << "Surface cell has " << surf_cell_npoints << " points, which is simply too many" << endl;
        insertNPrism(full_mesh,surf_cell_npoints,new_cell_point_ids);
        // full_mesh->InsertNextCell(VTK_HEXAGONAL_PRISM, 12, new_cell_point_ids);
      }
      
      // update the mapping of cell ids:
      // the order of the data file is still a bit unclear, so hopefully
      // this is the same as in the data file...
      // as long as cell_kmin,cell_kmax are correct, there's not a compelling
      // reason to additionally store this mapping, afaict.

      // int expected_id = surf_cell_id + k*n_cells2d;
      // full_cell2valid[ expected_id ] = real_cell_id;
      // real_cell_id++;
    }
  }
  debug1 << "Done constructing 3D cells" << endl;
  
  // hopefully it's okay to set the points here, *after* defining the cells
  full_mesh->SetPoints(all_points);
  all_points->Delete();
  // surface->Delete(); // caller does this now.
   
  // this causes a leak.
  // full_mesh->Register(NULL); 

  // this is the counterpart to outv->Register() call in ZoneToNode2D, and
  // successfully cleans up the leak without introducing invalid reads
  if(eta!=NULL) eta->Delete();
  if(bedlevel!=NULL) bedlevel->Delete();
  
  debug1 << "Return 3D mesh" << endl;
  return full_mesh;
}


//////////---- VarInfo ----//////////

VarInfo::VarInfo(std::string n) 
  : name(n) 
{
  init();
}

VarInfo::VarInfo(void)
{
  init();
}

void
VarInfo::init(void) 
{
  time_dimi=cell_dimi=layer_dimi=node_dimi=-1;
  var_id=-1;
  ncid=-1;
  pseudo=P_REAL;
}

VarInfo::VarInfo(const VarInfo &a)
{
  // this really shouldn't be necessary, but it's acting up...
  name=a.name;
  mesh_name=a.mesh_name;
  time_dimi=a.time_dimi;
  cell_dimi=a.cell_dimi;
  layer_dimi=a.layer_dimi;
  node_dimi=a.node_dimi;
  time_dim=a.time_dim;
  pseudo=a.pseudo;

  var_id=a.var_id;
  ncid=a.ncid;
}


// ****************************************************************************
//  Method: avtUGRIDSingle constructor
//
//  Programmer: rusty -- generated by xml2avt
//  Creation:   Sat Mar 5 18:27:45 PST 2016
//
// ****************************************************************************

avtUGRIDSingle::avtUGRIDSingle(const char *filename,int _domain)
  : avtMTSDFileFormat(&filename, 1)
{
  // INITIALIZE DATA MEMBERS
  int retval;
  domain=_domain;
  
  if ( nc_open(filename, NC_NOWRITE, &ncid) ) {
    debug1 << "Failed to open " << filename << endl;
    return;
  } 
  debug1 << "UGRID: opened " << filename << endl;

  if ( nc_inq_dimid(ncid,"time", &time_dim) ) {
    debug1 << "Failed to read time dimensions" << endl;
  }

  if ( nc_inq_varid(ncid, "time", &time_var) ) {
    debug1 << "Couldn't read time variable id" << endl;
    return;
  }

  cell_kmin=cell_kmax=NULL;
  active_timestate=-1;

  int nvars;
  if ( nc_inq_nvars(ncid,&nvars) ) {
    debug1 << "Failed. How could nc_inq_nvars fail??" << endl;
    return;
  }

  // Scan the variables, making a note of any cf_role='mesh_topology'
  // variables.
  default_ugrid_mesh=""; // defaults to first found

  for(int varid=0;varid<nvars;varid++) {
    std::string cf_role=get_att_as_string(ncid,varid,"cf_role");
    if( cf_role=="mesh_topology" ) {
      MeshInfo minfo(ncid,varid);
      if(default_ugrid_mesh=="") {
        default_ugrid_mesh=minfo.name;
      }
      minfo.parent=this;
      mesh_table[minfo.name]=minfo;
    }
  }
  if ( default_ugrid_mesh=="" ) {
    debug1 << "Couldn't find any mesh variables with cf_role=='mesh_topology'" << endl;
    return;
  }

  initialize_metadata();
}

avtUGRIDSingle::~avtUGRIDSingle() {
  delete[] cell_kmin;
  delete[] cell_kmax;
}

std::string avtUGRIDSingle::var_mesh2d(std::string varname) {
  int varid;
  nc_inq_varid(ncid,varname.c_str(),&varid);
  return var_mesh2d(varid);
}

std::string avtUGRIDSingle::var_mesh2d(int varid) {
  // any reason not to use this same convenience util?
  std::string result=get_att_as_string(ncid,varid,"mesh");

  if( result=="" ) {
    result=default_ugrid_mesh;
    debug1 << "No mesh attribute - using default" << endl;
  } else if( result != default_ugrid_mesh ) {
    debug1 << "Found another ugrid mesh.  Hold on tight." << endl;
  }
  return result;
}

// ****************************************************************************
//  Method: avtEMSTDFileFormat::GetNTimesteps
//
//  Purpose:
//      Tells the rest of the code how many timesteps there are in this file.
//
//  Programmer: rusty -- generated by xml2avt
//  Creation:   Sat Mar 5 18:27:45 PST 2016
//
// ****************************************************************************

int
avtUGRIDSingle::GetNTimesteps(void)
{
  size_t length;
  char recname[NC_MAX_NAME+1];
  nc_inq_dim(ncid, time_dim, recname, &length);

  return length;
}

void avtUGRIDSingle::GetCycles(std::vector<int> &cycles)
{
  int nsteps=GetNTimesteps();

  if ( nsteps == 0 ) {
    // put in a fake 0 timestep
    cycles.push_back(0);
  } else {
    for(int i=0;i<nsteps;i++) {
      cycles.push_back(i);
    }
  }
}


void avtUGRIDSingle::GetTimes(std::vector<double> &times)
{
  // use days to make it a bit friendlier
  size_t startp[1]; // safely assume that time is one-dimensional
  size_t countp[1];

  countp[0]=GetNTimesteps();
  startp[0]=0;

  double *result=new double[countp[0]];

  if( nc_get_vara_double(ncid,time_var,startp,countp,result) ) {
    debug1 << "Failed to read double array for time " << endl;
    delete[] result;
  }

  double to_days=1;

  std::string units=get_att_as_string(ncid,time_var,"units");

  if ( units.find("second")== 0) {
    to_days=1./86400;
  } else if ( units.find("minute")==0 ) {
    to_days=1./(24*60);
  } else if ( units.find("hour")==0) {
    to_days=1./24;
  } else if ( units.find("day")==0 ) {
    to_days=1.;
  } else {
    debug1 << "UGRID::GetTimes: couldn't understand units " << units << endl;
  }

  if ( countp[0] == 0 ) {
    // put in a fake 0 timestep
    times.push_back(0.0);
  } else {
    for(int i=0;i<countp[0];i++) {
      times.push_back(result[i]*to_days);
    }
  }
  delete[] result;
}



// ****************************************************************************
//  Method: avtUGRIDSingle::FreeUpResources
//
//  Purpose:
//      When VisIt is done focusing on a particular timestep, it asks that
//      timestep to free up any resources (memory, file descriptors) that
//      it has associated with it.  This method is the mechanism for doing
//      that.
//
//  Programmer: rusty -- generated by xml2avt
//  Creation:   Sat Mar 5 18:27:45 PST 2016
//
// ****************************************************************************

void
avtUGRIDSingle::FreeUpResources(void)
{
}


// ****************************************************************************
//  Method: avtUGRIDSingle::initialize_metadata
//
//  Purpose:
//      Scan input file, figure out variables, meshes, dimensions etc.
//      Gather all of the information required for PopulateDatabaseMetaData,
//      but that is a separate call.  This is called by the constructor
//      if the basics of the netcdf check out
//
//  Programmer: rusty -- generated by xml2avt
//  Creation:   Sat Mar 5 18:27:45 PST 2016
//
// ****************************************************************************
void avtUGRIDSingle::initialize_metadata(void) 
{
  std::string ugrid_mesh=default_ugrid_mesh;

  // Loop through variables in file, trying to match them to grids
  // and register with metadata
  int nvars;
  if ( nc_inq_nvars(ncid,&nvars) ) {
    debug1 << "Failed. How could nc_inq_nvars fail??" << endl;
    return;
  }

  debug1 << "UGRID: Scanning variable definitions" << endl;

  char var_scan[NC_MAX_NAME];
  // char dim_name[NC_MAX_NAME];

  for(int var_num=0;var_num<nvars;var_num++) {
    if ( nc_inq_varname(ncid,var_num,var_scan) ) {
      debug1 << "Failed to read name of variable at index " << var_num << endl;
      continue;
    }

    debug1 << "Variable: "<<var_scan << " [" << var_num << "]" <<endl;

    // read dimensions - if it has either cell/node and a layer dimension, assign to
    // 3D mesh.
    // if just cell/node, put it on the 2d mesh.
    
    VarInfo var_inf(var_scan);
    var_inf.var_id=var_num;

    if( nc_inq_varndims(ncid,var_num,&(var_inf.ndims)) || 
        nc_inq_vardimid(ncid,var_num,var_inf.dims) ) {
      debug1 << "Failed.  How did that happen?" << endl;
      continue;
    }

    var_inf.ncid=ncid;

    // in the future, this will create the mesh on the fly...
    if ( !setMeshInfo(var_inf) )
      continue;

    var_table[var_scan] = var_inf;

    debug1 << "initialize_metadata: var '" << var_inf.name << "' => mesh " 
           << var_inf.mesh_name << endl;
  }

  // For any meshes that were found along the way, create a .domain variable
  // There does not appear to be function in VisIt expressions which returns
  // the domain ID.  Add that as a pseudo-variable for each mesh.
  for (std::map<std::string,MeshInfo>::iterator it=mesh_table.begin();
       it!=mesh_table.end();
       ++it) {
    std::string name=it->second.name + ".domain";
      
    VarInfo var_inf(name);
    var_inf.pseudo=VarInfo::P_DOMAIN;
    var_inf.ncid=-1; // no netcdf variable associated
    var_inf.time_dimi=-1; // no time
    // copy some dimensions from the mesh
    // Need to rethink this
    if ( it->second.cell_dim >= 0)
      var_inf.cell_dimi=1;
    if ( it->second.node_dim >= 0)
      var_inf.node_dimi=1;
    if ( it->second.layer_dim >= 0)
      var_inf.layer_dimi=1;
    var_inf.mesh_name=it->second.name;
      
    var_table[name]=var_inf;
    debug1 << "initialize_metadata: added mesh domain var " << var_inf.name << endl;
  }
}

// ****************************************************************************
//  Method: avtUGRIDSingle::PopulateDatabaseMetaData
//
//  Purpose:
//      This database meta-data object is like a table of contents for the
//      file.  By populating it, you are telling the rest of VisIt what
//      information it can request from you.
//
//  Programmer: rusty -- generated by xml2avt
//  Creation:   Sat Mar 5 18:27:45 PST 2016
//
// ****************************************************************************


void
avtUGRIDSingle::PopulateDatabaseMetaData(avtDatabaseMetaData *md, int timeState)
{
  std::string svar_scan;
  int nblocks = 1;  // <-- this must be 1 for MTSD
  int block_origin = 0;
  int spatial_dimension = 2;
  int topological_dimension = 2;
  double *extents = NULL;

  for (std::map<std::string,VarInfo>::iterator it=var_table.begin();
       it!=var_table.end();
       ++it) {
    avtScalarMetaData *smd=new avtScalarMetaData;

    smd->name=it->first;

    VarInfo &var_inf=it->second;

    std::string units=get_att_as_string(ncid,var_inf.var_id,"units");
    if (units != "" ) {
      smd->hasUnits=true;
      smd->units=units;
    }
    
    smd->meshName=mesh_table[var_inf.mesh_name].name;

    if ( var_inf.cell_dimi>=0 ) 
      smd->centering=AVT_ZONECENT;
    else if ( var_inf.node_dimi>=0 )
      smd->centering=AVT_NODECENT;
    else {
      debug1 << "var looked okay, but has neither cell nor node centering" << endl;
      delete smd; smd=NULL;
      continue;
    }

    md->Add(smd);
  }
    
  // We should have all of the meshes now - any 2D meshes defined
  // by a cf_role='mesh_topology' attribute, plus the 3D meshes
  // inferred from variable dimensions
  avtMeshType mt = AVT_UNSTRUCTURED_MESH;
  
  for (std::map<std::string,MeshInfo>::iterator it=mesh_table.begin();
       it!=mesh_table.end();
       ++it) {
    if ( it->second.layer_dim >= 0 ) {
      spatial_dimension=topological_dimension=3;
    } else {
      spatial_dimension=topological_dimension=2;
    }

    AddMeshToMetaData(md, it->second.name, mt, extents, nblocks, block_origin,
                      spatial_dimension, topological_dimension);
    debug1 << "PopulateDatabaseMetadata: adding mesh " << it->second.name
           << " spatial_dimension=" << spatial_dimension << endl;
    // AddMeshToMetaData(md,it->second.name+".nodes",AVT_POINT_MESH,NULL,1,0,2,0);

  }
  
  // CODE TO ADD A VECTOR VARIABLE
  //
  // string mesh_for_this_var = meshname; // ??? -- could be multiple meshes
  // string varname = ...
  // int vector_dim = 2;
  //
  // AVT_NODECENT, AVT_ZONECENT, AVT_UNKNOWN_CENT
  // avtCentering cent = AVT_NODECENT;
  //
  //
  // Here's the call that tells the meta-data object that we have a var:
  //
  // AddVectorVarToMetaData(md, varname, mesh_for_this_var, cent,vector_dim);
  //
}


// With spatial_dim_names populated and sorted, figure out which 
// mesh to give this thing
bool
avtUGRIDSingle::setMeshInfo(VarInfo &var_inf)
{
  int z_var=-1;
  char dim_name[NC_MAX_NAME];
  MeshInfo &mesh=mesh_table[default_ugrid_mesh];

  // presumably a 3D field (not yet dealing with structured grids!)
  for(int d=0;d<var_inf.ndims;d++){
    // iterate over known 2D meshes to see if there is a match on
    // cell or node dimension:
    if( var_inf.dims[d] == time_dim ) {
      var_inf.time_dimi=d;
      var_inf.time_dim=time_dim;
      continue;
    }
    
    // read dimension name for some kludgy tests below
    if ( nc_inq_dimname(ncid,var_inf.dims[d],dim_name) ) {
      debug1 << "  Failed to read name of dimension" << endl;
      return false;
    }
    bool found_dim_match=false;
    
    for (std::map<std::string,MeshInfo>::iterator it=mesh_table.begin();
         it!=mesh_table.end();
         ++it) {
      if ( var_inf.dims[d] == it->second.cell_dim ) {
        var_inf.cell_dimi=d; 
        var_inf.mesh_name=it->second.name;
        found_dim_match=true;
        break;
      } else if (var_inf.dims[d] == it->second.node_dim ) {
        var_inf.node_dimi=d; 
        var_inf.mesh_name=it->second.name;
        found_dim_match=true;
        break;
      }
    }
    if( found_dim_match ) continue;
    
    // extra checks based on the name. KLUDGE for DFM
    if ( std::string(dim_name) == "nFlowElem" ) {
      var_inf.cell_dimi=d;
      var_inf.mesh_name=default_ugrid_mesh;
      continue;
    }
    
    // Is this a vertical layer / level dimension?
    // how do we infer what is a layer dimension?
    //   could scan the variables, looking for variables
    //   with positive:up or positive:down.
    // if those are 1-D variables, then it's a
    // vertical dimension.
    // that might be enough to get us started, and just have to
    // come back at some point and deal with the case of mixed
    // layer grids - where layer definition varies by 2D cell.
    if ( z_var < 0 ) {
      z_var = vertical_coordinate_for_dimension(var_inf.dims[d]);
      if( z_var >= 0 ) {
        // code below will handle replacing var_inf.mesh_name with
        // a 3D mesh
        var_inf.layer_dimi=d;
        continue;
      }
    }
    debug5 << "  Failed to find any sort of match for dimension "<< dim_name <<endl;
    return false;
  }

  // at this point, must have found a 2D mesh
  if ( var_inf.mesh_name=="" )  {
    debug1 << "  setMeshInfo(" << var_inf.name << "): no 2D grid identified" << endl;
    return false;
  }
  
  if ( z_var >= 0 ) {
    debug1<<"  Emerged, and variable appears to be 3D." << endl;
    debug1<<"  Current mesh name is " << var_inf.mesh_name << endl;
    var_inf.mesh_name=create_3d_mesh(var_inf.mesh_name,
                                     var_inf.dims[var_inf.layer_dimi],
                                     z_var);
  }

  debug1 << "End of setMeshInfo, and variable has mesh name " << var_inf.mesh_name << endl;
  return true;
}

std::string
avtUGRIDSingle::create_3d_mesh(std::string mesh2d,int z_dim,int z_var) {

  char z_var_name[NC_MAX_NAME];
  nc_inq_varname(ncid,z_var,z_var_name);

  std::string result=mesh2d+"."+std::string(z_var_name);

  debug1 << "Combining 2D mesh " << mesh2d << " with vertical coordinate " << z_var_name << endl;
  debug1 << "  Result is " << result << endl;

  mesh_table[result] = MeshInfo(ncid,mesh_table[mesh2d].varid,z_var);
  MeshInfo &mesh3d=mesh_table[result];
  mesh3d.parent=this;
  mesh3d.layer_dim = z_dim;
  mesh3d.layer_z_var = z_var;

  if( nc_inq_dimlen(ncid,z_dim,&(mesh3d.n_layers)) ) {
    debug1 << "  failed to find layer dimension.  will punt with n_layers=1" << endl;
    mesh3d.n_layers=1;
  }

  return result;
}


// evolving semantics...
// return a var_id which appears to be a vertical coordinate along
// the given dimension.
// at this point only deal with completely separate 1-D coordinates,
// so the return var_id must be a 1-D variable.
// return -1 if none is found.
int
avtUGRIDSingle::vertical_coordinate_for_dimension(int dim) 
{
  // scan variables, see if there is a corresponding vertical coordinate 
  // variable.
  char var_scan[NC_MAX_NAME];
  int nvars;

  if ( nc_inq_nvars(ncid,&nvars) ) {
    debug1 << "  Failed. How could nc_inq_nvars fail??" << endl;
    return -1;
  }

  char dim_name[NC_MAX_NAME];
  int ndims;
  int dims[MAX_DIMS];
  std::string positive;

  for(int var_num=0;var_num<nvars;var_num++) {
    debug5 << " Is variable " << var_num << " a vertical coordinate for dim " << dim << endl;
    if( nc_inq_varndims(ncid,var_num,&ndims) || 
        nc_inq_vardimid(ncid,var_num,dims) ) {
      debug1 << "  Failed.  How did that happen?" << endl;
      continue;
    }
    // does it have any common dimensions with the requested dimension?
    int d;
    for(d=0;d<ndims;d++) {
      if( dim==dims[d] ) {
        break;
      }
    }
    if( d==ndims ) 
      continue ; // no match

    positive = get_att_as_string(ncid,var_num,"positive");
    if ( positive == "" ) {
      debug5 << "  Variable does not have positive attribute" << endl;
      continue; // not a vertical coordiante
    }
    if ( ndims!=1 ) {
      debug1 << "  Variable is a vertical coordinate, but not 1-D.  Keep looking" << endl;
      continue;
    }
    debug1 << "  Found a 1-D vertical variable " << var_num << endl;
    return var_num;
  }
  return -1;
}

// ****************************************************************************
//  Method: avtUGRIDSingle::GetMesh
//
//  Purpose:
//      Gets the mesh associated with this file.  The mesh is returned as a
//      derived type of vtkDataSet (ie vtkRectilinearGrid, vtkStructuredGrid,
//      vtkUnstructuredGrid, etc).
//
//  Arguments:
//      timestate   The index of the timestate.  If GetNTimesteps returned
//                  'N' time steps, this is guaranteed to be between 0 and N-1.
//      meshname    The name of the mesh of interest.  This can be ignored if
//                  there is only one mesh.
//
//  Programmer: rusty -- generated by xml2avt
//  Creation:   Sat Mar 5 18:27:45 PST 2016
//
// ****************************************************************************


vtkPoints *avtUGRIDSingle::GetNodes(const std::string ugrid_mesh) 
{
  return mesh_table[ugrid_mesh].GetNodes();
}


vtkUnstructuredGrid *avtUGRIDSingle::GetMeshNodes(const std::string ugrid_mesh) 
{
  vtkPoints *points=GetNodes(ugrid_mesh);

  int n_nodes=points->GetNumberOfPoints();

  vtkIdType onevertex;

  vtkUnstructuredGrid *ugrid = vtkUnstructuredGrid::New();
  ugrid->SetPoints(points);
  points->Delete();

  for(int n=0;n<n_nodes;n++) {
    onevertex = n;
    ugrid->InsertNextCell(VTK_VERTEX,1,&onevertex);
  }

  return ugrid;
}

vtkDataSet *
avtUGRIDSingle::GetMesh(int timestate, const char *meshname)
{
  std::string requested(meshname);
  std::string ugrid_mesh;
  int ndim;
  int retval;

  // meh - I guess this goes here
  activateTimestate(timestate);

  // this should probably be removed in favor of routing everything through
  // the mesh_table:
  if ( requested.find(".nodes") != std::string::npos) {
    ugrid_mesh = requested.substr(0,requested.length()-6);
    debug1 << "UGRID:GetMesh: requested '"<< requested << "' - will hand to GetMeshNodes" << endl;
    return GetMeshNodes(ugrid_mesh);
  } 

  // if(  requested.find(".2d",requested.length()-3) != std::string::npos ) {
  // ugrid_mesh = requested.substr(0,requested.length()-3);
  // ndim=2;
  if ( mesh_table.find(meshname) != mesh_table.end() ) {
    debug1 << "Found mesh name in mesh_table" << endl;
    return mesh_table[meshname].GetMesh(timestate);
  } else {
    debug1 << "In mesh name '"<<meshname<<"' couldn't decipher ugrid mesh" << endl;
    return NULL;
  }
  // code that was here moved to MeshInfo::GetMesh()
}



void avtUGRIDSingle::activateTimestate(int timestate) {
  if (timestate == active_timestate )
    return;
}


float * VarInfo::read_cell_at_time(int timestate, MeshInfo &mesh)
{
  size_t n_cells = mesh.n_cells3d; // or should that be n_cells2d??
  debug1 << "read_cell_at_time: allocating " << n_cells << endl;

  float * result=new float[n_cells];
  size_t startp[2];
  size_t countp[2];
 
  // assume that the dimensions are time,cell,layer
  if( time_dimi >= 0 ){
    startp[time_dimi] = timestate;
    countp[time_dimi]=1;
  }
  startp[cell_dimi] = 0;
  countp[cell_dimi]=n_cells;

  if( nc_get_vara_float(ncid,var_id,startp,countp,result) ) {
    debug1 << "Failed to read float array for " << name << endl;
    debug1 << "startp " << startp[0] << "," << startp[1] << endl;
    debug1 << "countp " << countp[0] << "," << countp[1] << endl;
    debug1 << "time_dimi " << time_dimi << endl;
    debug1 << "cell_dimi " << cell_dimi << endl;
    
    delete[] result;
    return NULL;
  }
  debug1 << "Finishing read_cell_at_time of " << name << endl;
  return result;
}

float * VarInfo::read_node_at_time(int timestate, MeshInfo &mesh)
{
  size_t n_nodes = mesh.n_nodes;
  debug1 << "read_node_at_time: allocating " << n_nodes << endl;

  float * result=new float[n_nodes];
  size_t startp[2];
  size_t countp[2];
 
  // assume that the dimensions are time,cell,layer
  if( time_dimi >= 0 ){
    startp[time_dimi] = timestate;
    countp[time_dimi]=1;
  }
  startp[node_dimi] = 0;
  countp[node_dimi]=n_nodes;

  if( nc_get_vara_float(ncid,var_id,startp,countp,result) ) {
    debug1 << "Failed to read float array for " << name << endl;
    debug1 << "startp " << startp[0] << "," << startp[1] << endl;
    debug1 << "countp " << countp[0] << "," << countp[1] << endl;
    debug1 << "time_dimi " << time_dimi << endl;
    debug1 << "node_dimi " << node_dimi << endl;
    
    delete[] result;
    return NULL;
  }
  debug1 << "Finishing read_node_at_time of " << name << endl;
  return result;
}

float *VarInfo::read_cell_z_at_time(int timestate, MeshInfo &mesh) {
  size_t n_cells,n_layers;

  if ( nc_inq_dimlen(ncid,mesh.cell_dim,&n_cells) ||
       nc_inq_dimlen(ncid,mesh.layer_dim,&n_layers) ) {
    debug1 << "Somehow late read of n_cells or n_layers failed for var " << name << endl;
    return NULL;
  }

  debug1 << "read_cell_z_at_time: allocating " << n_cells << "x" << n_layers << endl;

  float * result=new float[n_cells*n_layers];
  size_t startp[3];
  size_t countp[3];

  // assume that the dimensions are time,cell,layer
  if( time_dimi>=0 ) {
    startp[time_dimi] = timestate;
    countp[time_dimi]=1;
  } 

  debug1 << "read_cell_z_at_time: time_dimi="<< time_dimi << endl;
  debug1 << "                     cell_dimi="<< cell_dimi << endl;
  debug1 << "                    layer_dimi="<< layer_dimi << endl;

  startp[cell_dimi] = 0;
  countp[cell_dimi]=n_cells;

  startp[layer_dimi]=0;
  countp[layer_dimi]=n_layers;

  for(int i=0;i<3;i++){
    debug1 << "dim " << i << ": start=" << startp[i] << "  count="<<countp[i] << endl;
  }

  char var_scan[NC_MAX_NAME];
  if ( nc_inq_varname(ncid,var_id,var_scan) ) {
    debug1 << "Failed to read varname for var_id=" << var_id << endl;
    delete[] result;
    return NULL;
  } else {
    debug1 << "Expected var name " << name << " and queried to get " << var_scan << endl;
  }

  // this is failing, but not clear why...
  // maybe double check to make sure that var_id is correct?
  // could also print startp and countp, and query number of dimensions?

  if( nc_get_vara_float(ncid,var_id,startp,countp,result) ) {
    debug1 << "Failed to read 3D float array for " << name << endl;
    delete[] result;
    return NULL;
  } else {
    debug1 << "Successfully read 3D float array for " << name << endl;
  }

  // May have to transpose result:
  if ( cell_dimi > layer_dimi ) {
    debug1 << "In-memory transpose as dimensions are (layer,cell)" << endl;
    float * trans_result=new float[n_cells*n_layers];
    for(int cell=0;cell<n_cells;cell++) {
      for(int layer=0;layer<n_layers;layer++) {
        // in result, the fastest changing index is cell
        // in trans_result, the fastest changing index is layer
        trans_result[cell*n_layers + layer] = result[layer*n_cells + cell];
      }
    }
    delete[] result;
    result=trans_result;
    debug1 << "Done with transpose" << endl;
  }
  
  return result;
}

// Read a 3D netcdf variable by name and timestate.  
// Caller is responsible for freeing the result
float *
avtUGRIDSingle::read_cell_z_full(std::string varname,int timestate) 
{
  VarInfo &var=var_table[varname];
  return var.read_cell_z_at_time(timestate,mesh_table[var.mesh_name]);
}

// ****************************************************************************
//  Method: avtUGRIDSingle::GetVar
//
//  Purpose:
//      Gets a scalar variable associated with this file.  Although VTK has
//      support for many different types, the best bet is vtkFloatArray, since
//      that is supported everywhere through VisIt.
//
//  Arguments:
//      timestate  The index of the timestate.  If GetNTimesteps returned
//                 'N' time steps, this is guaranteed to be between 0 and N-1.
//      varname    The name of the variable requested.
//
//  Programmer: rusty -- generated by xml2avt
//  Creation:   Sat Mar 5 18:27:45 PST 2016
//
// ****************************************************************************

vtkDataArray *
avtUGRIDSingle::GetVar(int timestate, const char *varname)
{
  std::string svarname(varname);
  VarInfo vi=var_table[svarname];

  if( vi.name == "" ) {
    debug1 << "Case-mismatch in variable names?" << endl;
    
    // Iterate and look for case-insensitive matches
    std::map<std::string,VarInfo>::iterator it = var_table.begin();
 
    // Iterate over the map using Iterator till end.
    while (it != var_table.end()) {
      // Accessing KEY from element pointed by it.
      std::string key=std::string(it->first);
      std::transform(key.begin(),key.end(),key.begin(),::tolower);
      std::transform(svarname.begin(),svarname.end(),svarname.begin(),::tolower);

      if ( key==svarname ) {
        debug1 << "Found case-insensitive match to " << varname << endl;
        vi=it->second;
        break;
      }
      it++; // move to next entry
    }
  }

  if( vi.name == "" ) {
    debug1 << "BAD: appears that varname " << varname << " is not in the table" << endl;
    return NULL;
  }

  if ( vi.layer_dimi >= 0 ) {
    debug1 << "GetVar(" << timestate << "," << varname << ") => GetVar3D" << endl;
    return GetVar3D(timestate,vi);
  } else {
    debug1 << "GetVar(" << timestate << "," << varname << ") => GetVar2D" << endl;
    // debug1 << "in GetVar, vi.name is " << vi.name << endl;
    // // maybe need to make it into a string??
    // debug1 << "but direct from table? " << var_table[varname].name << endl;
    return GetVar2D(timestate,vi);
  }
}

vtkDataArray *
avtUGRIDSingle::GetVar3D(int timestate,VarInfo &vi)
{ 
  activateTimestate(timestate); // this updates n_cells3d if needed

  MeshInfo &mesh=mesh_table[vi.mesh_name];

  float *full;

  vtkFloatArray *rv = vtkFloatArray::New();

  debug1 << "GetVar3D(" << timestate << "," << vi.name << ")" << endl;
  debug1 << "GetVar3D: var.mesh_name is " << vi.mesh_name << endl;
  debug1 << "GetVar3D: mesh is " << mesh.name << endl;
  debug1 << "GetVar3D: pseudo is " << vi.pseudo << endl;
  debug1 << "GetVar3D: P_DOMAIN " << VarInfo::P_DOMAIN << endl;

  if ( vi.cell_dimi>=0 ) {
    rv->SetNumberOfTuples(mesh.n_cells3d);
    rv->SetNumberOfComponents(1);

    debug1 << " allocating " << mesh.n_cells3d << endl;

    if ( vi.pseudo==VarInfo::P_DOMAIN ) {
      for(int i=0;i<mesh.n_cells3d;i++) {
        rv->SetTuple1(i,(float)domain);
      }
    } else {
      full=vi.read_cell_z_at_time(timestate,mesh);
      // almost certainly some quick shortcut for this.
      for(int i=0;i<mesh.n_cells3d;i++) {
        rv->SetTuple1(i,full[i]);
      }
      delete[] full;
    }
  } else if ( vi.node_dimi>=0 ) {
    debug1 << "Not ready for 3D node variables" << endl;
    return NULL;
  }

  debug1 << "Done with GetVar3D for " << vi.name << endl;

  return rv;
}

vtkDataArray *
avtUGRIDSingle::GetVar2D(int timestate,VarInfo &vi)
{ 
  float *full;

  activateTimestate(timestate);
  vtkFloatArray *rv = vtkFloatArray::New();

  MeshInfo &mesh=mesh_table[vi.mesh_name];

  debug1 << "GetVar2D(" << timestate << "," << vi.name << ") allocating " << mesh.n_cells2d << endl;
  debug1 << "GetVar2D: var.mesh_name is " << vi.mesh_name << endl;
  debug1 << "GetVar2D: mesh is " << mesh.name << endl;

  if ( vi.cell_dimi>=0 ) {
    rv->SetNumberOfTuples(mesh.n_cells2d);
    rv->SetNumberOfComponents(1);

    if ( vi.pseudo==VarInfo::P_DOMAIN ) {
      for(int i=0;i<mesh.n_cells2d;i++) {
        rv->SetTuple1(i,domain);
      }
    } else {
      full=vi.read_cell_at_time(timestate,mesh);
      if( full==NULL) 
        return rv;

      for(int i=0;i<mesh.n_cells2d;i++) {
        rv->SetTuple1(i,full[i]);
      }
      delete[] full;
    }
  } else if ( vi.node_dimi>=0 ) {
    rv->SetNumberOfTuples(mesh.n_nodes);
    rv->SetNumberOfComponents(1);

    full=vi.read_node_at_time(timestate,mesh);
    if(full==NULL)
      return rv;

    for(int i=0;i<mesh.n_nodes;i++) {
      rv->SetTuple1(i,full[i]);
    }
    delete[] full;
  }

  debug1 << "Done with GetVar2D for " << vi.name << endl;
  return rv;
}


// ****************************************************************************
//  Method: avtUGRIDSingle::GetVectorVar
//
//  Purpose:
//      Gets a vector variable associated with this file.  Although VTK has
//      support for many different types, the best bet is vtkFloatArray, since
//      that is supported everywhere through VisIt.
//
//  Arguments:
//      timestate  The index of the timestate.  If GetNTimesteps returned
//                 'N' time steps, this is guaranteed to be between 0 and N-1.
//      varname    The name of the variable requested.
//
//  Programmer: rusty -- generated by xml2avt
//  Creation:   Sat Mar 5 18:27:45 PST 2016
//
// ****************************************************************************

vtkDataArray *
avtUGRIDSingle::GetVectorVar(int timestate, const char *varname)
{
  // YOU MUST IMPLEMENT THIS
  return NULL;
  
    //
    // If you have a file format where variables don't apply (for example a
    // strictly polygonal format like the STL (Stereo Lithography) format,
    // then uncomment the code below.
    //
    // EXCEPTION1(InvalidVariableException, varname);
    //

    //
    // If you do have a vector variable, here is some code that may be helpful.
    //
    // int ncomps = YYY;  // This is the rank of the vector - typically 2 or 3.
    // int ntuples = XXX; // this is the number of entries in the variable.
    // vtkFloatArray *rv = vtkFloatArray::New();
    // int ucomps = (ncomps == 2 ? 3 : ncomps);
    // rv->SetNumberOfComponents(ucomps);
    // rv->SetNumberOfTuples(ntuples);
    // float *one_entry = new float[ucomps];
    // for (int i = 0 ; i < ntuples ; i++)
    // {
    //      int j;
    //      for (j = 0 ; j < ncomps ; j++)
    //           one_entry[j] = ...
    //      for (j = ncomps ; j < ucomps ; j++)
    //           one_entry[j] = 0.;
    //      rv->SetTuple(i, one_entry); 
    // }
    //
    // delete [] one_entry;
    // return rv;
    //
}




// Multi-domain support
avtUGRIDFileFormat::avtUGRIDFileFormat(const char *filename)
  : avtMTMDFileFormat(filename)
{
  // placeholder code for minimal wrapper test
  populate_filenames(filename);
  domain_count=filenames.size();
  
  domain_cache.resize(domain_count,NULL);
  metadata_cache.resize(domain_count,NULL);
}

// Some machinery for reading directories the visit way...
void save_file_name(void *data, const std::string &filename, bool a, bool b, long c) {
  std::vector<std::string> *files;
  files = (std::vector<std::string>*)data;
  files->push_back(filename);

  debug1 << "Found file: " << filename << endl;
}

// Based on the path provided by the GUI, look for similar
// filenames to make up a multi-processor run. Sets
// filenames[i].
void
avtUGRIDFileFormat::populate_filenames(const char *filename)
{
  std::vector<std::string> all_files; // everything in the directory

  std::string my_path = FileFunctions::Dirname(filename);
  std::string basename=FileFunctions::Basename(filename);

  // the specified file is always first
  filenames.push_back(filename);

  // Somehow this chunk of code is problematic:

  // The quest for subdomains:
  // trouble linking this one:
  // exact prototype is (include/visit/common/misc/FileFunctions.h)
  //   bool        MISC_API ReadAndProcessDirectory(const std::string &,
  //                                         ProcessDirectoryCallback *,
  //                                         void * = 0,
  //                                         bool = false);
  
       
  FileFunctions::ReadAndProcessDirectory( my_path, (FileFunctions::ProcessDirectoryCallback*)save_file_name, (void*)&all_files, false );
  debug1 << "my_path: " << my_path << endl;
  debug1 << "base: " << basename << endl;
  // now all of the files, including ., .., are in the filenames vector.
  // those are full, absolute, paths.
  // my_path has the directory, without trailing /
  // base has our basename.

  // match against a hopeful pattern of multidomain output
  regex_t cre;
  regmatch_t pm[2];
  regmatch_t pm_sub[2]; // for matching to other subdomains

  // need to match both
  // This is a single output file, prefixed only by the name of the run.
  // short_test_08_0000_map.nc
  // and
  // This is the name of the run, the proce, but then a timestamp
  // wy2013c_0000_20120801_000000_map.nc
  // The greediness is making that more difficult.
  if ( regcomp(&cre, ".*_([0-9][0-9][0-9][0-9])_(.*_)?map\\.nc$", REG_EXTENDED)) {
    debug1 << "Couldn't compile pattern for finding raw data files" << endl;
    return;
  } else {
    debug1 << "Compiled the regex" << endl;
  }

  // does the basename match?
  if ( regexec(&cre, basename.c_str(), 2, pm, 0) != 0 ) {
    regfree(&cre);
    return;
  }

  debug1 << "Matched, pm[1] " << pm[1].rm_so << " to " << pm[1].rm_eo << endl;

  std::string proc_str=std::string(basename.c_str()+pm[1].rm_so,
                                   pm[1].rm_eo - pm[1].rm_so);
  int proc_num=atoi(proc_str.c_str());

  // only look for subdomains if it looks like we were given domain 0.
  if( proc_num!=0 ) {
    regfree(&cre);
    return;
  }

  debug1 << "Looks like proc 0: " << basename << endl;

  debug1 << "About to loop through candidate filenames" << endl;
  for(int proc=1;proc<MAX_SUBDOMAINS;proc++) {
    int file_idx=0;
    for(;file_idx<all_files.size();file_idx++) {
      std::string basename_sub=FileFunctions::Basename(all_files[file_idx]);

      if ( regexec(&cre, basename_sub.c_str(), 2, pm_sub, 0) != 0 ) 
        continue;

      // only okay if the everything except the domain is the same.
      if( basename.substr(0,pm[1].rm_so) != 
          basename_sub.substr(0,pm_sub[1].rm_so) )
        continue;
      if( basename.substr(pm[1].rm_eo,std::string::npos) != 
          basename_sub.substr(pm_sub[1].rm_eo,std::string::npos) )
        continue;

      debug5 << "Maybe " << basename_sub << endl;
      
      std::string proc_str=std::string(basename_sub.c_str()+pm[1].rm_so,
                                       pm[1].rm_eo - pm[1].rm_so);
      int proc_num=atoi(basename_sub.substr(pm_sub[1].rm_so,
                                            pm_sub[1].rm_eo).c_str());
      if( proc_num!= proc ) 
        continue;

      debug1 << "Yep! " << basename_sub << endl;
      filenames.push_back(all_files[file_idx]);
      break;
    }
    if( file_idx==all_files.size() ) {
      regfree(&cre);
      return;
    }
  }
  regfree(&cre);
}

// How to handle avtDatabaseMetaData for top-level vs. 
// subdomains?
// Any tricks to allocating one of those ourselves, such
// that the toplevel instance can get each subdomain to
// create a metadata object, which takes care of their
// internal initialization, then maybe copy/modify one of
// those to handle the toplevel call to get metadata?

avtUGRIDSingle *avtUGRIDFileFormat::subdomain(int domain) {
  if( domain_cache[domain] == NULL ) {
    debug1 << "Allocating new UGRIDSingle instance" << endl;

    // placeholder handling of filepath
    domain_cache[domain] = new avtUGRIDSingle(filenames[domain].c_str(),domain);
    // this will initialize some internal state, too.
    metadata_cache[domain] = new avtDatabaseMetaData();
    // punt with timestate=0
    domain_cache[domain]->PopulateDatabaseMetaData(metadata_cache[domain],0);
  }
  return domain_cache[domain];
}

int
avtUGRIDFileFormat::GetNTimesteps(void)
{
  return subdomain(0)->GetNTimesteps();
}

void
avtUGRIDFileFormat::GetCycles(std::vector<int> &cycles) 
{
  subdomain(0)->GetCycles(cycles);
}

void
avtUGRIDFileFormat::GetTimes(std::vector<double> &times)
{
  subdomain(0)->GetTimes(times);
}

void
avtUGRIDFileFormat::FreeUpResources(void)
{
}

vtkDataSet *
avtUGRIDFileFormat::GetMesh(int timestate, int domain, const char *meshname)
{
  debug1 << "UGRID: request for domain " << domain << " mesh " << meshname << endl;
  return subdomain(domain)->GetMesh(timestate,meshname);
}

vtkDataArray *
avtUGRIDFileFormat::GetVar(int timestate, int domain, const char *varname)
{
  return subdomain(domain)->GetVar(timestate,varname);
}

vtkDataArray *
avtUGRIDFileFormat::GetVectorVar(int timestate, int domain,const char *varname)
{
  return subdomain(domain)->GetVectorVar(timestate,varname);
}


void
avtUGRIDFileFormat::PopulateDatabaseMetaData(avtDatabaseMetaData *md, int timeState)
{
  // these seem to work only when debug is 5
  // even with some updates to PopulateDatabaseMetaData
  if( 0 ) {
    subdomain(0); // force allocation and loading metadata
    avtDatabaseMetaData *sub_meta=metadata_cache[0];
   
    // Should copy info over:
    *md = *sub_meta;
  } else {
    // state in the subdomain, and duplicate some work.
    subdomain(0)->PopulateDatabaseMetaData(md,timeState);
  }

  // That was leading to a "invalid variable name Mesh2D" error.
  // maybe the copy doesn't really work?
  // This works only when debug level is high.  It fails for debug=1, okay
  // for 5.

  debug1 << "UGRID: overwriting domain count in metadata" << endl;
  for(int mesh_num=0;mesh_num<md->GetNumMeshes();mesh_num++) {
    md->SetBlocksForMesh(mesh_num,domain_count);
  }
}

