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
//                            avtUGRIDFileFormat.h                           //
// ************************************************************************* //

#ifndef AVT_UGRID_FILE_FORMAT_H
#define AVT_UGRID_FILE_FORMAT_H

#include <avtMTSDFileFormat.h>
#include <avtMTMDFileFormat.h>
#include <vtkUnstructuredGrid.h>

#include <vector>
#include <map>

#define MAX_SIDES 50
#define MAX_DIMS 6
#define MAX_SUBDOMAINS 256

// forward declarations:
class avtUGRIDFileFormat;
class avtUGRIDSingle;

// corresponds 1:1 (I think) with meshes in the eyes of visit.  
// so whether a variable is on the nodes or cells is up to the variable.
class MeshInfo {
public:
  std::string name;

  // formula terms in case of sigma layers:
  std::string sigma_sigma,sigma_eta,sigma_bedlevel;

  int ncid,varid;
  int active_timestate;
  std::vector<int> cell_kmin,cell_kmax;
  avtUGRIDSingle *parent;

  // dimension ids - set to -1 if this mesh doesn't have those
  // dimensions
  int node_dim, cell_dim, layer_dim;
  int node_x_var,node_y_var, layer_z_var;
  size_t n_nodes, n_cells2d, n_layers, n_cells3d;

  MeshInfo(int ncid,int varid,int z_var=-1);
  MeshInfo() {ncid=-1 ; varid=-1; layer_z_var=-1; active_timestate=-1; parent=NULL; };
  vtkPoints *GetNodes(void);
  vtkDataSet *GetMesh(int timestate);
  
  vtkUnstructuredGrid *GetMesh2D(int timestate);

  vtkUnstructuredGrid *ExtrudeTo3D(int timestate,vtkUnstructuredGrid *surface);
  void set_layer_bounds(std::string z_std_name,
                        std::vector<float> &layer_bottom,
                        std::vector<float> &layer_top);
  bool set_bounds_from_variable(int bounds_var,
                                std::vector<float> &layer_bottom,
                                std::vector<float> &layer_top);
  bool set_bounds_from_dfm_flowelem_zw(std::vector<float> &layer_bottom,
                                       std::vector<float> &layer_top);

  void activateTimestate(int);

  vtkDataArray *ZoneToNode2D(vtkDataArray *,vtkUnstructuredGrid *);
};

class VarInfo {
public:
  std::string name;
  std::string mesh_name;
  std::vector<std::string> spatial_dim_names;
  int ndims;
  int dims[MAX_DIMS];

  // these are -1 if not present, and otherwise
  // give where the respective dimension falls in the variable definition.
  // node_dimi vs. cell_dimi control the staggering of the variable.
  int time_dimi,cell_dimi,layer_dimi,node_dimi;
  // cell, layer, node are stored in a mesh object - shared among variables
  // which live on the same grid.  time_dim is trickier - for a static grid,
  // doesn't necessarily need to be shared.
  int time_dim; // ,cell_dim,layer_dim,node_dim;
  int var_id,ncid;

  VarInfo(std::string);
  VarInfo(void);
  VarInfo(const VarInfo &);
  void init(void);
  
  float *read_cell_at_time(int,MeshInfo &);
  float *read_node_at_time(int,MeshInfo &);
  float *read_cell_z_at_time(int,MeshInfo &);
};


// ****************************************************************************
//  Class: avtUGRIDSingle
//
//  Purpose:
//      Reads in UGRID files for a single subdomain.
//
//  Programmer: rusty -- generated by xml2avt
//  Creation:   Sat Mar 5 18:27:45 PST 2016
//
// ****************************************************************************

class avtUGRIDSingle : public avtMTSDFileFormat
{
  public:
  avtUGRIDSingle(const char *);

  virtual ~avtUGRIDSingle();

  virtual int            GetNTimesteps(void);
  virtual void           GetTimes(std::vector<double> &);
  virtual void           GetCycles(std::vector<int> &);

  virtual const char    *GetType(void)   { return "UGRIDsingle"; };
  virtual void           FreeUpResources(void); 

  virtual vtkDataSet    *GetMesh(int, const char *);
  virtual vtkDataArray  *GetVar(int, const char *);
  virtual vtkDataArray  *GetVectorVar(int, const char *);

  vtkUnstructuredGrid *GetMeshNodes(const std::string);
  // vtkDataSet *ExtrudeTo3D(const std::string,int,vtkUnstructuredGrid *);
  vtkPoints *GetNodes(const std::string);

  void initialize_metadata(void);

  // DATA MEMBERS

  virtual void PopulateDatabaseMetaData(avtDatabaseMetaData *, int);

protected:

  int ncid; // handle for netcdf file
  int time_dim; // some variables/meshes may have a different time - that's not tested, tho.
  int time_var; 
  // int node_x_var,node_y_var;
  std::string default_ugrid_mesh;
  std::map<std::string,VarInfo> var_table;
  std::map<std::string,MeshInfo> mesh_table;

  // basic dimensions 

  // map 3D cell ids to real cells, because cells that are
  // underground are not output.
  // probably this needs to move to mesh
  std::map<int,int> full_cell2valid;

  // netcdf helpers:
  // read a full field of values at a given timestate
  // layer is the fastest changing index
  float *read_cell_z_full(std::string,int);

  std::string var_mesh2d(std::string);
  std::string var_mesh2d(int);
  bool setMeshInfo(VarInfo &);

  void activateTimestate(int);
  // various geometry-related data specific to the active timestate
  int active_timestate;
  int *cell_kmin;
  int *cell_kmax;
  int ncells_3d;

  int vertical_coordinate_for_dimension(int dim); 
  std::string create_3d_mesh(std::string mesh2d,int z_dim,int z_var);

  vtkDataArray *GetVar3D(int, VarInfo &);
  vtkDataArray *GetVar2D(int, VarInfo &);
};

// ****************************************************************************
//  Class: avtUGRIDFileFormat
//
//  Purpose:
//      Wraps one more single-domains into a multi-domain database
//
//  Programmer: rusty -- generated by xml2avt
//  Creation:   Sat Mar 5 18:27:45 PST 2016
//
// ****************************************************************************

class avtUGRIDFileFormat : public avtMTMDFileFormat
{
  public:
  avtUGRIDFileFormat(const char *);

  // virtual ~avtUGRIDFileFormat();

  virtual int            GetNTimesteps(void);
  virtual void           GetTimes(std::vector<double> &);
  virtual void           GetCycles(std::vector<int> &);

  virtual const char    *GetType(void)   { return "UGRID"; };
  virtual void           FreeUpResources(void); 

  virtual vtkDataSet    *GetMesh(int, int, const char *);
  virtual vtkDataArray  *GetVar(int, int, const char *);
  virtual vtkDataArray  *GetVectorVar(int, int, const char *);

protected:
  // DATA MEMBERS
  virtual void PopulateDatabaseMetaData(avtDatabaseMetaData *, int);

  std::vector<avtUGRIDSingle*> domain_cache;
  std::vector<avtDatabaseMetaData*> metadata_cache;
  std::vector<std::string> filenames;
  int domain_count;

  avtUGRIDSingle *subdomain(int domain);
  void populate_filenames(const char *);

  // I think this can be handled within the subdomains.
  // void activateTimestate(int);
  // various geometry-related data specific to the active timestate
  // int active_timestate;
};



#endif
