#include <output.h>
#include <iostream>
#include <fstream>
// #include <stdio.h>

void output::output_data(char* name, double *v, double *r, int size, double time){
  #ifdef USE_SDF
  gft_out_full(name,time,&size,"r",1,r,v);
  #else
  static bool first = true;
  if(first){
    printf("[0;34mSDF output currently disabled. Further warnings will be suppressed.\n[0m");
    first = false;
  }
  #endif
}

void writeVTKScalarHelper(std::ofstream &file, double *f, int *shp, char *fname)
{
  int iter = 10;
  file << "        <DataArray Name=\"" << fname << "\" TimeStep=\"" << iter << "\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">" << std::endl;

  unsigned int c = 0;
  unsigned int nd = shp[0] * shp[1] * shp[2];
  while (c < nd) {
    unsigned int ncol = 0;
    while (ncol < 10 && c < nd) {
      file << f[c] << " ";
      ncol++;
      c++;
    }
    file << std::endl;
  }

  file << "        </DataArray>" << std::endl;
}

void writeVTKCoordsHelper(std::ofstream &file, const double *x, int nx, char *cname)
{
  //--- write coordinates ---
  file << "        <DataArray Name=\"" << cname << "\" type=\"Float64\" format=\"ascii\">" << std::endl;
  int c = 0;
  while (c < nx) {
    int cc = 0;
    while (cc < 10 && c < nx) {
      file << x[c] << " ";
      cc++;
      c++;
    }
    file << std::endl;
  }

  file << "        </DataArray>" << std::endl;
  //--- end write coordinates ---
}

void output::output_vtk(char *name, double *f, const double *x, int nx, double time, int iter)
{
  char filename[128];
  sprintf(filename, "%s.%04d.vtr",name,iter);
  
  std::ofstream file;
  file.open(filename);

  // define default coordinates for y and z
  int ny = 1;
  int nz = 1;
  double y[1] = {0.0};
  double z[1] = {0.0};
  int shp[3] = {nx, ny, nz};

  //--- write preamble to the file ---

  file << "<?xml version=\"1.0\"?>" << std::endl;

  file << "<VTKFile type=\"RectilinearGrid\" version=\"1.0\" byte_order=\"LittleEndian\">" << std::endl;
  file << "  <RectilinearGrid WholeExtent=\"";
  int bbox[6] = {0, nx-1, 0, ny-1, 0, nz-1};
  for (int i = 0; i < 5; i++) {
    file << bbox[i] << " ";
  }
  file << bbox[5] << "\">" << std::endl;

  file << "    <Piece Extent=\"";
  for (int i = 0; i < 5; i++) {
    file << bbox[i] << " ";
  }
  file << bbox[5] << "\">" << std::endl;

  //--- write time ---
  file << "      <FieldData>" << std::endl;
  file << "        <DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\">" << std::endl;
  file << "          " << time << std::endl;
  file << "        </DataArray>" << std::endl;
  file << "      </FieldData>" << std::endl;
  //--- end write time ---


  //--- write data ---
  file << "      <PointData>" << std::endl;
  writeVTKScalarHelper(file, f, shp, name);
  file << "      </PointData>" << std::endl;

  //--- write coordinates ---
  file << "      <Coordinates>" << std::endl;
  char cname[32]; 
  sprintf(cname, "x");
  writeVTKCoordsHelper(file, x, nx, cname);
  sprintf(cname, "y");
  writeVTKCoordsHelper(file, y, ny, cname);
  sprintf(cname, "z");
  writeVTKCoordsHelper(file, z, nz, cname);

  file << "      </Coordinates>" << std::endl;
  //--- end write coordinates ---

  file << "    </Piece>" << std::endl;
  file << "  </RectilinearGrid>" << std::endl;
  file << "</VTKFile>" << std::endl;

  file.close();
}


