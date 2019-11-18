#include <stdio.h>
#include <stdlib.h>

/* Ricardo-Francisco LUIS MARTINEZ*/



void saveOVF(const float* ptr,const int cells[3], const float size[3],const char* filepath){
  int volume = cells[0]*cells[1]*cells[2];
  FILE *archivo;
  archivo = fopen(filepath,"w");
  unsigned char signle4[]={0x38,0xB4,0x96,0x49};
  fprintf(archivo,"%s \n","# OOMMF OVF 2.0" );
  fprintf(archivo,"%s \n","# Segment count: 1" );
  fprintf(archivo,"%s \n","# Begin: Segment" );
  fprintf(archivo,"%s \n","# Begin: Header" );
  fprintf(archivo,"%s \n","# Title: m" );
  fprintf(archivo,"%s \n","# meshtype: rectangular" );
  fprintf(archivo,"%s \n","# meshunit: m" );
  fprintf(archivo,"%s \n","# xmin: 0" );
  fprintf(archivo,"%s \n","# ymin: 0" );
  fprintf(archivo,"%s \n","# zmin: 0" );
  fprintf(archivo,"%s %g\n","# xmax: ", cells[0]*size[0] );
  fprintf(archivo,"%s %g\n","# ymax: ", cells[1]*size[1] );
  fprintf(archivo,"%s %g\n","# zmax: ", cells[2]*size[2] );
  fprintf(archivo,"%s \n","# valuedim: 3" );
  fprintf(archivo,"%s \n","# valuelabels: m_x m_y m_z" );
  fprintf(archivo,"%s \n","# valueunits: 1 1 1" );
  fprintf(archivo,"%s \n","# Desc: " );
  fprintf(archivo,"%s %g\n","# xbase: " , 0.5*size[0] );
  fprintf(archivo,"%s %g\n","# ybase: " , 0.5*size[1] );
  fprintf(archivo,"%s %g\n","# zbase: " , 0.5*size[2] );
  fprintf(archivo,"%s %d\n","# xnodes: ", cells[0] );
  fprintf(archivo,"%s %d\n","# ynodes: ", cells[1] );
  fprintf(archivo,"%s %d\n","# znodes: ", cells[2] );
  fprintf(archivo,"%s %g\n","# xstepsize: ", size[0] );
  fprintf(archivo,"%s %g\n","# ystepsize: ", size[1] );
  fprintf(archivo,"%s %g\n","# zstepsize: ", size[2] );
  fprintf(archivo,"%s \n","# End: Header" );
  fprintf(archivo,"%s \n","# Begin: Data Binary 4" );
  fwrite(signle4,1,4,archivo);
  fwrite((char*) ptr,12,volume,archivo);
  fprintf(archivo,"%s \n","# End: Data Binary 4" );
  fprintf(archivo,"%s \n","# End: Segment" );
  fclose(archivo);
}
