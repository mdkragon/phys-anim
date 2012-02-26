// Written by Peter Kutz.
// See Bridson's Fluids Notes, section 4.3.1 Putting Them In Matrix-Vector Form (on page 31, or PDF page 43)
// for an explanation of this symmetric sparse matrix data structure.

#ifndef GRID_DATA_MATRIX_H
#define GRID_DATA_MATRIX_H

#include "grid_data.h"

class GridDataMatrix {
private:
protected:
public:
	GridDataMatrix() {
		diag.initialize();
		plusI.initialize();
		plusJ.initialize();
		plusK.initialize();
	}
	GridData diag;
	GridData plusI;
	GridData plusJ;
	GridData plusK;

  void print() {
    vec3 dim = diag.getDim(); 
    //printf("%f %f %f\n", dim[0], dim[1], dim[2]);
    for (int k = 0; k < dim[2]; k++) {
      for (int j = 0; j < dim[1]; j++) {
        for (int i = 0; i < dim[0]; i++) {
          for (int kk = 0; kk < dim[2]; kk++) {
            for (int jj = 0; jj < dim[1]; jj++) {
              for (int ii = 0; ii < dim[0]; ii++) {
                //printf("(%d,%d,%d)", ii,jj,kk);
                if (i == ii && j == jj && k == kk) {
                  //printf("d% 8.2f  ", diag(i,j,k));
                  printf("% 8.2f  ", diag(i,j,k));
                } else {
                  if (i-1 == ii && j == jj && k == kk) {
                    printf("% 8.2f  ", plusI(ii,jj,kk));
                  } else if (i+1 == ii && j == jj && k == kk) {
                    printf("% 8.2f  ", plusI(i,jj,kk));
                  } else if (i == ii && j-1 == jj && k == kk) {
                    printf("% 8.2f  ", plusJ(ii,jj,kk));
                  } else if (i == ii && j+1 == jj && k == kk) {
                    printf("% 8.2f  ", plusJ(ii,j,kk));
                  } else if (i == ii && j == jj && k-1 == kk) {
                    printf("% 8.2f  ", plusK(ii,jj,kk));
                  } else if (i == ii && j == jj && k+1 == kk) {
                    printf("% 8.2f  ", plusK(ii,jj,k));
                  } else {
                    printf("% 8.2f  ", 0.0);
                  }
                    
                    

                  /*
                  // TODO: fix this logic so the neighbors are printed correctly
                  if (ii+1 <= dim[0] && plusI(ii,jj,kk)) {
                    //printf("i% 8.2f  ", plusI(ii,jj,kk));
                    printf("% 8.2f  ", plusI(ii,jj,kk));
                  } else if (jj+1 <= dim[1] && plusJ(ii,jj,kk)) {
                    //printf("j% 8.2f  ", plusJ(ii,jj,kk));
                    printf("% 8.2f  ", plusJ(ii,jj,kk));
                  } else if (kk+1 <= dim[2] && plusK(ii,jj,kk)) {
                    //printf("k% 8.2f  ", plusK(ii,jj,kk));
                    printf("% 8.2f  ", plusK(ii,jj,kk));
                  } else {
                    //printf("u% 8.2f  ", 0.0);
                    printf("% 8.2f  ", 0.0);
                  }
                  */
                }
              }
            }
          }
          printf("\n");
        }
      }
    }
  }
};

#endif // GRID_DATA_MATRIX_H
