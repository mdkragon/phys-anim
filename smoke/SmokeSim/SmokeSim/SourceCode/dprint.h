#ifndef __DPRINT_H__
#define __DPRINT_H__

// enable debug printing?
#define __DPRINT__

// print 3d matrix
/*
 * Prints a 3D matrix to screen if the display grid order
 *  X-left, Y-up, Z-into screen
 *  each row contains | X1 | X2 | X3 | ... | Xn  (where n is the z dimension)
 *  
 * ex: 2x2x2 grid (with values 1:8)
 *  | 4 3 | 8 7 |
 *  | 2 1 | 6 5 |
 *
 * the coordinates:
 *  | (1,1,0) (0,1,0) | (1,1,1) (0,1,1) |
 *  | (1,0,0) (0,0,0) | (1,0,1) (0,0,1) |
 */
// 
void print_grid_data(GridData &grid) {
  // matrix dimensions
  vec3 dim = grid.getDim();
  int xdim = dim[0];
  int ydim = dim[1];
  int zdim = dim[2];

  // loop over rows (starting from the top)
  for (int r = ydim-1; r >= 0; r--) {
    // loop over stacks
    printf("| ");
    for (int s = 0; s < zdim; s++) {
      // loop over columns (starting from end)
      //for (int c = 0; c < xdim; c++) {
      for (int c = xdim-1; c >= 0; c--) {
        printf("% 8.2f  ", grid(c, r, s));
      }
      printf("  |  ");
    }
    printf("\n");
  }

  fflush(stdout);
}

void print_grid_data_as_column(GridData &grid) {
  // matrix dimensions
  vec3 dim = grid.getDim();

  for (int k = 0; k < dim[2]; k++) {
    for (int j = 0; j < dim[1]; j++) {
      for (int i = 0; i < dim[0]; i++) {
        printf("% 8.2f\n", grid(i,j,k));
      }
    }
  }
  fflush(stdout);
}

void print_grid_data_as_row(GridData &grid) {
  // matrix dimensions
  vec3 dim = grid.getDim();

  for (int k = 0; k < dim[2]; k++) {
    for (int j = 0; j < dim[1]; j++) {
      for (int i = 0; i < dim[0]; i++) {
        printf("% 8.2f  ", grid(i,j,k));
      }
    }
  }
  fflush(stdout);
}

/*
 * prints the sparse grid
 */
void print_grid_data_matrix(GridDataMatrix &grid) {
  // matrix dimensions
}



#endif

