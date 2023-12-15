#!/usr/bin/env python3

import numpy as np
import os

def write_matrix(name, rowptr, colidx, values):
    if not os.path.exists(name):
        os.mkdir(f'{name}')
    np.array(rowptr, dtype=np.int32).tofile(f"{name}/rowptr.raw")
    np.array(colidx, dtype=np.int32).tofile(f"{name}/colidx.raw")
    np.array(values, dtype=np.float64).tofile(f"{name}/values.raw")

def write_vector(name, values):
    np.array(values, dtype=np.float64).tofile(f"{name}.raw")

def print_matrix(rowptr, colidx, values):
    print(len(rowptr) - 1)
    print(len(colidx))

    for i in range(0, len(rowptr)-1):
        b = rowptr[i]
        e = rowptr[i+1]

        print(f'{i}) ')
        for k in range(b, e):
            c = colidx[k]
            v = values[k]
            print(f'({c}, {v}) ')


def same_grid():
    rowptr = [ 0, 1, 4, 7, 10, 13, 14 ]
    colidx = [0, 
    		  0, 1, 2, 
    		  	 1, 2, 3,
    		  	 	2, 3, 4,
    		  	 	   3, 4, 5,
    		  	 	   		 5 ]
    values = [1, 
    		  -1, 2, -1, 
    		  	  -1, 2, -1, 
    		  	  	  -1, 2, -1,
    		  	  	  	  -1, 2, -1,
    		  	  	  	  		  1 ]

    print_matrix(rowptr, colidx, values)

    write_matrix("A", rowptr, colidx, values)

    rowptr = [ 0, 1, 2, 3, 4, 5, 6 ]
    colidx = [ 0, 1, 2, 3, 4, 5 ]
    values = [ 1, 1, 1, 1, 1, 1 ]

    write_matrix("T", rowptr, colidx, values)
    write_matrix("O", rowptr, colidx, values)

    is_contact = [ 0, 1, 1, 1, 1, 0 ]
    write_vector("is_contact", is_contact)

    v = 10
    rhs = [ 0, v, v, v, v, 0 ]
    write_vector("rhs", rhs)

    upper_bound = [ 0.3, 0.3, 0.3, 0.3, 0.3, 0.3]
    write_vector("upper_bound", upper_bound)

def smaller_grid():
    rowptr = [ 0, 1, 4, 7, 10, 13, 14 ]
    colidx = [ 0, 
               0, 1, 2, 
                  1, 2, 3,
                     2, 3, 4,
                        3, 4, 5,
                              5 ]
    values = [ 1, 
              -1, 2, -1, 
                  -1, 2, -1, 
                      -1, 2, -1,
                          -1, 2, -1,
                                  1 ]

    print_matrix(rowptr, colidx, values)

    write_matrix("A", rowptr, colidx, values)

    v = 10
    rhs = [ 0, v, v, v, v, 0 ]
    write_vector("rhs", rhs)

    # Contact

    rowptr = [ 0, 1, 2, 4 ]
    colidx = [ 2, 3, 4, 5 ]
    values = [ 1, 1, 1, 0 ]

    write_matrix("T", rowptr, colidx, values)

    rowptr = [ 0, 1, 2, 3]
    colidx = [ 0, 1, 2 ]
    values = [ 1, 1, 1 ]

    write_matrix("O", rowptr, colidx, values)

    is_contact = [ 1, 1, 1 ]
    write_vector("is_contact", is_contact)

    upper_bound = [ 0.3, 0.3, 0.3 ]
    write_vector("upper_bound", upper_bound)

def noncoforming_grid():
    rowptr = [ 0, 1, 4, 7, 10, 13, 14 ]
    colidx = [ 0, 
               0, 1, 2, 
                  1, 2, 3,
                     2, 3, 4,
                        3, 4, 5,
                             5 ]
    values = [ 1, 
              -1, 2, -1, 
                  -1, 2, -1, 
                      -1, 2, -1,
                          -1, 2, -1,
                                  1 ]

    print_matrix(rowptr, colidx, values)

    write_matrix("A", rowptr, colidx, values)

    v = 10
    rhs = [ 0, v, v, v, v, 0 ]
    write_vector("rhs", rhs)

    # Contact
    rowptr = [ 0, 2, 4, 7 ]
    colidx = [ 1, 2, 2, 3, 3, 4, 5 ]
    values = [ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0 ]

    write_matrix("T", rowptr, colidx, values)

    rowptr = [ 0, 1, 2, 3]
    colidx = [ 0, 1, 2 ]
    values = [ 1, 1, 1 ]

    write_matrix("O", rowptr, colidx, values)

    is_contact = [ 1, 1, 1 ]
    write_vector("is_contact", is_contact)

    upper_bound = [ 0.3, 0.3, 0.3 ]
    write_vector("upper_bound", upper_bound)


# same_grid()
# smaller_grid()
noncoforming_grid()
