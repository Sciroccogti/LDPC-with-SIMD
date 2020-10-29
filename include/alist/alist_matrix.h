/*
 * @Author: Yifan Zhang
 * @Date: 2020-10-29 15:02:21
 * @Last Modified by: Yifan Zhang
 * @Last Modified time: 2020-10-29 18:19:27
 */

#ifndef ALIST_MATRIX_H
#define ALIST_MATRIX_H

// http://www.inference.org.uk/mackay/codes/alist.html
// N: number of rows
// M: number of cols
// num_n: weight of a row
// num_m: weight of a col

#include <malloc.h>

typedef struct {
    int N, M;       /* size of the matrix */
    int **mlist;    /* list of integer coordinates in the m direction where the
                       non-zero entries are */
    int **nlist;    /* list of integer coordinates in the n direction where the
                       non-zero entries are */
    int *num_mlist; /* weight of each row, m */
    int *num_nlist; /* weight of each column n */
    int *l_up_to;   // todo
    int *u_up_to;   // todo
    int *norder;    // todo
    int biggest_num_m; /* actual biggest sizes */
    int biggest_num_n;
    int biggest_num_m_alloc; /* sizes used for memory allocation */
    int biggest_num_n_alloc;
    int tot;         // todo
    int same_length; /* whether all vectors in mlist and nlist have same length
                      */
} alist_matrix;

void write_alist(FILE *, alist_matrix *);
int read_alist(FILE *, alist_matrix *);

// void init_alist(alist_matrix *);
#endif
