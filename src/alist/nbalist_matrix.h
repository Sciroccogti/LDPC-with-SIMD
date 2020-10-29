/*
 * @Author: Yifan Zhang
 * @Date: 2020-10-29 15:04:15
 * @Last Modified by: Yifan Zhang
 * @Last Modified time: 2020-10-29 15:31:36
 */

#ifndef NBALIST_H
#define NBALIST_H

// Non-Binary Alist
//
// N: number of rows
// M: number of cols
//
// num_n: weight of a row
// num_m: weight of a col
//
// m/nlist contains pairs of position and values for each non-zero element

#include <fstream>

typedef struct {
    int N, M;       /* size of the matrix */
    int GF;         /* GF of the matrix */
    int **mlist;    /* list of integer coordinates in the m direction where the
                       non-zero entries are and what values they take */
    int **nlist;    /* list of integer coordinates in the n direction where the
                       non-zero entries are and what values they take */
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
} nbalist_matrix;

// void write_nbalist(FILE *, nbalist_matrix *);
// int read_nbalist(const FILE *, nbalist_matrix *);
#endif