/**
 * @file alist_matrix.h
 * @author Yifan Zhang (scirocco_gti@yeah.net)
 * @brief
 * http://www.inference.org.uk/mackay/codes/alist.html
 * N: number of rows
 * M: number of cols
 * num_n: weight of a row
 * num_m: weight of a col
 * @date 2020-10-29 13:00:39
 * @modified: 2021-03-26 20:35:42
 */

#ifndef ALIST_MATRIX_H
#define ALIST_MATRIX_H

#include <malloc.h>

// a struct to store spare matrix
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

// alist_matrix support both C and CXX
#ifdef __cplusplus
extern "C" void write_alist(FILE *, alist_matrix *);
extern "C" int read_alist(FILE *, alist_matrix *);
#else
void write_alist(FILE *, alist_matrix *);
int read_alist(FILE *, alist_matrix *);
#endif

int read_ivector(FILE *fp, int *num_list, int unk, const int length);
void write_ivector(FILE *fp, const int *num_list, int unk, const int length);

#endif
