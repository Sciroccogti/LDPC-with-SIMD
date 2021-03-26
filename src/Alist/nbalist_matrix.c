#include "Alist/nbalist_matrix.h"

#include "Alist/alist_matrix.h"

/**
 * @brief
 *
 * @param fp
 * @param list
 * @param GFlist the list of GF values for each none-zero element
 * @param length length of num_list
 * @param GF the GF order
 * @param num_list
 * @return int return 0 if ok
 */
int read_nbimatrix(FILE *fp, int **list, int **GFlist, const int length,
                   const int GF, const int *num_list) {
    int i, j;
    for (i = 0; i < length; i++) {
        list[i] = (int *)malloc(num_list[i] * sizeof(int));
        GFlist[i] = (int *)malloc(num_list[i] * sizeof(int));
        for (j = 0; j < num_list[i]; j++) {
            if (fscanf(fp, "%d %d", &list[i][j], &GFlist[i][j]) == EOF ||
                GFlist[i][j] >= GF) {
                return EOF;
            }
            // some alist contains 0
            if (!list[i][j]) {
                j--;
            }
        }
    }
    return 0;
}

// read alist from file, return 0 if ok
int read_nbalist(FILE *fp, nbalist_matrix *a) {
    /* this assumes that mlist and nlist have the form of a rectangular
       matrix in the file; if lists have unequal lengths, then the
       entries should be present (eg zero values) but are ignored
       */

    fscanf(fp, "%d %d %d\n", &a->N, &a->M, &a->GF);
    fscanf(fp, "%d %d\n", &a->biggest_num_n, &a->biggest_num_m);
    a->biggest_num_n_alloc = a->N * a->biggest_num_n;
    a->biggest_num_m_alloc = a->M * a->biggest_num_m;

    a->num_nlist = (int *)malloc(a->N * sizeof(int));
    a->num_mlist = (int *)malloc(a->M * sizeof(int));

    a->nlist = (int **)malloc(a->N * sizeof(int *));
    a->mlist = (int **)malloc(a->M * sizeof(int *));
    a->nGFlist = (int **)malloc(a->N * sizeof(int *));
    a->mGFlist = (int **)malloc(a->M * sizeof(int *));

    if (read_ivector(fp, a->num_nlist, 1, a->N) == EOF) {
        return EOF;
    }
    if (read_ivector(fp, a->num_mlist, 1, a->M) == EOF) {
        return EOF;
    }
    if (read_nbimatrix(fp, a->nlist, a->nGFlist, a->N, a->GF, a->num_nlist) ==
        EOF) {
        return EOF;
    }
    if (read_nbimatrix(fp, a->mlist, a->mGFlist, a->M, a->GF, a->num_mlist) ==
        EOF) {
        return EOF;
    }

    return 0;
}

void write_nbimatrix(FILE *fp, int *const *list, int *const *GFlist, const int length,
                     const int GF, const int *num_list) {
    int i, j;
    for (i = 0; i < length; i++) {
        for (j = 0; j < num_list[i]; j++) {
            fprintf(fp, "%d %d\t", list[i][j], GFlist[i][j]);
        }
        fprintf(fp, "\n");
    }
}

// write alist to file
void write_nbalist(FILE *fp, nbalist_matrix *a) {
    /* this assumes that mlist and nlist have the form of a rectangular
       matrix in the file; if lists have unequal lengths, then the
       entries should be present (eg zero values) but are ignored
       */
    int N = a->N, M = a->M, GF = a->GF;

    fprintf(fp, "%d %d %d\n", N, M, GF);
    fprintf(fp, "%d %d\n", a->biggest_num_n, a->biggest_num_m);
    write_ivector(fp, a->num_nlist, 1, N);
    write_ivector(fp, a->num_mlist, 1, M);
    write_nbimatrix(fp, a->nlist, a->nGFlist, N, GF, a->num_nlist);
    write_nbimatrix(fp, a->mlist, a->mGFlist, M, GF, a->num_mlist);
}