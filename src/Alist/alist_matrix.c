#include "Alist/alist_matrix.h"

// return 0 if ok
int read_ivector(FILE *fp, int *num_list, int unk, const int length) {
    int i;
    for (i = 0; i < length; i++) {
        if (fscanf(fp, "%d ", &num_list[i]) == EOF) {
            return EOF;
        }
    }
    return 0;
}

// return 0 if ok
int read_imatrix(FILE *fp, int **list, int unk1, const int length, int unk2,
                 const int *num_list) {
    int i, j;
    for (i = 0; i < length; i++) {
        list[i] = (int *)malloc(num_list[i] * sizeof(int));
        for (j = 0; j < num_list[i]; j++) {
            if (fscanf(fp, "%d ", &list[i][j]) == EOF) {
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
int read_alist(FILE *fp, alist_matrix *a) {
    /* this assumes that mlist and nlist have the form of a rectangular
       matrix in the file; if lists have unequal lengths, then the
       entries should be present (eg zero values) but are ignored
       */

    fscanf(fp, "%d %d\n", &a->N, &a->M);
    fscanf(fp, "%d %d\n", &a->biggest_num_n, &a->biggest_num_m);
    a->biggest_num_n_alloc = a->N * a->biggest_num_n;
    a->biggest_num_m_alloc = a->M * a->biggest_num_m;

    a->num_nlist = (int *)malloc(a->N * sizeof(int));
    a->num_mlist = (int *)malloc(a->M * sizeof(int));

    a->nlist = (int **)malloc(a->N * sizeof(int *));
    a->mlist = (int **)malloc(a->M * sizeof(int *));

    if (read_ivector(fp, a->num_nlist, 1, a->N) == EOF) {
        return EOF;
    }
    if (read_ivector(fp, a->num_mlist, 1, a->M) == EOF) {
        return EOF;
    }
    if (read_imatrix(fp, a->nlist, 1, a->N, 1, a->num_nlist) == EOF) {
        return EOF;
    }
    if (read_imatrix(fp, a->mlist, 1, a->M, 1, a->num_mlist) == EOF) {
        return EOF;
    }

    return 0;
}

void write_ivector(FILE *fp, const int *num_list, int unk, const int length) {
    int i;
    for (i = 0; i < length; i++) {
        fprintf(fp, "%d ", num_list[i]);
    }
    fprintf(fp, "\n");
}

void write_imatrix(FILE *fp, int *const *list, int unk1, const int length,
                   int unk2, const int *num_list) {
    int i, j;
    for (i = 0; i < length; i++) {
        for (j = 0; j < num_list[i]; j++) {
            fprintf(fp, "%d ", list[i][j]);
        }
        fprintf(fp, "\n");
    }
}

// write alist to file
void write_alist(FILE *fp, alist_matrix *a) {
    /* this assumes that mlist and nlist have the form of a rectangular
       matrix in the file; if lists have unequal lengths, then the
       entries should be present (eg zero values) but are ignored
       */
    int N = a->N, M = a->M;

    fprintf(fp, "%d %d\n", N, M);
    fprintf(fp, "%d %d\n", a->biggest_num_n, a->biggest_num_m);
    write_ivector(fp, a->num_nlist, 1, N);
    write_ivector(fp, a->num_mlist, 1, M);
    write_imatrix(fp, a->nlist, 1, N, 1, a->num_nlist);
    write_imatrix(fp, a->mlist, 1, M, 1, a->num_mlist);
}