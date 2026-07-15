#ifndef MATRIX_H
#define MATRIX_H

#include <stddef.h>
#include <stdbool.h>

typedef struct {
    size_t rows;
    size_t cols;
    double *data;
} Matrix;

[[nodiscard]]
Matrix *matrix_new(size_t rows, size_t cols);

void matrix_free(Matrix *m);

[[nodiscard]]
double matrix_get(const Matrix *m, size_t row, size_t col);

bool matrix_set(Matrix *m,
                size_t row,
                size_t col,
                double value);

[[nodiscard]]
Matrix *matrix_add(const Matrix *a,
                   const Matrix *b);

[[nodiscard]]
Matrix *matrix_sub(const Matrix *a,
                   const Matrix *b);

[[nodiscard]]
Matrix *matrix_mul(const Matrix *a,
                   const Matrix *b);

[[nodiscard]]
Matrix *matrix_transpose(const Matrix *m);

[[nodiscard]]
double matrix_det(const Matrix *m);

[[nodiscard]]
Matrix *matrix_inverse(const Matrix *m);

void matrix_print(const Matrix *m);
bool matrix_is_invertible(const Matrix *m);

#define MAT_INDEX(m, r, c) ((size_t)(r) * (m)->cols + (size_t)(c))
#define MAT_AT(m, r, c)    ((m)->data[MAT_INDEX((m), (r), (c))])
#define MAT_ROWS(m)        ((m)->rows)
#define MAT_COLS(m)        ((m)->cols)
#define MAT_SIZE(m)        ((m)->rows * (m)->cols)
#endif
