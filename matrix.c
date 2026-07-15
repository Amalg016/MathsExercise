#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "matrix.h"

[[nodiscard]]
Matrix *matrix_new(size_t rows, size_t cols)
{
    Matrix *m = malloc(sizeof(*m));
    if (!m)
        return nullptr;

    m->rows = rows;
    m->cols = cols;
    m->data = calloc(rows * cols, sizeof(double));

    if (!m->data) {
        free(m);
        return nullptr;
    }

    return m;
}

void matrix_free(Matrix *m)
{
    if (!m)
        return;

    free(m->data);
    free(m);
}

void matrix_print(const Matrix *m)
{
    for (size_t i = 0; i < MAT_ROWS(m); ++i) {
        for (size_t j = 0; j < MAT_COLS(m); ++j)
            printf("%10.4f ", MAT_AT(m, i, j));

        putchar('\n');
    }
}

[[nodiscard]]
Matrix *matrix_add(const Matrix *a, const Matrix *b)
{
    if (MAT_ROWS(a) != MAT_ROWS(b) ||
        MAT_COLS(a) != MAT_COLS(b))
        return nullptr;

    Matrix *c = matrix_new(MAT_ROWS(a), MAT_COLS(a));
    if (!c)
        return nullptr;

    for (size_t i = 0; i < MAT_SIZE(a); ++i)
        c->data[i] = a->data[i] + b->data[i];

    return c;
}

[[nodiscard]]
Matrix *matrix_sub(const Matrix *a, const Matrix *b)
{
    if (MAT_ROWS(a) != MAT_ROWS(b) ||
        MAT_COLS(a) != MAT_COLS(b))
        return nullptr;

    Matrix *c = matrix_new(MAT_ROWS(a), MAT_COLS(a));
    if (!c)
        return nullptr;

    for (size_t i = 0; i < MAT_SIZE(a); ++i)
        c->data[i] = a->data[i] - b->data[i];

    return c;
}

[[nodiscard]]
Matrix *matrix_mul(const Matrix *a, const Matrix *b)
{
    if (MAT_COLS(a) != MAT_ROWS(b))
        return nullptr;

    Matrix *c = matrix_new(MAT_ROWS(a), MAT_COLS(b));
    if (!c)
        return nullptr;

    for (size_t i = 0; i < MAT_ROWS(a); ++i)
        for (size_t k = 0; k < MAT_COLS(a); ++k) {
            const double aik = MAT_AT(a, i, k);

            for (size_t j = 0; j < MAT_COLS(b); ++j)
                MAT_AT(c, i, j) +=
                    aik * MAT_AT(b, k, j);
        }

    return c;
}

[[nodiscard]]
Matrix *matrix_transpose(const Matrix *m)
{
    Matrix *t = matrix_new(MAT_COLS(m), MAT_ROWS(m));
    if (!t)
        return nullptr;

    for (size_t i = 0; i < MAT_ROWS(m); ++i)
        for (size_t j = 0; j < MAT_COLS(m); ++j)
            MAT_AT(t, j, i) = MAT_AT(m, i, j);

    return t;
}

[[nodiscard]]
double matrix_det(const Matrix *m)
{
//Should always be a square matrix. As determinant exists only for a square matrix
    if (MAT_ROWS(m) != MAT_COLS(m))
        return NAN;

    if (MAT_ROWS(m) == 1)
        return MAT_AT(m, 0, 0);

    if (MAT_ROWS(m) == 2)
        return MAT_AT(m, 0, 0) * MAT_AT(m, 1, 1) -
               MAT_AT(m, 0, 1) * MAT_AT(m, 1, 0);

    const size_t n = MAT_ROWS(m);
    double det = 0.0;

// Using cofactor expansion for size greater than 2
    for (size_t col = 0; col < n; ++col) {

        Matrix *minor = matrix_new(n - 1, n - 1);
        if (!minor)
            return NAN;

        size_t r2 = 0;

        for (size_t r = 1; r < n; ++r) {

            size_t c2 = 0;

            for (size_t c = 0; c < n; ++c) {

                if (c == col)
                    continue;

                MAT_AT(minor, r2, c2++) = MAT_AT(m, r, c);
            }

            ++r2;
        }

        double sign = (col & 1u) ? -1.0 : 1.0;

        det += sign * MAT_AT(m, 0, col) * matrix_det(minor);

        matrix_free(minor);
    }

    return det;
}

[[nodiscard]]
Matrix *matrix_inverse(const Matrix *m)
{
    if (MAT_ROWS(m) != MAT_COLS(m))
        return nullptr;

    const size_t n = MAT_ROWS(m);

    Matrix *aug = matrix_new(n, 2 * n);
    if (!aug)
        return nullptr;

    for (size_t i = 0; i < n; ++i) {

        for (size_t j = 0; j < n; ++j)
            MAT_AT(aug, i, j) = MAT_AT(m, i, j);

        MAT_AT(aug, i, n + i) = 1.0;
    }

    for (size_t i = 0; i < n; ++i) {

        double pivot = MAT_AT(aug, i, i);

        if (fabs(pivot) < 1e-12) {
            matrix_free(aug);
            return nullptr;
        }

        for (size_t j = 0; j < 2 * n; ++j)
            MAT_AT(aug, i, j) /= pivot;

        for (size_t r = 0; r < n; ++r) {

            if (r == i)
                continue;

            double factor = MAT_AT(aug, r, i);

            for (size_t c = 0; c < 2 * n; ++c)
                MAT_AT(aug, r, c) -=
                    factor * MAT_AT(aug, i, c);
        }
    }

    Matrix *inv = matrix_new(n, n);
    if (!inv) {
        matrix_free(aug);
        return nullptr;
    }

    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            MAT_AT(inv, i, j) = MAT_AT(aug, i, n + j);

    matrix_free(aug);

    return inv;
}

bool matrix_is_invertible(const Matrix *m)
{
    return fabs(matrix_det(m)) > 1e-12;
}
