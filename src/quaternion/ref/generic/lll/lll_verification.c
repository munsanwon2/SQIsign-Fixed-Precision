#include "lll_internals.h"

// functions to verify lll
void
quat_lll_set_ibq_parameters(ibq_t *delta, ibq_t *eta)
{
    ibz_t num, denom;
    ibz_init(&num);
    ibz_init(&denom);
    ibq_set(delta, &ibz_const_one, &ibz_const_two);
    ibz_set(&num, EPSILON_NUM);
    ibz_set(&denom, EPSILON_DENOM);
    ibq_set(eta, &num, &denom);
    ibq_add(eta, eta, delta);
    ibz_set(&num, DELTA_NUM);
    ibz_set(&denom, DELTA_DENOM);
    ibq_set(delta, &num, &denom);
    ibz_finalize(&num);
    ibz_finalize(&denom);
}

void
ibq_vec_4_copy_ibz(ibq_vec_4_t *vec, const ibz_t *coeff0, const ibz_t *coeff1, const ibz_t *coeff2, const ibz_t *coeff3)
{
    ibz_t one;
    ibz_init(&one);
    ibz_set(&one, 1);
    ibq_set(&((*vec)[0]), coeff0, &one);
    ibq_set(&((*vec)[1]), coeff1, &one);
    ibq_set(&((*vec)[2]), coeff2, &one);
    ibq_set(&((*vec)[3]), coeff3, &one);
    ibz_finalize(&one);
}

void
quat_lll_bilinear(ibq_t *b, const ibq_vec_4_t *vec0, const ibq_vec_4_t *vec1, const ibz_t *q)
{
    ibq_t sum, prod, norm_q;
    ibz_t one;
    ibz_init(&one);
    ibz_set(&one, 1);
    ibq_init(&sum);
    ibq_init(&prod);
    ibq_init(&norm_q);
    ibq_set(&norm_q, q, &one);

    ibq_mul(&sum, &((*vec0)[0]), &((*vec1)[0]));
    ibq_mul(&prod, &((*vec0)[1]), &((*vec1)[1]));
    ibq_add(&sum, &sum, &prod);
    ibq_mul(&prod, &((*vec0)[2]), &((*vec1)[2]));
    ibq_mul(&prod, &prod, &norm_q);
    ibq_add(&sum, &sum, &prod);
    ibq_mul(&prod, &((*vec0)[3]), &((*vec1)[3]));
    ibq_mul(&prod, &prod, &norm_q);
    ibq_add(b, &sum, &prod);

    ibz_finalize(&one);
    ibq_finalize(&sum);
    ibq_finalize(&prod);
    ibq_finalize(&norm_q);
}

void
quat_lll_gram_schmidt_transposed_with_ibq(ibq_mat_4x4_t *orthogonalised_transposed,
                                          const ibz_mat_4x4_t *mat,
                                          const ibz_t *q)
{
    ibq_mat_4x4_t work;
    ibq_vec_4_t vec;
    ibq_t norm, b, coeff, prod;
    ibq_init(&norm);
    ibq_init(&coeff);
    ibq_init(&prod);
    ibq_init(&b);
    ibq_mat_4x4_init(&work);
    ibq_vec_4_init(&vec);
    // transpose the input matrix to be able to work on vectors
    for (int i = 0; i < 4; i++) {
        ibq_vec_4_copy_ibz(&(work[i]), &((*mat)[0][i]), &((*mat)[1][i]), &((*mat)[2][i]), &((*mat)[3][i]));
    }

    for (int i = 0; i < 4; i++) {
        quat_lll_bilinear(&norm, &(work[i]), &(work[i]), q);
        ibq_inv(&norm, &norm);
        for (int j = i + 1; j < 4; j++) {
            ibq_vec_4_copy_ibz(&vec, &((*mat)[0][j]), &((*mat)[1][j]), &((*mat)[2][j]), &((*mat)[3][j]));
            quat_lll_bilinear(&b, &(work[i]), &vec, q);
            ibq_mul(&coeff, &norm, &b);
            for (int k = 0; k < 4; k++) {
                ibq_mul(&prod, &coeff, &(work[i][k]));
                ibq_sub(&(work[j][k]), &(work[j][k]), &prod);
            }
        }
    }

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ibq_copy(&((*orthogonalised_transposed)[i][j]), &(work[i][j]));
        }
    }

    ibq_finalize(&norm);
    ibq_finalize(&coeff);
    ibq_finalize(&prod);
    ibq_finalize(&b);
    ibq_mat_4x4_finalize(&work);
    ibq_vec_4_finalize(&vec);
}

// int
// quat_lll_verify(const ibz_mat_4x4_t *mat, const ibq_t *delta, const ibq_t *eta, const quat_alg_t *alg)
// {
//     int res = 1;
//     ibq_mat_4x4_t orthogonalised_transposed;
//     ibq_vec_4_t tmp_vec;
//     ibq_t div, tmp, mu, two, norm, b;
//     ibz_t mu2_floored, num, denom;
//     ibq_mat_4x4_init(&orthogonalised_transposed);
//     ibq_vec_4_init(&tmp_vec);
//     ibq_init(&div);
//     ibq_init(&tmp);
//     ibq_init(&norm);
//     ibq_init(&b);
//     ibq_init(&mu);
//     ibq_init(&two);
//     ibz_init(&mu2_floored);
//     ibz_init(&num);
//     ibz_init(&denom);
//     ibz_set(&num, 2);
//     ibz_set(&denom, 1);
//     ibq_set(&two, &num, &denom);

//     quat_lll_gram_schmidt_transposed_with_ibq(&orthogonalised_transposed, mat, &(alg->p));
//     // check small bilinear products/norms
//     for (int i = 0; i < 4; i++) {
//         for (int j = 0; j < i; j++) {
//             ibq_vec_4_copy_ibz(&tmp_vec, &((*mat)[0][i]), &((*mat)[1][i]), &((*mat)[2][i]), &((*mat)[3][i]));
//             quat_lll_bilinear(&b, &(orthogonalised_transposed[j]), &tmp_vec, &(alg->p));
//             quat_lll_bilinear(&norm, &(orthogonalised_transposed[j]), &(orthogonalised_transposed[j]), &(alg->p));
//             ibq_inv(&tmp, &norm);
//             ibq_mul(&mu, &b, &tmp);
//             ibq_abs(&mu, &mu);
//             // compare to eta
//             res = res && (ibq_cmp(&mu, eta) <= 0);
//         }
//     }
//     for (int i = 1; i < 4; i++) {
//         ibq_vec_4_copy_ibz(&tmp_vec, &((*mat)[0][i]), &((*mat)[1][i]), &((*mat)[2][i]), &((*mat)[3][i]));
//         quat_lll_bilinear(&b, &(orthogonalised_transposed[i - 1]), &tmp_vec, &(alg->p));
//         quat_lll_bilinear(&norm, &(orthogonalised_transposed[i - 1]), &(orthogonalised_transposed[i - 1]), &(alg->p));
//         ibq_inv(&tmp, &norm);
//         ibq_mul(&mu, &b, &tmp);
//         // tmp is mu^2
//         ibq_mul(&tmp, &mu, &mu);
//         // mu is delta-mu^2
//         ibq_sub(&mu, delta, &tmp);
//         quat_lll_bilinear(&tmp, &(orthogonalised_transposed[i]), &(orthogonalised_transposed[i]), &(alg->p));
//         // get (delta-mu^2)norm(i-1)
//         ibq_mul(&div, &norm, &mu);
//         res = res && (ibq_cmp(&tmp, &div) >= 0);
//     }
//     ibq_mat_4x4_finalize(&orthogonalised_transposed);
//     ibq_vec_4_finalize(&tmp_vec);
//     ibq_finalize(&div);
//     ibq_finalize(&norm);
//     ibq_finalize(&b);
//     ibq_finalize(&tmp);
//     ibq_finalize(&mu);
//     ibq_finalize(&two);
//     ibz_finalize(&mu2_floored);
//     ibz_finalize(&num);
//     ibz_finalize(&denom);
//     return (res);
// }

// int
// quat_lll_verify(const ibz_mat_4x4_t *mat,
//                       const ibq_t *delta,
//                       const ibq_t *eta,
//                       const quat_alg_t *alg)
// {
//     int res = 1;
//     ibq_mat_4x4_t orthogonalised_transposed;
//     ibq_vec_4_t tmp_vec;
//     ibq_t div, tmp, mu, two, norm, b;
//     ibz_t mu2_floored, num, denom;

//     ibq_mat_4x4_init(&orthogonalised_transposed);
//     ibq_vec_4_init(&tmp_vec);
//     ibq_init(&div);
//     ibq_init(&tmp);
//     ibq_init(&norm);
//     ibq_init(&b);
//     ibq_init(&mu);
//     ibq_init(&two);
//     ibz_init(&mu2_floored);
//     ibz_init(&num);
//     ibz_init(&denom);

//     ibz_set(&num, 2);
//     ibz_set(&denom, 1);
//     ibq_set(&two, &num, &denom);

//     // Gram–Schmidt (transposed) with ibq
//     quat_lll_gram_schmidt_transposed_with_ibq(&orthogonalised_transposed, mat, &(alg->p));

//     /********************************************************************
//      * 1. 사이즈 리덕션 조건: |mu_{i,j}| <= eta  인지 체크
//      ********************************************************************/
//     for (int i = 0; i < 4; i++) {
//         for (int j = 0; j < i; j++) {
//             ibq_vec_4_copy_ibz(&tmp_vec,
//                                &((*mat)[0][i]),
//                                &((*mat)[1][i]),
//                                &((*mat)[2][i]),
//                                &((*mat)[3][i]));

//             // b = <b_j^*, b_i>
//             quat_lll_bilinear(&b, &(orthogonalised_transposed[j]), &tmp_vec, &(alg->p));
//             // norm = <b_j^*, b_j^*>
//             quat_lll_bilinear(&norm, &(orthogonalised_transposed[j]), &(orthogonalised_transposed[j]), &(alg->p));

//             ibq_inv(&tmp, &norm);
//             ibq_mul(&mu, &b, &tmp);   // mu = b / norm
//             ibq_abs(&mu, &mu);        // |mu|

//             if (ibq_cmp(&mu, eta) > 0) {
//                 // 여기서 사이즈 리덕션 조건이 깨짐
//                 printf("[LLL verify] size-reduction FAIL at (i=%d, j=%d)\n", i, j);
//                 // 여기에 ibq/ibz 출력 함수가 있다면 같이 찍어준다.
//                 // 예: ibq_print가 있다면
//                 // ibq_print("  |mu|  = ", &mu);
//                 // ibq_print("  eta  = ", eta);

//                 res = 0;
//             }
//         }
//     }

//     /********************************************************************
//      * 2. Lovász 조건: 
//      *    ||b_i^*||^2 >= (delta - mu^2) * ||b_{i-1}^*||^2
//      ********************************************************************/
//     for (int i = 1; i < 4; i++) {
//         ibq_vec_4_copy_ibz(&tmp_vec,
//                            &((*mat)[0][i]),
//                            &((*mat)[1][i]),
//                            &((*mat)[2][i]),
//                            &((*mat)[3][i]));

//         // b = <b_{i-1}^*, b_i>
//         quat_lll_bilinear(&b, &(orthogonalised_transposed[i - 1]), &tmp_vec, &(alg->p));
//         // norm = <b_{i-1}^*, b_{i-1}^*>
//         quat_lll_bilinear(&norm, &(orthogonalised_transposed[i - 1]), &(orthogonalised_transposed[i - 1]), &(alg->p));

//         // mu = b / norm
//         ibq_inv(&tmp, &norm);
//         ibq_mul(&mu, &b, &tmp);

//         // tmp = mu^2
//         ibq_mul(&tmp, &mu, &mu);
//         // mu = delta - mu^2
//         ibq_sub(&mu, delta, &tmp);

//         // tmp = ||b_i^*||^2
//         quat_lll_bilinear(&tmp, &(orthogonalised_transposed[i]), &(orthogonalised_transposed[i]), &(alg->p));
//         // div = (delta - mu^2) * ||b_{i-1}^*||^2
//         ibq_mul(&div, &norm, &mu);

//         if (ibq_cmp(&tmp, &div) < 0) {
//             // 여기서 Lovász 조건이 깨짐
//             printf("[LLL verify] Lovasz FAIL at i=%d\n", i);
//             // 있으면 찍어보기
//             // ibq_print("  ||b_i^*||^2            = ", &tmp);
//             // ibq_print("  (delta - mu^2)||b_{i-1}^*||^2 = ", &div);
//             // ibq_print("  delta = ", delta);

//             res = 0;
//         }
//     }

//     ibq_mat_4x4_finalize(&orthogonalised_transposed);
//     ibq_vec_4_finalize(&tmp_vec);
//     ibq_finalize(&div);
//     ibq_finalize(&norm);
//     ibq_finalize(&b);
//     ibq_finalize(&tmp);
//     ibq_finalize(&mu);
//     ibq_finalize(&two);
//     ibz_finalize(&mu2_floored);
//     ibz_finalize(&num);
//     ibz_finalize(&denom);

//     return res;
// }

static void
ibz_mat_4x4_print_s(const ibz_mat_4x4_t *M)
{
    printf("ibz_mat_4x4:\n");
    for (int i = 0; i < 4; i++) {
        printf("  ");
        for (int j = 0; j < 4; j++) {
            ibz_print(&(*M)[i][j], 10);
            printf(" ");
        }
        printf("\n");
    }
    printf("\n");
}

static void
quat_lll_debug_gram_from_basis(const ibz_mat_4x4_t *mat,
                               const quat_alg_t *alg)
{
    ibq_mat_4x4_t Gq;
    ibq_vec_4_t vi, vj;
    ibq_t b;
    ibq_mat_4x4_init(&Gq);
    ibq_vec_4_init(&vi);
    ibq_vec_4_init(&vj);
    ibq_init(&b);

    printf("=== Basis (columns) ===\n");
    ibz_mat_4x4_print_s(mat);

    printf("=== Gram matrix from quat_lll_bilinear (as rationals) ===\n");
    for (int i = 0; i < 4; i++) {
        // b_i = column i
        ibq_vec_4_copy_ibz(&vi,
                           &(*mat)[0][i],
                           &(*mat)[1][i],
                           &(*mat)[2][i],
                           &(*mat)[3][i]);
        for (int j = 0; j < 4; j++) {
            ibq_vec_4_copy_ibz(&vj,
                               &(*mat)[0][j],
                               &(*mat)[1][j],
                               &(*mat)[2][j],
                               &(*mat)[3][j]);
            quat_lll_bilinear(&b, &vi, &vj, &(alg->p));
            printf("G[%d,%d] = ", i, j);
            ibz_print(&(b[0]), 10);
            printf("/");
            ibz_print(&(b[1]), 10);
            printf("\n");
        }
    }

    ibq_mat_4x4_finalize(&Gq);
    ibq_vec_4_finalize(&vi);
    ibq_vec_4_finalize(&vj);
    ibq_finalize(&b);
}


int
quat_lll_verify(const ibz_mat_4x4_t *mat,
                      const ibq_t *delta,
                      const ibq_t *eta,
                      const quat_alg_t *alg)
{
    int res = 1;
    ibq_mat_4x4_t orthogonalised_transposed;
    ibq_vec_4_t tmp_vec;
    ibq_t div, tmp, mu, norm, b;
    ibq_t mu_sq, lhs, rhs;

    ibq_mat_4x4_init(&orthogonalised_transposed);
    ibq_vec_4_init(&tmp_vec);
    ibq_init(&div);
    ibq_init(&tmp);
    ibq_init(&norm);
    ibq_init(&b);
    ibq_init(&mu);

    ibq_init(&mu_sq);
    ibq_init(&lhs);
    ibq_init(&rhs);

    /* Gram–Schmidt (transposed) */
    quat_lll_gram_schmidt_transposed_with_ibq(&orthogonalised_transposed,
                                              mat, &(alg->p));

    /* 필요하면 전체 직교화 행렬을 한 번 찍어볼 수도 있음 */
    // printf("=== GS-transposed ===\n");
    // ibq_mat_4x4_print(&orthogonalised_transposed);

    /********************************************************************
     * 1. size reduction 조건: |mu_{i,j}| <= eta
     ********************************************************************/
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < i; j++) {
            ibq_vec_4_copy_ibz(&tmp_vec,
                               &((*mat)[0][i]),
                               &((*mat)[1][i]),
                               &((*mat)[2][i]),
                               &((*mat)[3][i]));

            /* b = <b_j^*, b_i> */
            quat_lll_bilinear(&b, &(orthogonalised_transposed[j]),
                              &tmp_vec, &(alg->p));
            /* norm = <b_j^*, b_j^*> */
            quat_lll_bilinear(&norm, &(orthogonalised_transposed[j]),
                              &(orthogonalised_transposed[j]), &(alg->p));

            ibq_inv(&tmp, &norm);
            ibq_mul(&mu, &b, &tmp);   /* mu = b / norm */
            ibq_abs(&mu, &mu);        /* |mu| */

            if (ibq_cmp(&mu, eta) > 0) {
                printf("[LLL verify] size-reduction FAIL at (i=%d, j=%d)\n", i, j);
                ibq_print_scalar("    |mu| =", &mu);
                ibq_print_scalar("    eta =", eta);
                res = 0;
            }
        }
    }

    /********************************************************************
     * 2. Lovász 조건: ||b_i^*||^2 >= (delta - mu^2) * ||b_{i-1}^*||^2
     ********************************************************************/
    for (int i = 1; i < 4; i++) {
        ibq_vec_4_copy_ibz(&tmp_vec,
                           &((*mat)[0][i]),
                           &((*mat)[1][i]),
                           &((*mat)[2][i]),
                           &((*mat)[3][i]));

        /* b = <b_{i-1}^*, b_i> */
        quat_lll_bilinear(&b, &(orthogonalised_transposed[i - 1]),
                          &tmp_vec, &(alg->p));
        /* norm = <b_{i-1}^*, b_{i-1}^*> */
        quat_lll_bilinear(&norm, &(orthogonalised_transposed[i - 1]),
                          &(orthogonalised_transposed[i - 1]), &(alg->p));

        /* mu = b / norm */
        ibq_inv(&tmp, &norm);
        ibq_mul(&mu, &b, &tmp);

        /* mu_sq = mu^2 */
        ibq_mul(&mu_sq, &mu, &mu);

        /* tmp = delta - mu^2 */
        ibq_sub(&tmp, delta, &mu_sq);

        /* lhs = ||b_i^*||^2 */
        quat_lll_bilinear(&lhs, &(orthogonalised_transposed[i]),
                          &(orthogonalised_transposed[i]), &(alg->p));
        /* rhs = (delta - mu^2) * ||b_{i-1}^*||^2 */
        ibq_mul(&rhs, &norm, &tmp);

        if (ibq_cmp(&lhs, &rhs) < 0) {
            printf("[LLL verify] Lovasz FAIL at i=%d\n", i);

            ibq_print_scalar("    mu           =", &mu);
            ibq_print_scalar("    mu^2         =", &mu_sq);
            ibq_print_scalar("    delta        =", delta);
            ibq_print_scalar("    norm_prev    =", &norm);
            ibq_print_scalar("    lhs=||b_i^*||^2   =", &lhs);
            ibq_print_scalar("    rhs=(delta-mu^2)||b_{i-1}^*||^2 =", &rhs);

            /* difference: lhs - rhs */
            ibq_sub(&tmp, &lhs, &rhs);
            ibq_print_scalar("    lhs - rhs    =", &tmp);
            
            quat_lll_debug_gram_from_basis(mat, alg);

            res = 0;
        }
    }

    ibq_mat_4x4_finalize(&orthogonalised_transposed);
    ibq_vec_4_finalize(&tmp_vec);
    ibq_finalize(&div);
    ibq_finalize(&norm);
    ibq_finalize(&b);
    ibq_finalize(&tmp);
    ibq_finalize(&mu);
    ibq_finalize(&mu_sq);
    ibq_finalize(&lhs);
    ibq_finalize(&rhs);

    return res;
}
