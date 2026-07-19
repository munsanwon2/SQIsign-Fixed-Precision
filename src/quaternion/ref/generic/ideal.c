#include <quaternion.h>
#include <stdlib.h>
#include "internal.h"

// assumes parent order and lattice correctly set, computes and sets the norm
void
quat_lideal_norm(quat_left_ideal_t *lideal)
{
    quat_lattice_index(&(lideal->norm), &(lideal->lattice), (lideal->parent_order));
    if (!ibz_sqrt(&(lideal->norm), &(lideal->norm))) {
        // A left ideal of a maximal order has square relative index.  Treat a
        // violation as malformed arithmetic state in both debug and release.
        abort();
    }
}

// assumes parent order and lattice correctly set, recomputes and verifies its norm
static int
quat_lideal_norm_verify(const quat_left_ideal_t *lideal)
{
    int res;
    ibz_t index;
    ibz_init(&index);
    quat_lattice_index(&index, &(lideal->lattice), (lideal->parent_order));
    res = ibz_sqrt(&index, &index) && (ibz_cmp(&(lideal->norm), &index) == 0);
    ibz_finalize(&index);
    return (res);
}

void
quat_lideal_copy(quat_left_ideal_t *copy, const quat_left_ideal_t *copied)
{
    copy->parent_order = copied->parent_order;
    ibz_copy(&copy->norm, &copied->norm);
    ibz_copy(&copy->lattice.denom, &copied->lattice.denom);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ibz_copy(&copy->lattice.basis[i][j], &copied->lattice.basis[i][j]);
        }
    }
}

void
quat_lideal_create_principal(quat_left_ideal_t *lideal,
                             const quat_alg_elem_t *x,
                             const quat_lattice_t *order,
                             const quat_alg_t *alg)
{
    assert(quat_order_is_maximal(order, alg));
    assert(quat_lattice_contains(NULL, order, x));
    ibz_t norm_n, norm_d;
    ibz_init(&norm_n);
    ibz_init(&norm_d);

    // Multiply order on the right by x
    quat_lattice_alg_elem_mul(&(lideal->lattice), order, x, alg);

    // Reduce denominator. This conserves HNF
    quat_lattice_reduce_denom(&lideal->lattice, &lideal->lattice);

    // Compute norm and check it's integral
    quat_alg_norm(&norm_n, &norm_d, x, alg);
    assert(ibz_is_one(&norm_d));
    ibz_copy(&lideal->norm, &norm_n);

    // Set order
    lideal->parent_order = order;
    ibz_finalize(&norm_n);
    ibz_finalize(&norm_d);
}

/* Construct O<x,N> in coordinates relative to O.
 *
 * The former implementation first built the principal ideal O*x in ambient
 * coordinates and canonicalized it before adding N*O.  In the fixed-width
 * implementation that creates HNF row-reduction values much larger than the
 * final ideal (for example, about 7911 bits in the level-1 signing path while
 * an ibz_t has only 7039 signed magnitude bits).
 *
 * If phi: Z^4 -> O maps standard coordinates to the order basis and C contains
 * the O-coordinates of e_i*x, then
 *
 *     O*x + N*O = phi(C Z^4 + N Z^4).
 *
 * Hence C may first be reduced modulo |N| and the relative HNF can be taken on
 * [ C | |N| I ].  The latter contains |N| Z^4, so |N|^4 is a valid HNF
 * modulus.  This keeps the HNF input tied to N instead of to the much larger
 * determinant of an ambient principal-ideal basis.
 */
int
quat_lideal_create_from_generator_mod_N(quat_left_ideal_t *lideal,
                                         const quat_alg_elem_t *x,
                                         const ibz_t *N,
                                         const quat_lattice_t *order,
                                         const quat_alg_t *alg)
{
    int ok = 0;
    ibz_vec_4_t generators[8], absolute_generators[8], x_relative;
    ibz_mat_4x4_t relative_hnf;
    quat_alg_elem_t x_reduced, order_basis_elem, product;
    quat_left_ideal_t candidate;
    ibz_t abs_N, hnf_mod, absolute_hnf_mod, order_basis_det, relative_index;

    for (int i = 0; i < 8; i++) {
        ibz_vec_4_init(&generators[i]);
        ibz_vec_4_init(&absolute_generators[i]);
    }
    ibz_vec_4_init(&x_relative);
    ibz_mat_4x4_init(&relative_hnf);
    quat_alg_elem_init(&x_reduced);
    quat_alg_elem_init(&order_basis_elem);
    quat_alg_elem_init(&product);
    quat_left_ideal_init(&candidate);
    ibz_init(&abs_N);
    ibz_init(&hnf_mod);
    ibz_init(&absolute_hnf_mod);
    ibz_init(&order_basis_det);
    ibz_init(&relative_index);

    if (quat_alg_elem_is_zero(x) || !quat_lattice_contains(&x_relative, order, x)) {
        goto cleanup;
    }

    ibz_abs(&abs_N, N);
    if (ibz_is_zero(&abs_N)) {
        // O<x,0> is the principal ideal O*x; preserve the public API's old
        // behavior for this non-protocol edge case.
        quat_lideal_create_principal(&candidate, x, order, alg);
        quat_lideal_copy(lideal, &candidate);
        ok = 1;
        goto cleanup;
    }

    const size_t norm_bits = (size_t)ibz_bitsize(&abs_N);
    const size_t algebra_bits = (size_t)ibz_bitsize(&alg->p);
    if (8 * norm_bits + 4 * algebra_bits + 13 > (size_t)(IBZ_BITS - 1)) {
        goto cleanup;
    }

    // O<x,N> depends only on the O-coordinates of x modulo N.  Reducing x
    // before multiplying by the order basis also bounds the algebra products.
    for (int row = 0; row < 4; row++) {
        ibz_mod(&x_relative[row], &x_relative[row], &abs_N);
    }
    ibz_mat_4x4_eval(&x_reduced.coord, &order->basis, &x_relative);
    ibz_copy(&x_reduced.denom, &order->denom);
    quat_alg_normalize(&x_reduced);

    for (int column = 0; column < 4; column++) {
        for (int row = 0; row < 4; row++) {
            ibz_copy(&order_basis_elem.coord[row], &order->basis[row][column]);
        }
        ibz_copy(&order_basis_elem.denom, &order->denom);
        quat_alg_mul(&product, &order_basis_elem, &x_reduced, alg);
        if (!quat_lattice_contains(&generators[column], order, &product)) {
            goto cleanup;
        }
        for (int row = 0; row < 4; row++) {
            ibz_mod(&generators[column][row], &generators[column][row], &abs_N);
        }
        ibz_copy(&generators[4 + column][column], &abs_N);
    }

    ibz_mul(&hnf_mod, &abs_N, &abs_N);
    ibz_mul(&hnf_mod, &hnf_mod, &hnf_mod);
    ibz_mat_4xn_hnf_mod_core(&relative_hnf, 8, generators, &hnf_mod);

    // For a left ideal of a maximal order, its relative index is the square
    // of its norm.  Do not assume the norm is N: non-primitive x may yield a
    // proper divisor (and x=N yields norm N^2).
    ibz_mat_4x4_inv_with_det_as_denom(NULL, &relative_index, &relative_hnf);
    ibz_abs(&relative_index, &relative_index);
    if (!ibz_sqrt(&candidate.norm, &relative_index)) {
        goto cleanup;
    }

    /* Mapping the already-computed relative HNF through order->basis and then
     * running a second generic HNF still permits very large exact column
     * reductions: a level-1 signing case reached 7911 bits.  Instead, map the
     * original, N-sized generators and take the ambient HNF directly:
     *
     *   B (C Z^4 + N Z^4) = (B C) Z^4 + (N B) Z^4.
     *
     * The second summand has determinant |det(B)| N^4, which is therefore a
     * valid HNF modulus and is covered by the budget guard above.
     */
    for (int column = 0; column < 8; column++) {
        ibz_mat_4x4_eval(&absolute_generators[column], &order->basis, &generators[column]);
    }
    ibz_mat_4x4_inv_with_det_as_denom(NULL, &order_basis_det, &order->basis);
    ibz_abs(&order_basis_det, &order_basis_det);
    ibz_mul(&absolute_hnf_mod, &order_basis_det, &hnf_mod);
    ibz_mat_4xn_hnf_mod_core(&candidate.lattice.basis, 8, absolute_generators, &absolute_hnf_mod);
    ibz_copy(&candidate.lattice.denom, &order->denom);
    quat_lattice_reduce_denom(&candidate.lattice, &candidate.lattice);
    candidate.parent_order = order;
    if (!quat_lideal_norm_verify(&candidate)) {
        goto cleanup;
    }

    quat_lideal_copy(lideal, &candidate);
    ok = 1;

cleanup:
    ibz_finalize(&relative_index);
    ibz_finalize(&order_basis_det);
    ibz_finalize(&absolute_hnf_mod);
    ibz_finalize(&hnf_mod);
    ibz_finalize(&abs_N);
    quat_left_ideal_finalize(&candidate);
    quat_alg_elem_finalize(&product);
    quat_alg_elem_finalize(&order_basis_elem);
    quat_alg_elem_finalize(&x_reduced);
    ibz_mat_4x4_finalize(&relative_hnf);
    ibz_vec_4_finalize(&x_relative);
    for (int i = 0; i < 8; i++) {
        ibz_vec_4_finalize(&absolute_generators[i]);
        ibz_vec_4_finalize(&generators[i]);
    }
    return ok;
}

void
quat_lideal_create(quat_left_ideal_t *lideal,
                   const quat_alg_elem_t *x,
                   const ibz_t *N,
                   const quat_lattice_t *order,
                   const quat_alg_t *alg)
{
    assert(quat_order_is_maximal(order, alg));
    assert(!quat_alg_elem_is_zero(x));
    assert(quat_lattice_contains(NULL, order, x));

    // Failure means that the requested inputs exceed the configured fixed
    // precision budget or violate the documented ideal-construction
    // preconditions.  Do not continue with a silently truncated lattice in a
    // release build.
    if (!quat_lideal_create_from_generator_mod_N(lideal, x, N, order, alg)) {
        abort();
    }
}

int
quat_lideal_generator(quat_alg_elem_t *gen, const quat_left_ideal_t *lideal, const quat_alg_t *alg)
{
    ibz_t norm_int, norm_n, gcd, r, q, norm_denom;
    ibz_vec_4_t vec;
    ibz_vec_4_init(&vec);
    ibz_init(&norm_denom);
    ibz_init(&norm_int);
    ibz_init(&norm_n);
    ibz_init(&r);
    ibz_init(&q);
    ibz_init(&gcd);
    int a, b, c, d;
    int found = 0;
    int int_norm = 0;
    while (1) {
        int_norm++;
        for (a = -int_norm; a <= int_norm; a++) {
            for (b = -int_norm + abs(a); b <= int_norm - abs(a); b++) {
                for (c = -int_norm + abs(a) + abs(b); c <= int_norm - abs(a) - abs(b); c++) {
                    d = int_norm - abs(a) - abs(b) - abs(c);
                    ibz_vec_4_set(&vec, a, b, c, d);
                    ibz_vec_4_content(&gcd, &vec);
                    if (ibz_is_one(&gcd)) {
                        ibz_mat_4x4_eval(&(gen->coord), &(lideal->lattice.basis), &vec);
                        ibz_copy(&(gen->denom), &(lideal->lattice.denom));
                        quat_alg_norm(&norm_int, &norm_denom, gen, alg);
                        assert(ibz_is_one(&norm_denom));
                        ibz_div(&q, &r, &norm_int, &(lideal->norm));
                        assert(ibz_is_zero(&r));
                        ibz_gcd(&gcd, &(lideal->norm), &q);
                        found = (0 == ibz_cmp(&gcd, &ibz_const_one));
                        if (found)
                            goto fin;
                    }
                }
            }
        }
    }
fin:;
    ibz_finalize(&r);
    ibz_finalize(&q);
    ibz_finalize(&norm_denom);
    ibz_finalize(&norm_int);
    ibz_finalize(&norm_n);
    ibz_vec_4_finalize(&vec);
    ibz_finalize(&gcd);
    return (found);
}

void
quat_lideal_mul(quat_left_ideal_t *product,
                const quat_left_ideal_t *lideal,
                const quat_alg_elem_t *alpha,
                const quat_alg_t *alg)
{
    assert(quat_order_is_maximal((lideal->parent_order), alg));
    ibz_t norm, norm_d;
    ibz_init(&norm);
    ibz_init(&norm_d);
    quat_lattice_alg_elem_mul(&(product->lattice), &(lideal->lattice), alpha, alg);
    product->parent_order = lideal->parent_order;
    quat_alg_norm(&norm, &norm_d, alpha, alg);
    ibz_mul(&(product->norm), &(lideal->norm), &norm);
    assert(ibz_divides(&(product->norm), &norm_d));
    ibz_div(&(product->norm), &norm, &(product->norm), &norm_d);
    assert(quat_lideal_norm_verify(lideal));
    ibz_finalize(&norm_d);
    ibz_finalize(&norm);
}

void
quat_lideal_add(quat_left_ideal_t *sum, const quat_left_ideal_t *I1, const quat_left_ideal_t *I2, const quat_alg_t *alg)
{
    assert(I1->parent_order == I2->parent_order);
    assert(quat_order_is_maximal((I2->parent_order), alg));
    quat_lattice_add(&sum->lattice, &I1->lattice, &I2->lattice);
    sum->parent_order = I1->parent_order;
    quat_lideal_norm(sum);
}

void
quat_lideal_inter(quat_left_ideal_t *inter,
                  const quat_left_ideal_t *I1,
                  const quat_left_ideal_t *I2,
                  const quat_alg_t *alg)
{
    assert(I1->parent_order == I2->parent_order);
    assert(quat_order_is_maximal((I2->parent_order), alg));
    quat_lattice_intersect(&inter->lattice, &I1->lattice, &I2->lattice);
    inter->parent_order = I1->parent_order;
    quat_lideal_norm(inter);
}

int
quat_lideal_equals(const quat_left_ideal_t *I1, const quat_left_ideal_t *I2, const quat_alg_t *alg)
{
    assert(quat_order_is_maximal((I2->parent_order), alg));
    assert(quat_order_is_maximal((I1->parent_order), alg));
    return (I1->parent_order == I2->parent_order) & (ibz_cmp(&I1->norm, &I2->norm) == 0) &
           quat_lattice_equal(&I1->lattice, &I2->lattice);
}

void
quat_lideal_inverse_lattice_without_hnf(quat_lattice_t *inv, const quat_left_ideal_t *lideal, const quat_alg_t *alg)
{
    assert(quat_order_is_maximal((lideal->parent_order), alg));
    quat_lattice_conjugate_without_hnf(inv, &(lideal->lattice));
    ibz_mul(&(inv->denom), &(inv->denom), &(lideal->norm));
}

// following the implementation of ideal isomorphisms in the code of LearningToSQI's sage
// implementation of SQIsign
void
quat_lideal_right_transporter(quat_lattice_t *trans,
                              const quat_left_ideal_t *lideal1,
                              const quat_left_ideal_t *lideal2,
                              const quat_alg_t *alg)
{
    assert(quat_order_is_maximal((lideal1->parent_order), alg));
    assert(quat_order_is_maximal((lideal2->parent_order), alg));
    assert(lideal1->parent_order == lideal2->parent_order);
    quat_lattice_t inv;
    quat_lattice_init(&inv);
    quat_lideal_inverse_lattice_without_hnf(&inv, lideal1, alg);
    quat_lattice_mul(trans, &inv, &(lideal2->lattice), alg, &(lideal1->norm), &(lideal2->norm));
    quat_lattice_finalize(&inv);
}

void
quat_lideal_right_order(quat_lattice_t *order, const quat_left_ideal_t *lideal, const quat_alg_t *alg)
{
    assert(quat_order_is_maximal((lideal->parent_order), alg));
    quat_lideal_right_transporter(order, lideal, lideal, alg);
}

void
quat_lideal_class_gram(ibz_mat_4x4_t *G, const quat_left_ideal_t *lideal, const quat_alg_t *alg)
{
    quat_lattice_gram(G, &(lideal->lattice), alg);

    // divide by norm · denominator²
    ibz_t divisor, rmd;
    ibz_init(&divisor);
    ibz_init(&rmd);

    ibz_mul(&divisor, &(lideal->lattice.denom), &(lideal->lattice.denom));
    ibz_mul(&divisor, &divisor, &(lideal->norm));

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j <= i; j++) {
            ibz_div(&(*G)[i][j], &rmd, &(*G)[i][j], &divisor);
            assert(ibz_is_zero(&rmd));
        }
    }
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j <= i - 1; j++) {
            ibz_copy(&(*G)[j][i], &(*G)[i][j]);
        }
    }

    ibz_finalize(&rmd);
    ibz_finalize(&divisor);
}

void
quat_lideal_conjugate_without_hnf(quat_left_ideal_t *conj,
                                  quat_lattice_t *new_parent_order,
                                  const quat_left_ideal_t *lideal,
                                  const quat_alg_t *alg)
{
    quat_lideal_right_order(new_parent_order, lideal, alg);
    quat_lattice_conjugate_without_hnf(&(conj->lattice), &(lideal->lattice));
    conj->parent_order = new_parent_order;
    ibz_copy(&(conj->norm), &(lideal->norm));
}

int
quat_order_discriminant(ibz_t *disc, const quat_lattice_t *order, const quat_alg_t *alg)
{
    int ok = 0;
    ibz_t det, sqr, div;
    ibz_mat_4x4_t transposed, norm, prod;
    ibz_init(&det);
    ibz_init(&sqr);
    ibz_init(&div);
    ibz_mat_4x4_init(&transposed);
    ibz_mat_4x4_init(&norm);
    ibz_mat_4x4_init(&prod);
    ibz_mat_4x4_transpose(&transposed, &(order->basis));
    // multiply gram matrix by 2 because of reduced trace
    ibz_mat_4x4_identity(&norm);
    ibz_copy(&(norm[2][2]), &(alg->p));
    ibz_copy(&(norm[3][3]), &(alg->p));
    ibz_mat_4x4_scalar_mul(&norm, &ibz_const_two, &norm);
    ibz_mat_4x4_mul(&prod, &transposed, &norm);
    ibz_mat_4x4_mul(&prod, &prod, &(order->basis));
    ibz_mat_4x4_inv_with_det_as_denom(NULL, &det, &prod);
    ibz_mul(&div, &(order->denom), &(order->denom));
    ibz_mul(&div, &div, &div);
    ibz_mul(&div, &div, &div);
    ibz_div(&sqr, &div, &det, &div);
    ok = ibz_is_zero(&div);
    ok = ok & ibz_sqrt(disc, &sqr);
    ibz_finalize(&det);
    ibz_finalize(&div);
    ibz_finalize(&sqr);
    ibz_mat_4x4_finalize(&transposed);
    ibz_mat_4x4_finalize(&norm);
    ibz_mat_4x4_finalize(&prod);
    return (ok);
}

int
quat_order_is_maximal(const quat_lattice_t *order, const quat_alg_t *alg)
{
    int res;
    ibz_t disc;
    ibz_init(&disc);
    quat_order_discriminant(&disc, order, alg);
    res = (ibz_cmp(&disc, &(alg->p)) == 0);
    ibz_finalize(&disc);
    return (res);
}
