/// \file
///
/// Flint integral polynomials utilities
/// ====================================
/*
   Copyright (C)      2016 Vincent Delecroix
                 2021-2024 Julian Rüth

    This file is part of e-antic

    e-antic is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3.0 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef E_ANTIC_FMPZ_POLY_EXTRA_H
#define E_ANTIC_FMPZ_POLY_EXTRA_H

#include "local.h"

#include <flint/fmpz_poly.h>


#if __FLINT_RELEASE < 30000
#include <arb.h>
#else
#include <flint/arb.h>
#endif


#ifdef __cplusplus
extern "C" {
#endif

/// Set `p` to be a random irreducible polynomial
LIBEANTIC_API void fmpz_poly_randtest_irreducible(fmpz_poly_t p, flint_rand_t state, slong len, mp_bitcnt_t bits);

/// Set the polynomial `p` from the string `s` using `var` as variable name
LIBEANTIC_API int fmpz_poly_set_str_pretty(fmpz_poly_t p, const char * s, const char * var);

/// \name Real root isolation
///@{
/// Precompose the polynomial `pol` by an affine transformation so that the interval
/// [`a`, `b`] becomes the interval [0,1]
LIBEANTIC_API void _fmpz_poly_scale_0_1_fmpq(fmpz * pol, slong len, fmpq_t a, fmpq_t b);

/// Return an upper bound on the bitsize of largest real root of `pol`.
LIBEANTIC_API slong fmpz_poly_positive_root_upper_bound_2exp(const fmpz_poly_t pol);

/// Return an upper bound on the bitsize of largest real root of `(pol, len)`.
LIBEANTIC_API slong _fmpz_poly_positive_root_upper_bound_2exp(const fmpz * pol, slong len);

LIBEANTIC_API slong _fmpz_poly_positive_root_upper_bound_2exp_local_max(const fmpz * pol, slong len);

//// \name Root refinement
///@{
LIBEANTIC_API int _fmpz_poly_newton_step_arb(arb_t res, const fmpz * pol, const fmpz * der, slong len, arb_t a, slong prec);

LIBEANTIC_API int _fmpz_poly_bisection_step_arb(arb_t res, fmpz * pol, slong len, arb_t a, slong prec);

LIBEANTIC_API int fmpz_poly_newton_step_arb(arb_t res, const fmpz_poly_t pol, const fmpz_poly_t der, arb_t a, slong prec);

LIBEANTIC_API void _fmpz_poly_bisection_step_arf(arf_t l, arf_t r, const fmpz * pol, slong len, int sl, int sr, slong prec);

LIBEANTIC_API int fmpz_poly_bisection_step_arb(arb_t res, const fmpz_poly_t pol, arb_t a, slong prec);

LIBEANTIC_API void fmpz_poly_evaluate_at_one(fmpz_t res, fmpz * p, slong len);
///@}

/// \name FLINT, Arb extra
///@{
/// Set `res` to be the polynomial whose coefficients are the absolute values
/// of the ones in `p`
LIBEANTIC_API void fmpz_poly_abs(fmpz_poly_t res, fmpz_poly_t p);

/// Evaluate the polynomial `pol` at the ball `a` and set result in `res`
LIBEANTIC_API void _fmpz_poly_evaluate_arb(arb_t res, const fmpz * pol, slong len, const arb_t a, slong prec);

LIBEANTIC_API void fmpz_poly_evaluate_arb(arb_t b, const fmpz_poly_t pol, const arb_t a, slong prec);

LIBEANTIC_API void _fmpz_poly_evaluate_arf(arf_t res, const fmpz * pol, slong len, const arf_t a, slong prec);

LIBEANTIC_API void fmpz_poly_evaluate_arf(arf_t res, const fmpz_poly_t pol, const arf_t a, slong prec);

/* condition number of a_0 + a_1 x + ... + a_d x^d at a number t is: */
/*  |a_0| + |a_1 t| + ... + |a_d t^d|                                */
LIBEANTIC_API int _fmpz_poly_relative_condition_number_2exp(slong * cond, fmpz * p, slong len, arb_t x, slong prec);
LIBEANTIC_API int fmpz_poly_relative_condition_number_2exp(slong * cond, fmpz_poly_t p, arb_t x, slong prec);
///@}

#ifdef __cplusplus
}
#endif

#endif
