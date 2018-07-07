/*################################################################################
  ##
  ##   Copyright (C) 2017-2018 Odran Bonnet, 
  ##                           Alfred Galichon, 
  ##                           Keith O'Hara, and
  ##                           Matt Shum
  ##
  ##   This file is part of the yogurts package.
  ##
  ##   yogurts is free software: you can redistribute it and/or modify
  ##   it under the terms of the GNU General Public License as published by
  ##   the Free Software Foundation, either version 2 of the License, or
  ##   (at your option) any later version.
  ##
  ##   yogurts is distributed in the hope that it will be useful,
  ##   but WITHOUT ANY WARRANTY; without even the implied warranty of
  ##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ##   GNU General Public License for more details.
  ##
  ##   You should have received a copy of the GNU General Public License
  ##   along with yogurts. If not, see <http://www.gnu.org/licenses/>.
  ##
  ################################################################################*/

#ifndef _SPR_HPP
#define _SPR_HPP

#include <RcppArmadillo.h>
#include <chrono>

namespace yogurts
{
    using uint_t = unsigned int;

    #include "zeroin.hpp"
    #include "tictoc.hpp"

    #include "ipfp.hpp"
    #include "msa.hpp"
}

RcppExport SEXP msa_R(SEXP obs_choice_prob_R, SEXP eps_mat_R, SEXP tol_auction_R);

RcppExport SEXP ipfp_R(SEXP obs_choice_prob_R, SEXP eps_mat_R,
                       SEXP ipfp_temper_R, SEXP err_tol_R, SEXP max_iter_R);

RcppExport SEXP ipfp_lse_R(SEXP obs_choice_prob_R, SEXP eps_mat_R, 
                           SEXP ipfp_temper_R, SEXP err_tol_R, SEXP max_iter_R);

RcppExport SEXP ipfp_arum_R(SEXP obs_choice_prob_R, SEXP eps_mat_R, 
                            SEXP ipfp_temper_R, SEXP err_tol_R, SEXP max_iter_R);

#endif
