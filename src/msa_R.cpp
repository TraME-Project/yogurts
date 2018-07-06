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

#include "yogurts.hpp"
using namespace Rcpp;
 
SEXP msa_R(SEXP obs_choice_prob_R, SEXP eps_mat_R, SEXP tol_auction_R)
{
    try {
        arma::mat eps_mat( REAL(eps_mat_R), Rf_nrows(eps_mat_R), Rf_ncols(eps_mat_R), false, true );
        arma::vec obs_choice_prob( REAL(obs_choice_prob_R), Rf_length(obs_choice_prob_R), false, true );

        yogurts::comptime_t algo_runtime;
        arma::vec delta_j_lower, delta_j_upper;

        yogurts::msa_auction(obs_choice_prob, eps_mat, as<double>(tol_auction_R),
                             delta_j_lower, delta_j_upper, &algo_runtime);

        double runtime_out = algo_runtime.count();

        return Rcpp::List::create(Rcpp::Named("delta_j_lower") = delta_j_lower,
                                  Rcpp::Named("delta_j_upper") = delta_j_upper,
                                  Rcpp::Named("elapsed_time")  = runtime_out);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "yogurts: C++ exception (unknown reason)" );
    }
    return R_NilValue;
}
