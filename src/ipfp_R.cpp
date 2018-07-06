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

SEXP ipfp_R(SEXP obs_choice_prob_R, SEXP eps_mat_R, SEXP ipfp_temper_R, SEXP err_tol_R, SEXP max_iter_R)
{
    try {
        arma::mat eps_mat( REAL(eps_mat_R), Rf_nrows(eps_mat_R), Rf_ncols(eps_mat_R), false, true );
        arma::vec obs_choice_prob( REAL(obs_choice_prob_R), Rf_length(obs_choice_prob_R), false, true );

        yogurts::comptime_t algo_runtime;
        arma::vec delta_j_out;

        yogurts::ipfp(obs_choice_prob, eps_mat, as<double>(ipfp_temper_R),
                      as<double>(err_tol_R), as<yogurts::uint_t>(max_iter_R),
                      delta_j_out, &algo_runtime);

        double runtime_out = algo_runtime.count();

        return Rcpp::List::create(Rcpp::Named("delta_j") = delta_j_out,
                                  Rcpp::Named("elapsed_time") = runtime_out);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "yogurts: C++ exception (unknown reason)" );
    }
    return R_NilValue;
}

SEXP ipfp_lse_R(SEXP obs_choice_prob_R, SEXP eps_mat_R, SEXP ipfp_temper_R, SEXP err_tol_R, SEXP max_iter_R)
{
    try {
        arma::mat eps_mat( REAL(eps_mat_R), Rf_nrows(eps_mat_R), Rf_ncols(eps_mat_R), false, true );
        arma::vec obs_choice_prob( REAL(obs_choice_prob_R), Rf_length(obs_choice_prob_R), false, true );

        double temper_par = as<double>(ipfp_temper_R);
        double err_tol    = as<double>(err_tol_R);

        yogurts::uint_t max_iter = as<yogurts::uint_t>(max_iter_R);

        //

        yogurts::comptime_t algo_runtime;
        double runtime_out = 0.0;
        arma::vec delta_j_out;

        arma::vec delta_j_init = arma::zeros(obs_choice_prob.n_elem,1);

        if (temper_par > 0.1) {
            yogurts::ipfp_lse(delta_j_init, obs_choice_prob, eps_mat,
                              temper_par, err_tol, max_iter,
                              delta_j_out, &algo_runtime);

            runtime_out = algo_runtime.count();
        } else {
            std::cout << "ipfp: beginning test run with temperature value = 0.1." << std::endl;
            yogurts::ipfp_lse(delta_j_init, obs_choice_prob, eps_mat,
                              0.1, err_tol, max_iter,
                              delta_j_out, &algo_runtime);

            runtime_out = algo_runtime.count();

            //

            delta_j_init = std::move(delta_j_out);

            std::cout << "\nipfp: beginning main run with temperature value = " << temper_par << "." << std::endl;
            yogurts::ipfp_lse(delta_j_init, obs_choice_prob, eps_mat,
                              temper_par, err_tol, max_iter,
                              delta_j_out, &algo_runtime);

            runtime_out += algo_runtime.count();
        }


        // if (temper_par > 1E-04)
        // {
        //     yogurts::ipfp_lse(obs_choice_prob, eps_mat, temper_par,
        //                       err_tol, max_iter,
        //                       delta_j_out, &algo_runtime);
        // } else {
        //     yogurts::ipfp_lse_long(obs_choice_prob, eps_mat, temper_par,
        //                            err_tol, max_iter,
        //                            delta_j_out, &algo_runtime);
        // }

        return Rcpp::List::create(Rcpp::Named("delta_j") = delta_j_out,
                                  Rcpp::Named("elapsed_time") = runtime_out);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "yogurts: C++ exception (unknown reason)" );
    }
    return R_NilValue;
}

SEXP ipfp_arum_R(SEXP obs_choice_prob_R, SEXP eps_mat_R, SEXP ipfp_temper_R, SEXP err_tol_R, SEXP max_iter_R)
{
    try {
        arma::mat eps_mat( REAL(eps_mat_R), Rf_nrows(eps_mat_R), Rf_ncols(eps_mat_R), false, true );
        arma::vec obs_choice_prob( REAL(obs_choice_prob_R), Rf_length(obs_choice_prob_R), false, true );

        yogurts::comptime_t algo_runtime;
        arma::vec delta_j_out;

        yogurts::ipfp_arum(obs_choice_prob, eps_mat, as<double>(ipfp_temper_R),
                           as<double>(err_tol_R), as<yogurts::uint_t>(max_iter_R),
                           delta_j_out, &algo_runtime);

        double runtime_out = algo_runtime.count();

        return Rcpp::List::create(Rcpp::Named("delta_j") = delta_j_out,
                                  Rcpp::Named("elapsed_time") = runtime_out);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "yogurts: C++ exception (unknown reason)" );
    }
    return R_NilValue;
}
