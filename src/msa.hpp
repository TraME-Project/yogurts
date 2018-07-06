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

//
// auction algorithm

inline
double
max_all(uint_t i, const arma::vec& delta_j, const arma::mat& eps_mat_w0, uint_t& j_max)
{
    arma::vec U_i = delta_j + eps_mat_w0.row(i).t();

    j_max = U_i.index_max();

    return U_i(j_max);
}

inline
arma::vec
bid_round(const arma::vec& assignment_vec, arma::vec& delta_j, 
          const arma::vec& m, const arma::mat& eps_mat_w0,
          const double tq, const int bound_case = 1)
{
    uint_t nbJ0 = eps_mat_w0.n_cols;
    uint_t n_0 = eps_mat_w0.n_rows;

    arma::vec tmp_vec(n_0,1);

#ifdef YOG_USE_OPENMP
    if (n_0 > 100U) {
        #pragma omp parallel for
#endif
        for (uint_t k = 0; k < (n_0 & (~3U)); k+=4)
        {   // unroll loop
            uint_t jj;

            max_all(k,delta_j,eps_mat_w0,jj);
            tmp_vec(k) = jj;

            max_all(k+1,delta_j,eps_mat_w0,jj);
            tmp_vec(k+1) = jj;

            max_all(k+2,delta_j,eps_mat_w0,jj);
            tmp_vec(k+2) = jj;

            max_all(k+3,delta_j,eps_mat_w0,jj);
            tmp_vec(k+3) = jj;
        }

        if (n_0 % 4 > 0)
        {
            uint_t jj;

            max_all(n_0-3,delta_j,eps_mat_w0,jj);
            tmp_vec(n_0-3) = jj;

            max_all(n_0-2,delta_j,eps_mat_w0,jj);
            tmp_vec(n_0-2) = jj;

            max_all(n_0-1,delta_j,eps_mat_w0,jj);
            tmp_vec(n_0-1) = jj;
        }
        
#ifdef YOG_USE_OPENMP
    } else {
        for (uint_t k = 0; k < n_0; k++)
        {
            uint_t jj;

            max_all(k,delta_j,eps_mat_w0,jj);

            tmp_vec(k) = jj;
        }
    }
#endif

    for (uint_t j = 1; j < nbJ0; j++)
    {
        arma::uvec sum_find = arma::find( tmp_vec == j );

        int bin;
        if (bound_case == 0) { // lower
            bin = ( sum_find.n_elem >= m(j) ) ? 1 : 0;
        } else if (bound_case == 1) { // upper
            bin = ( sum_find.n_elem > m(j) ) ? 1 : 0;
        } else { // lower final: bound_case == -1
            bin = ( sum_find.n_elem < m(j) + 0.01 ) ? - 1 : 0;
        }

        delta_j(j) -= bin * tq;
    }

    return tmp_vec;
}

inline
void
msa_auction(const arma::vec& s_vec, const arma::mat& eps_mat, const double tol_auction,
            arma::vec& delta_j_lower, arma::vec& delta_j_upper, comptime_t* time_spent)
{
    clocktime_t start_time;
    if (time_spent) {
        start_time = tic();
    }

    //

    uint_t nbI = eps_mat.n_rows;
    uint_t nbJ = eps_mat.n_cols;

    uint_t nbJ0 = nbJ + 1;

    arma::mat eps_mat_w0 = eps_mat; // append a column of zeros
    eps_mat_w0.insert_cols(0,1);

    arma::vec delta_j_start = -100.0*arma::ones(nbJ0);
    delta_j_start(0) = 0.0;

    arma::vec m(nbJ0);
    m(0) = nbI - arma::accu(s_vec)*nbI;
    m.rows(1,nbJ0-1) = s_vec*nbI;

    // double tol_auction = 1E-06;
    double tq = 50;
    double tq_min = tol_auction;

    //
    // upper bound

    arma::vec delta_j = delta_j_start;
    arma::vec assignment = - arma::ones(nbI,1);
    arma::vec tmp_max_vec(nbI);

    while (tq > tq_min)
    {
        int tmp_check = 0;
        
        while (tmp_check + 0.01 < m(0))
        {
            tmp_max_vec = bid_round(assignment,delta_j,m,eps_mat_w0,tq,1);

            arma::uvec check_vec = arma::find( tmp_max_vec == 0 );
            tmp_check = check_vec.n_elem;
        }

        delta_j.rows(1,nbJ0-1) += tq*nbJ0;

        tq *= 0.25;
    }

    // final step for convergence

    double eta = tq;
    double test_decrease = 1.0;
    arma::vec delta_j_old = delta_j;

    while (test_decrease > 0)
    {
        delta_j.rows(1,nbJ0-1) += eta;
        delta_j_old = delta_j;

        double min_i = 1.0;

        while (min_i > 0.0)
        {
            tmp_max_vec = bid_round(assignment,delta_j,m,eps_mat_w0,tq,1);

            arma::vec sim_market_shares(nbJ0);

            for (uint_t j=0; j < nbJ0; j++)
            {
                arma::uvec check_vec = arma::find( tmp_max_vec == j );
                sim_market_shares(j) = check_vec.n_elem;
            }

            min_i = arma::max( sim_market_shares.rows(1,nbJ0-1) - m.rows(1,nbJ0-1) );
        }

        test_decrease = arma::min( arma::abs( delta_j - delta_j_old ) );

        eta *= 2.0;
    }

    delta_j_upper = delta_j.rows(1,nbJ0-1);

    //
    // lower bound

    delta_j = delta_j_start;

    tq = 50.0;

    while (tq > tq_min)
    {
        int sbloc = 0;
        arma::vec bloc = arma::zeros(nbJ,1);

        delta_j_old = delta_j;
        arma::vec delta_j_tmp = delta_j;

        while (sbloc < nbJ)
        {
            delta_j_old = delta_j;
            bid_round(assignment,delta_j,m,eps_mat_w0,tq,0);

            for (uint_t j=1; j < nbJ0; j++)
            {
                delta_j(j) = (bloc(j-1) > 0) ? delta_j_old(j) : delta_j(j);
            }

            delta_j_tmp = delta_j;
            tmp_max_vec = bid_round(assignment,delta_j_tmp,m,eps_mat_w0,tq,0);

            //

            arma::vec sim_market_shares(nbJ0);

            for (uint_t j=0; j < nbJ0; j++)
            {
                arma::uvec check_vec = arma::find( tmp_max_vec == j );
                sim_market_shares(j) = check_vec.n_elem;
            }

            for (uint_t j=1; j < nbJ0; j++)
            {
                if (sim_market_shares(0) > m(0) && sim_market_shares(j) < m(j))
                {
                    bloc(j-1) = 1;
                }
            }

            sbloc = arma::accu(bloc);
        }

        delta_j.rows(1,nbJ0-1) += tq*nbJ0;

        tq *= 0.25;
    }

    // final step for convergence

    eta = 4.0*tq;

    double test_increase = -1;

    delta_j_old = delta_j;

    while (test_increase <= 0.0)
    {
        delta_j.rows(1,nbJ0-1) = delta_j_old.rows(1,nbJ0-1) - 3*eta;

        delta_j_old = delta_j;

        double min_i = - 1.0;

        while (min_i < 0.0)
        {
            tmp_max_vec = bid_round(assignment,delta_j,m,eps_mat_w0,tq,-1);

            arma::vec sim_market_shares(nbJ0);

            for (uint_t j=0; j < nbJ0; j++)
            {
                arma::uvec check_vec = arma::find( tmp_max_vec == j );
                sim_market_shares(j) = check_vec.n_elem;
            }

            min_i = arma::min( sim_market_shares.rows(1,nbJ0-1) - m.rows(1,nbJ0-1) );
        }

        test_increase = arma::min( delta_j.rows(1,nbJ0-1) - delta_j_old.rows(1,nbJ0-1) );

        eta *= 2.0;
    }

    delta_j_lower = delta_j.rows(1,nbJ0-1);

    //

    if (time_spent) {
        *time_spent = tic() - start_time;
    }
}
