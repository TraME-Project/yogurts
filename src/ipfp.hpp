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
// general ipfp

using uint_t = unsigned int;

inline
double
Ucal(uint_t i, uint_t j, const double delta_val, const arma::mat& eps_mat)
{
    return delta_val + eps_mat(i,j);
}

inline
arma::vec
Ucal(uint_t i, const double delta_val, const arma::mat& eps_mat,
     const bool row_ord = true)
{
    if (row_ord) {
        return delta_val + eps_mat.row(i).t();
    } else {
        return delta_val + eps_mat.col(i);
    }
}

inline
arma::vec
Ucal(uint_t i, const arma::vec& delta_vec, const arma::mat& eps_mat,
     const bool row_ord = true)
{
    if (row_ord) {
        return delta_vec + eps_mat.row(i).t();
    } else {
        return delta_vec + eps_mat.col(i);
    }
}

inline
arma::mat
Ucal(const arma::vec& delta_vec, const arma::mat& eps_mat)
{
    return arma::repmat(delta_vec.t(),eps_mat.n_rows,1) + eps_mat;
}

//

struct root_data_t
{
    uint_t j;
    uint_t N;
    double ipfp_temper_par;

    arma::vec s_vec;
    arma::vec U_i;
    arma::mat eps_mat;
};

inline
double
root_fn(const double val_inp, void* opt_data)
{
    root_data_t* d = reinterpret_cast<root_data_t*>(opt_data);

    return arma::accu( arma::exp( Ucal(d->j,val_inp,d->eps_mat,false) / d->ipfp_temper_par ) / d->U_i ) - d->s_vec(d->j);
}

inline
void
ipfp(const arma::vec& s_vec, const arma::mat& eps_mat,
     const double ipfp_temper_par, const double err_tol, const uint_t max_iter,
     arma::vec& delta_j_out, comptime_t* time_spent)
{
    clocktime_t start_time;
    if (time_spent) {
        start_time = tic();
    }

    //

    uint_t nbI = eps_mat.n_rows;
    uint_t nbJ = eps_mat.n_cols;

    arma::vec delta_j = arma::zeros(nbJ);
    arma::vec delta_j_next = delta_j;

    double delta_min = -10.0;
    double delta_max = 10.0;

    root_data_t root_data;
    root_data.ipfp_temper_par = ipfp_temper_par;
    root_data.s_vec = s_vec;
    root_data.eps_mat = eps_mat;

    //
    // begin loop

    double err = 2*err_tol;
    uint_t iter = 0;

    while (err > err_tol && iter < max_iter)
    {
        iter++;

        root_data.U_i = nbI * ( 1.0 + arma::sum(arma::exp(Ucal(delta_j,eps_mat)),1) / ipfp_temper_par );

        for (uint_t j=0; j < nbJ; j++)
        {
            root_data.j = j;

            delta_j_next(j) = zeroin(delta_min,delta_max,root_fn,&root_data);
        }

        err = arma::max( arma::abs(delta_j_next - delta_j) / ( arma::abs(delta_j_next) + arma::abs(delta_j) ) );

        delta_j = delta_j_next;
    }

    delta_j_out = delta_j_next;

    //

    if (time_spent) {
        *time_spent = tic() - start_time;
    }
}

//
// Log-Sum-Exp

inline
double
root_fn_lse(const double val_inp, void* opt_data)
{
    root_data_t* d = reinterpret_cast<root_data_t*>(opt_data);

    arma::vec diff_vec = Ucal(d->j,val_inp,d->eps_mat,false) - d->U_i;
    double max_diff = diff_vec.max();

    double term_1 = d->ipfp_temper_par * std::log( arma::accu(arma::exp( (diff_vec - max_diff) / d->ipfp_temper_par )) );
    double term_2 = max_diff - d->ipfp_temper_par*std::log(d->N*d->s_vec(d->j));

    return term_1 + term_2;
}

inline
void
ipfp_lse(const arma::vec& delta_j_init, const arma::vec& s_vec, const arma::mat& eps_mat,
         const double ipfp_temper_par, const double err_tol, const uint_t max_iter,
         arma::vec& delta_j_out, comptime_t* time_spent)
{
    clocktime_t start_time;
    if (time_spent) {
        start_time = tic();
    }

    //

    uint_t nbI = eps_mat.n_rows;
    uint_t nbJ = eps_mat.n_cols;

    arma::vec delta_j = delta_j_init;
    arma::vec delta_j_next = delta_j;

    arma::mat U_mat(nbI,nbJ+1);
    U_mat.col(0).zeros(); // add the outside option

    double delta_min = -10.0;
    double delta_max = 10.0;

    root_data_t root_data;
    root_data.N = nbI;
    root_data.ipfp_temper_par = ipfp_temper_par;
    root_data.s_vec = s_vec;
    root_data.eps_mat = eps_mat;

    //
    // begin loop

    double err = 2*err_tol;
    uint_t iter = 0;

    while (err > err_tol && iter < max_iter)
    {
        iter++;

        U_mat.cols(1,nbJ) = Ucal(delta_j,eps_mat);

        arma::vec maxU = arma::max(U_mat,1);

        root_data.U_i = maxU + ipfp_temper_par * arma::log( arma::sum(arma::exp((U_mat - arma::repmat(maxU,1,U_mat.n_cols)) / ipfp_temper_par),1) );

        for (uint_t j=0; j < nbJ; j++)
        {
            root_data.j = j;

            delta_j_next(j) = zeroin(delta_min,delta_max,root_fn_lse,&root_data);
        }

        err = arma::max( arma::abs(delta_j_next - delta_j) / ( arma::abs(delta_j_next) + arma::abs(delta_j) ) );

        delta_j = delta_j_next;
    }

    delta_j_out = delta_j_next;

    //

    std::cout << "ipfp: error = " << err << std::endl;
    std::cout << "      iter  = " << iter << std::endl;

    //

    if (time_spent) {
        *time_spent = tic() - start_time;
    }
}

//
// additional precision to deal with small temperature values

template<typename T>
struct root_data_long_t
{
    uint_t j;
    uint_t N;
    T ipfp_temper_par;

    arma::vec s_vec;
    arma::vec U_i;
    arma::mat eps_mat;
};

template<typename T>
inline
T
root_fn_lse_long(const T val_inp, void* opt_data)
{
    root_data_long_t<T>* d = reinterpret_cast<root_data_long_t<T>*>(opt_data);

    //

    uint_t n_vec_vals = d->U_i.n_elem;
    const T temper_par = d->ipfp_temper_par; // temperature parameter

    T* diff_vals = new T[n_vec_vals]; // Keith: fix memory alignment here

    for (uint_t jj=0U; jj < n_vec_vals; jj++)
    {
        diff_vals[jj] = val_inp + d->eps_mat(jj,d->j) - d->U_i(jj);
    }

    T max_diff = diff_vals[0];

    for (uint_t jj=1U; jj < n_vec_vals; jj++)
    {
        if (max_diff < diff_vals[jj])
        {
            max_diff = diff_vals[jj];
        }
    }

    for (uint_t jj=0U; jj < n_vec_vals; jj++)
    {
        diff_vals[jj] -= max_diff;
    }
    
    //

    T running_sum = T(0.0);

    for (uint_t jj=0U; jj < n_vec_vals; jj++)
    {
        running_sum += std::exp(diff_vals[jj] / temper_par);
    }

    //

    T term_1 = temper_par * std::log( running_sum );
    T term_2 = max_diff - temper_par*std::log(d->N*d->s_vec(d->j));

    //

    delete[] diff_vals;

    return term_1 + term_2;
}

inline
void
ipfp_lse_long(const arma::vec& delta_j_init, const arma::vec& s_vec, const arma::mat& eps_mat,
              const long double ipfp_temper_par, const long double err_tol, const uint_t max_iter,
              arma::vec& delta_j_out, comptime_t* time_spent)
{
    clocktime_t start_time;
    if (time_spent) {
        start_time = tic();
    }

    //
    printf("ipfp: using long precision values.");

    uint_t nbI = eps_mat.n_rows;
    uint_t nbJ = eps_mat.n_cols;

    arma::vec delta_j = arma::zeros(nbJ);
    arma::vec delta_j_next = delta_j;

    arma::mat U_mat(nbI,nbJ+1);
    U_mat.col(0).zeros(); // add the outside option

    // (long double)* exp_vals = new long double[nbI];

    long double delta_min = -10.0;
    long double delta_max = 10.0;

    root_data_long_t<long double> root_data;
    root_data.N = nbI;
    root_data.ipfp_temper_par = ipfp_temper_par;
    root_data.s_vec = s_vec;
    root_data.eps_mat = eps_mat;
    root_data.U_i = arma::zeros(nbI,1);

    //
    // begin loop

    double err = 2*err_tol;
    uint_t iter = 0;

    while (err > err_tol && iter < max_iter)
    {
        iter++;

        U_mat.cols(1,nbJ) = Ucal(delta_j,eps_mat);

        arma::vec maxU = arma::max(U_mat,1);

        //

        arma::mat diff_mat = U_mat - arma::repmat(maxU,1,U_mat.n_cols);

        for (uint_t i=0; i < nbI; i++)
        {
            long double sum_val = 0.0L;
            for (uint_t j=0; j < diff_mat.n_cols; j++)
            {
                sum_val += std::exp( (long double)(diff_mat(i,j)) / ipfp_temper_par );
            }

            // exp_vals[i] = ipfp_temper_par * std::log(sum_val);
            root_data.U_i(i) = maxU(i) + ipfp_temper_par * std::log(sum_val);
        }

        //

        // root_data.U_i = maxU + ipfp_temper_par * arma::log( arma::sum(arma::exp((U_mat - arma::repmat(maxU,1,U_mat.n_cols)) / ipfp_temper_par),1) );

        for (uint_t j=0; j < nbJ; j++)
        {
            root_data.j = j;

            delta_j_next(j) = zeroin_long<long double>(delta_min,delta_max,root_fn_lse_long<long double>,&root_data);
        }

        err = arma::max( arma::abs(delta_j_next - delta_j) / ( arma::abs(delta_j_next) + arma::abs(delta_j) ) );

        delta_j = delta_j_next;
    }

    delta_j_out = delta_j_next;

    //

    std::cout << "ipfp: error = " << err << std::endl;
    std::cout << "      iter  = " << iter << std::endl;

    //

    if (time_spent) {
        *time_spent = tic() - start_time;
    }
}


//
// IPFP-ARUM

inline
void
ipfp_arum(const arma::vec& s_vec, const arma::mat& eps_mat,
          const double ipfp_temper_par, const double err_tol, const uint_t max_iter,
          arma::vec& delta_j_out, comptime_t* time_spent)
{
    clocktime_t start_time;
    if (time_spent) {
        start_time = tic();
    }

    //

    uint_t nbI = eps_mat.n_rows;
    uint_t nbJ = eps_mat.n_cols;
    uint_t nbJ0 = nbJ + 1;

    //

    arma::vec s_vec0(nbJ0);
    s_vec0(0) = 1.0 - arma::accu(s_vec);               // outside option in row 0
    s_vec0.rows(1,nbJ0-1) = s_vec;

    arma::mat K_mat = arma::exp(eps_mat / ipfp_temper_par);
    K_mat = arma::join_rows(arma::ones(nbI,1),K_mat);  // outside option in column 0

    arma::vec A_vec = arma::ones(nbI,1);
    arma::vec B_vec(nbJ0,1);

    //
    // begin loop

    double err = 2*err_tol;
    uint_t iter = 0;

    while (err > err_tol && iter < max_iter)
    {
        iter++;

        B_vec = s_vec0 / (K_mat.t() * A_vec);
        A_vec = 1.0 / (nbI * K_mat*B_vec);

        err = arma::max( arma::abs(B_vec % (K_mat.t()*A_vec) - 1.0) ) ;
    }

    delta_j_out = ipfp_temper_par * arma::log( B_vec.rows(1,nbJ0-1) / B_vec(0) );

    //

    if (time_spent) {
        *time_spent = tic() - start_time;
    }
}
