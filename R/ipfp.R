################################################################################
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
################################################################################

ipfp <- function(obs_choice_prob,eps_mat,ipfp_temper=1.0,err_tol=1E-06,max_iter=2000)
{
    res <- .Call("ipfp_R", obs_choice_prob,eps_mat,ipfp_temper,err_tol,max_iter, PACKAGE = "yogurts")

    return(res)
}

ipfp_lse <- function(obs_choice_prob,eps_mat,ipfp_temper=1.0,err_tol=1E-06,max_iter=2000)
{
    res <- .Call("ipfp_lse_R", obs_choice_prob,eps_mat,ipfp_temper,err_tol,max_iter, PACKAGE = "yogurts")

    return(res)
}

ipfp_arum <- function(obs_choice_prob,eps_mat,ipfp_temper=1.0,err_tol=1E-06,max_iter=2000)
{
    res <- .Call("ipfp_arum_R", obs_choice_prob,eps_mat,ipfp_temper,err_tol,max_iter, PACKAGE = "yogurts")

    return(res)
}
