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

msa <- function(obs_choice_prob,eps_mat,tol_auction=1E-06)
{
    res <- .Call("msa_R", obs_choice_prob,eps_mat,tol_auction, PACKAGE = "yogurts")

    return(res)
}
