# -*- coding: utf-8 -*-
#
# Methods for downloading ERA-Interim data from the ECMWF server for limited
# areas and limited times.
#
# (C) Copyright Stephan Gruber (2013–2017)
#         
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

""" Python tools for downloading ERA-Interim data from ECMWF web services."""

from era_interim import ERAgeneric, ERApl, ERAsa, ERAsf, ERAto, ERAbatch  

__all__ = [
    "ERAgeneric",
    "ERApl",
    "ERAsa",
    "ERAsf",
    "ERAto",
    "ERAbatch",
    "ecmwfapi"        
]
