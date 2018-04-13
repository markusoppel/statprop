******************************************

    Stateprop

    Copyright (c) 1996 - 2018 University of Vienna


    This file is part of stateprop.

    stateprop is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    stateprop is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with stateprop.  If not, see <http://www.gnu.org/licenses/>.


******************************************

Propagation of a wave function on one potential energy surface 
under the action of a laserfield in the state space representation
(i.e. in the basis ofthe eigenfunctions of that potential) 

The time-evolution-operator   -i*H*dt    -i*(H0 - E(t)*mue)dt
                                 e        = e
is evaluated using the split-operator-technique
(see subroutine propagate)


