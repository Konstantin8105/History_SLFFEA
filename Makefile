# SLFFEA Makefile
# Version:  1.3
# Copyright (C) 1999, 2000, 2001, 2002  San Le (slffea@yahoo.com)
# http://slffea.sourceforge.net/
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Library General Public License for more details.
#
# You should have received a copy of the GNU Library General Public
# License along with this library; if not, write to the Free
# Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.



default:
	@echo ""
	@echo "The options to build SLFFEA are:"
	@echo ""
	@echo "  make everything        You have Mesa-3.0 or OpenGL and want to build everything"
	@echo "  make science           You just want to build the scientific codes "
	@echo "  make gui               After building the scientific codes,"
	@echo "                                 you now want to build the GUI"
	@echo "  make clean             Remove every executable"
	@echo "  make clean_science     Remove science executables"
	@echo "  make clean_gui         Remove GUI executables"
	@echo ""

everything:
	cd common ;  make
	cd beam/beam ; make 
	cd brick/brick ;  make
	cd brick/brick2 ;  make
	cd brick/brick_nl ;  make 
	cd plate/plate ; make 
	cd quad/quad ;  make
	cd shell/shell ; make
	cd truss/truss ;  make 
	cd truss/truss_nl ;  make 
	cd tetra/tetra ;  make 
	cd tri/tri ;  make 
	cd wedge/wedge ;  make 
#	cd truss/fortran ;  make
#	cd truss/geotruss ;  make
	cd common_gr ;  make
	cd beam/beam_gr ; make 
	cd brick/brick_gr ;  make 
	cd brick/brick2gr ;  make 
	cd plate/plate_gr ; make 
	cd quad/quad_gr ;  make
	cd shell/shell_gr ; make
	cd truss/truss_gr ;  make 
	cd tetra/tetra_gr ;  make 
	cd tri/tri_gr ;  make 
	cd wedge/wedge_gr ;  make 

# This makes the scientific FEM executables only.

science:
	cd common ;  make
	cd beam/beam ; make 
	cd brick/brick ;  make
	cd brick/brick2 ;  make
	cd brick/brick_nl ;  make 
	cd plate/plate ; make 
	cd quad/quad ;  make
	cd shell/shell ; make
	cd truss/truss ;  make 
	cd truss/truss_nl ;  make 
	cd tetra/tetra ;  make 
	cd tri/tri ;  make 
	cd wedge/wedge ;  make 

# This makes the graphics FEM executables only, but
# the science must have been made already.

gui:
	cd common_gr ;  make
	cd beam/beam_gr ; make 
	cd brick/brick_gr ;  make 
	cd brick/brick2gr ;  make 
	cd plate/plate_gr ; make 
	cd quad/quad_gr ;  make
	cd shell/shell_gr ; make
	cd truss/truss_gr ;  make 
	cd tetra/tetra_gr ;  make 
	cd tri/tri_gr ;  make 
	cd wedge/wedge_gr ;  make 

clean:
	cd common ;  make clean
	cd beam/beam ; make clean
	cd brick/brick ;  make clean
	cd brick/brick2 ;  make clean
	cd brick/brick_nl ;  make clean
	cd plate/plate ; make clean
	cd quad/quad ;  make clean
	cd shell/shell ; make clean
	cd truss/truss ;  make clean
	cd truss/truss_nl ;  make clean
	cd tetra/tetra ;  make clean
	cd tri/tri ;  make clean
	cd wedge/wedge ;  make clean
#	cd truss/fortran ;  make clean
#	cd truss/geotruss ;  make clean
	cd common_gr ;  make clean
	cd beam/beam_gr ; make clean
	cd brick/brick_gr ;  make clean
	cd brick/brick2gr ;  make clean
	cd plate/plate_gr ; make clean
	cd quad/quad_gr ;  make clean
	cd shell/shell_gr ; make clean
	cd truss/truss_gr ;  make clean
	cd tetra/tetra_gr ; make clean
	cd tri/tri_gr ; make clean
	cd wedge/wedge_gr ;  make clean

# The lines below remove the science or the graphics or all
# of the executables.

clean_science:
	cd common ; make clean
	cd beam/beam ; make clean
	cd brick/brick ; make clean
	cd brick/brick2 ; make clean
	cd brick/brick_nl ; make clean
	cd plate/plate ; make clean
	cd quad/quad ; make clean
	cd shell/shell ; make clean
	cd truss/truss ; make clean
	cd truss/truss_nl ; make clean
	cd tetra/tetra ; make clean
	cd tri/tri ; make clean
	cd wedge/wedge ; make clean

clean_gui:
	cd common_gr ; make clean
	cd beam/beam_gr ; make clean
	cd brick/brick_gr ; make clean
	cd brick/brick2gr ; make clean
	cd plate/plate_gr ; make clean
	cd quad/quad_gr ; make clean
	cd shell/shell_gr ; make clean
	cd truss/truss_gr ; make clean
	cd tetra/tetra_gr ; make clean
	cd tri/tri_gr ; make clean
	cd wedge/wedge_gr ; make clean
