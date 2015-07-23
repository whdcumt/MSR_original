% Multiscale Retinex
# ABOUT

* Author    : Catalina Sbert  <catalina.sbert@uib.es>
              Ana Belén Petro <anabelen.petro@uib.es>
* Copyright : (C) 2009-2013 IPOL Image Processing On Line http://www.ipol.im/
* License   :

- MSR_original_lib.c, MSR_original.cpp and MSR_original_lib.h
  This file implements an algorithm possibly linked to the patents

  US 5991456, "Method of improving a digital image," Issued Nov 23, 1999
  US 6834125, "Method of improving a digital image as a function of its
  dynamic range," Issued Dec 21, 2004
  US 6842543 B2, "Method of improving a digital image having white
  zones," Issued Jan 11, 2005
  US 8111943, "Smart Image Enhancement Process," Issued Feb 7, 2012
  EP 0901671, "Method of improving a digital image," Issued September 3, 2003
  AUS 713076, "Method of improving a digital image," Issued February 26, 1998
  WO 1997045809 A1, "Method of improving a digital image," July 4, 2006
  JPO 4036391 B2, "Method of improving a digital image"
  These codes are distributed under the terms of the BSD license 
  (see file license.txt)

- All the other files are distributed under the terms of the
  GPLv3 license.


# OVERVIEW

Given an image f, the algorithm applies the Multiscale Retinex with color restoration, based on the work of J.J. Jobson, Z. Rahman and G. A. Woodell 
" A Multiscale Retinex for Bridging the Gap between Color Images and the Human Observation of Scenes" in IEEE Trans. on Image Processing, vol. 6(7) pp. 965-976, 1997.

# REQUIREMENTS

The code is written in  C, the code requires getopt.h and should compile on any system with a C compiler.

The libpng header and libraries are required on the system for
compilation and execution. See http://www.libpng.org/pub/png/libpng.html

The fftw3 header and libraries are required on the system for
compilation and execution. See http://www.fftw.org/

# COMPILATION

Simply use the provided makefile, with the command `make`.

# USAGE

MSR_original [-S S] [-L,L] [-M M] [-H H] [-N N] [-l l] [-R R]  input  MSR_rgb  MSR_gray
	-S  S	number of scales. 
		if scales number is 1 then low scale is the scale to use.
		If the number of scales is 2 then low and medium scales are used 
		(Default: 3)
	-L  L	Low scale (Default: 15)
	-M  M	Medium scale (Default: 80)
	-H  H	High scale (Default: 250)
	-N  N	If 0 final 'canonical' gain/offset; if 1 with final simplest color
		balance (Default: 1)
	-l  l	percentage of saturation on the left (simplest color balance) (Default:1)
	-R  R	percentage  of saturation on the right (simplest color balance) 
		(Default: 1)
	input	 input file
	MSR_rgb	 Multiscale retinex on rgb channels
	MSR_gray	 Multiscale retinex on the gray 




	# ABOUT THIS FILE

Copyright 2009-2013 IPOL Image Processing On Line http://www.ipol.im/
Author: Catalina Sbert <catalina.sbert@uib.es>

Copying and distribution of this file, with or without modification,
are permitted in any medium without royalty provided the copyright
notice and this notice are preserved.  This file is offered as-is,
without any warranty.
