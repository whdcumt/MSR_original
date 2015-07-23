/*
 *  MSR_original_lib.h
 *  
 *
 *  Created by Catalina Sbert Juan on 12/07/13.
 *
 * Copyright 2013 IPOL Image Processing On Line http://www.ipol.im/
 *
 *
 * This file implements an algorithm possibly linked to the patents:
 *
 *  - US 5991456, "Method of improving a digital image," Issued Nov 23, 1999
 *  - US 6834125, "Method of improving a digital image as a function of its
 *  dynamic range," Issued Dec 21, 2004
 *  - US 6842543 B2, "Method of improving a digital image having white
 *  zones," Issued Jan 11, 2005
 *  - US 8111943, "Smart Image Enhancement Process," Issued Feb 7, 2012
 *  - EP 0901671, "Method of improving a digital image,"
 *  Issued September 3, 2003
 *  - AUS 713076, "Method of improving a digital image,"
 *  Issued February 26, 1998
 *  - WO 1997045809 A1, "Method of improving a digital image," July 4, 2006
 *  - JPO 4036391 B2, "Method of improving a digital image"
 *
 * This file is made available for the exclusive aim of serving as
 * scientific tool to verify the soundness and completeness of the
 * algorithm description. Compilation, execution and redistribution of
 * this file may violate patents rights in certain countries. The
 * situation being different for every country and changing
 * over time, it is your responsibility to determine which patent rights
 * restrictions apply to you before you compile, use, modify, or
 * redistribute this file. A patent lawyer is qualified to make this
 * determination. If and only if they don't conflict with any patent
 * terms, you can benefit from the following license terms attached to this
 * file.
 *
 */

#ifndef _MSR_ORIGINAL_LIB_H
#define _MSR_ORIGINAL_LIB_H

#ifdef __cplusplus
extern "C" {
#endif


double *convolution(double *input, double scale, double *output, size_t nx, size_t ny);


double *MSRetinex(double *out, double *input,double *scale, int nscales, double w, size_t nx, size_t ny);

double *Color_Restoration(double *out,  double *input, double *gray, size_t dim);

double *Gain_offset(double *out, double *input, double G, double b , size_t dim);

#ifdef __cplusplus
}
#endif

#endif /* _MSR_ORIGINAL_LIB_H*/