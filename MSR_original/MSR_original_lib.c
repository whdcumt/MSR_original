/*
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

/**
 *  @file MSR_original_lib.c
 *
 *  @brief Libraries using in the MSR_original.cpp
 *
 *
 *  @author Catalina Sbert Juan
 */



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include "MSR_original_lib.h"


#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif
#define PI2  6.283185307179586  /* 2*pi*/


/**
 * @brief Convolution with a Gaussian kernel using FFT.
 *
 *
 * @param input double array
 * @param scale the size  of the gaussian kernel
 * @param nx x-size of the array
 * @param ny y-size of the array
 *
 * @return output the convolved array
 */

double *convolution(double *input, double scale, double *output,
                    size_t nx, size_t ny)
{
    double *out;
    fftw_plan p;
    int image_size, image_size4;
    int i,j,index;
    double sigma,normx, normy;

    out = (double*) fftw_malloc(sizeof(double) * (nx*ny));

    /*compute the Fourier transform of the input data*/

    p= fftw_plan_r2r_2d((int)ny, (int)nx, input, out, FFTW_REDFT10,
                        FFTW_REDFT10,FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);

    /*define the gaussian constants for the convolution*/

    sigma=scale*scale/2.;
    normx=M_PI/(double)nx;
    normy=M_PI/(double) ny;
    normx*=normx;
    normy*=normy;

    image_size=(int)nx * (int)ny;
    image_size4=4*image_size;

    for(j=0; j<(int)ny; j++)
    {
        index=j*(int)nx;
        for(i=0; i<(int)nx; i++)
            out[i+index]*=exp((double)(-sigma)*(normx*i*i+normy*j*j));
    }

    /*compute the Inverse Fourier transform */

    p=fftw_plan_r2r_2d((int)ny, (int)nx, out, output, FFTW_REDFT01,
                       FFTW_REDFT01, FFTW_ESTIMATE);
    fftw_execute(p);

    for(index=0; index<image_size; index++)
        output[index]/=image_size4;

    fftw_destroy_plan(p);
    fftw_free(out);

    return output;
}

/**
 * @brief The main part of the Multiscale Retinex
 *
 * @f$ MSRout= \sum w (\log(input)-\log(input* G_\sigma))\f$
 *
 * @param input input color channel
 * @param scale[nscales] the  scales for the convolution
 * @param nscales number of scales
 * @param w the weight for each scale
 * @param nx x-size of the image
 * @param ny y-size of the image
 *
 * @return out output of the multiscale retinex
 */

double *MSRetinex(double *out, double *input, double *scale, int nscales,
                  double w, size_t nx, size_t ny)
{
    int i, image_size, n;
    double *pas;

    image_size=(int)nx* (int)ny;

    pas=(double*) malloc(image_size*sizeof(double));

    /* initialization of the output*/
    
    for(i=0; i<image_size; i++)
        out[i]=0.;

    /* Compute Retinex output*/

    for(n=0; n<nscales; n++)
    {
        convolution(input,scale[n],pas,nx,ny);
        for(i=0; i<(int)image_size; i++)
            out[i]+=w*(log(input[i])-log(pas[i]));
    }

    free(pas);

    return out;
}

/**
 * @brief Color restoration for the multiscale Retinex
 *
 * Consists of multiplying the output of the multiscale Retinex by a
 * color restoration function (CRF).
 *
 * @f$ CRF= (\log(125 input)-\log(gray))\f$
 *
 * @f$ MSRCRout=CRF MSRout\f$
 *
 * @param input input color channel
 * @param gray intensity channel of the input color image
 * @param image_size size of the image
 *
 * @return out output of the multiscale retinex with color restoration
 */


double *Color_Restoration(double *out,  double *input, double *gray,
                          size_t image_size)
{

    int i;
    double A;

    for(i=0; i<(int)image_size; i++)
    {
        A=log(3*gray[i]);
        out[i]*=(log(125.*input[i])-A);
    }

    return out;

}

/**
 * @brief Gain/offset
 * is a linear transformation to transform the logarithmic domain into the
 * display domain [0,255]
 *
 * @f$ out=G(input-b) \f$
 *
 * @param input input color channel
 * @param G the gain
 * @param b the offset
 * @param image_size size of the image
 *
 * @return out  output color channel
 */


double *Gain_offset(double *out, double *input, double G, double b,
                    size_t image_size)
{
    int i;

    for(i=0; i<(int)image_size; i++)
    {
        out[i]=G*(input[i]-b);
        if(out[i] < 0.) out[i]=0.;
        if(out[i] > 255.) out[i]=255.;
    }
    return out;
}
