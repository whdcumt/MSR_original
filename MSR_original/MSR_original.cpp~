/*
 *
 * Copyright 2013 IPOL Image Processing On Line http://www.ipol.im/
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
 * @file MSR_original.cpp
 * @brief Multiscale Retinex with color restoration
 *
 * Algorithm based on the original work of Jobson et al.
 * "A multiscale Retinex for bridging the gap between color images
 * and the human observations of scenes"
 *
 * Read/write operations (png format) make use
 * of io_png.c and io_png.h, by Nicolas Limare
 *
 * @author Catalina Sbert <catalina.sbert@uib.es/>
 * @author Ana Belén Petro <anabelen.petro@uib.es/>
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "io_png.h"
#include "MSR_original_lib.h"
#include "parser.h"
#include "auxiliary_lib.h"

#define min(a,b) ((a) < (b) ? (a):(b))

using namespace std;

int main(int argc,  char **argv)
{


    int nscales, channels;
    double w, scale[3];
    unsigned char *data_in, *data_outC, *data_outG;
    size_t nx,ny, image_size, nc;
    double *R, *G,*B;
    double *Rout, *Gout, *Bout;
    double *gray, *grayout;
    float s1, s2;
    unsigned char option;



    std::vector <OptStruct *> options;
    OptStruct oS  = {"S:", 0, "3", NULL, "number of scales. if scales number is 1 then low scale is the scale to use with the value you want. If the number of scales is 2 then low and medium scales are used"};
    options.push_back(&oS);
    OptStruct oL  = {"L:", 0, "15", NULL, "Low scale"};
    options.push_back(&oL);
    OptStruct oM  = {"M:", 0, "80", NULL, "Medium scale"};
    options.push_back(&oM);
    OptStruct oH  = {"H:", 0, "250",NULL, "High scale"};
    options.push_back(&oH);
    OptStruct oN  = {"N:", 0, "1", NULL, "If 0 final 'canonical' gain/offset; if 1 with final simplest color balance"};
    options.push_back(&oN);
    OptStruct ol  = {"l:", 0, "1", NULL, "percentage of saturation on the left (simplest color balance)"};
    options.push_back(&ol);
    OptStruct oR  = {"R:", 0, "1", NULL, "percentage  of saturation on the right (simplest color balance)"};
    options.push_back(&oR);

    std::vector<ParStruct *> pparameters;
    ParStruct pinput = {"input", NULL, "input file"};
    pparameters.push_back(&pinput);
    ParStruct pMSR_channels= {"MSR_channels", NULL,
                              "Multiscale retinex on rgb channels"
                             };
    pparameters.push_back(&pMSR_channels);
    ParStruct pMSR_gray = {"MSR_gray", NULL,
                           "Multiscale retinex on the gray "
                          };
    pparameters.push_back(&pMSR_gray);

    if (!parsecmdline("MSR_original","Multiscale Retinex with color restoration",
                      argc, argv, options, pparameters))
        return EXIT_FAILURE;
    

    /* read the PNG input image and the options parameters*/

    nscales=atoi(oS.value);

    if(nscales <= 0 || nscales > 3)
    {
        printf("nscales must be 1, 2 or 3\n");
        return EXIT_FAILURE;
    }
    else if(nscales == 1) scale[0]=atof(oL.value);
    else if(nscales == 2)
    {
        scale[0]=atof(oL.value);
        scale[1]=atof(oM.value);
    }
    else
    {
        scale[0]=atof(oL.value);
        scale[1]=atof(oM.value);
        scale[2]=atof(oH.value);
    }

    option=(unsigned char) atoi(oN.value);

    s1=atof(ol.value);
    s2=atof(oR.value);

    data_in=io_png_read_u8(pinput.value, &nx, &ny,&nc);

    if(nc >= 3) channels=3;
    else channels=1;

    image_size=nx*ny;
    
    
    /* memory for the output*/
    
    R=(double*) malloc(image_size*sizeof(double));
    G=(double*) malloc(image_size*sizeof(double));
    B=(double*) malloc(image_size*sizeof(double));
    Rout=(double*) malloc(image_size*sizeof(double));
    Gout=(double*) malloc(image_size*sizeof(double));
    Bout=(double*) malloc(image_size*sizeof(double));
    gray=(double*) malloc(image_size*sizeof(double));
    grayout=(double*) malloc(image_size*sizeof(double));
    data_outC=(unsigned char*)malloc(3*image_size*sizeof(unsigned char));
    data_outG=(unsigned char*)malloc(3*image_size*sizeof(unsigned char));

    
    /* define the three color channels and intensity as double arrays with a translation of 1 for the log function*/

    input_rgb(data_in, R, G, B, image_size);
    gray_intensity(gray, R, G, B, image_size);
    
    
    /* weight of each scale*/

    w=(double)(1./nscales);

    
    /* MSR on RGB channels*/

    MSRetinex( Rout,  R,  scale, nscales, w, nx,ny);
    MSRetinex( Gout,  G,  scale, nscales, w, nx,ny);
    MSRetinex( Bout,  B,  scale, nscales, w, nx,ny);
    

    /* Color restoration  */

    Color_Restoration(Rout, R, gray, image_size);
    Color_Restoration(Gout, G, gray, image_size);
    Color_Restoration(Bout, B, gray, image_size);


    if(option == 0)  /* Gain/Offset*/
    {
        Gain_offset(Rout, Rout, 30, -6, image_size);
        Gain_offset(Gout, Gout, 30, -6, image_size);
        Gain_offset(Bout, Bout, 30, -6, image_size);
    }
    else  /* Simplest color balance with s1% and s2% of saturation*/
    {
        simplest_color_balance(Rout, Rout, image_size, s1, s2);
        simplest_color_balance(Gout, Gout, image_size, s1, s2);
        simplest_color_balance(Bout, Bout, image_size, s1, s2);
    }
    
    
    /*write the output*/

    rgb_output(Rout, Gout, Bout, data_outC, image_size);
    io_png_write_u8(pMSR_channels.value, data_outC, nx,ny, nc);
    

    /* MSR on the gray*/

    MSRetinex(grayout, gray, scale, nscales, w, nx, ny);
    simplest_color_balance(grayout, grayout, image_size, s1, s2);
    compute_color_from_grayscale(Rout, Gout, Bout, R, G, B, gray, grayout,
                                 image_size);
    rgb_output(Rout, Gout, Bout, data_outG, image_size);
    io_png_write_u8(pMSR_gray.value, data_outG, nx,ny, nc);
    

    /*Free the memory*/
    
    free(R); free(G); free(B);
    free(Rout); free(Gout); free(Bout);
    free(gray); free(grayout);
    free(data_outG); free(data_outC);
    
    return EXIT_SUCCESS;

}

