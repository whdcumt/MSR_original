/*
 *  auxiliary_lib.c
 *
 *
 *  Created by Catalina Sbert Juan on 13/11/12.
 *  Copyright 2012 Universitat de les Illes Balears. All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MIN(a,b) ((a) < (b) ? (a):(b))
#define MAX(a,b) ((a) < (b) ? (b):(a))


/**
 * @brief Transform an unsigned char array into three double arrays
 *
 * @param data_in input unsigned char array
 * @return R red double array
 * @return G green double array
 * @return B blue double array
 */

double *input_rgb(unsigned char *input, double *R, double *G, double *B,
                  size_t size)
{

    int n, image_size;

    image_size=(int)size;

    for (n=0; n<image_size; n++)
    {
        R[n]=(double)input[n]+1.;
        G[n]=(double)input[image_size+n]+1.;
        B[n]=(double)input[2*image_size+n]+1.;
    }

    return R;
    return G;
    return B;


}
/**
 *@brief Transform  three double arrays into an unsigned char array
 *
 *
 * @param R red double array
 * @param G green double array
 * @param B blue double array
 * @return output  unsigned char array
 */

unsigned char *rgb_output(double *R, double *G, double *B,
                          unsigned char *output, size_t image_size)
{
    int n;

    for (n=0; n <(int)image_size; n++)
    {
        output[n]=(unsigned char)(R[n]+0.5f);
        output[(int)image_size+n]=(unsigned char)(G[n]+0.5f);
        output[2*(int)image_size+n]=(unsigned char)(B[n]+0.5f);
    }

    return output;
}




/**
 * @brief Computes the gray intensity value of a color image
*
* @f$ gray= (R+ G+ B)/3 \f$
*
* @param data_in input color image
* @param dim size of the image
*
* @return gray output gray image
*
*/

double *gray_intensity(double *gray, double *R, double *G, double *B,
                       size_t image_size)
{
    int index;

    for(index=0; index<(int)image_size; index++)
        gray[index]=(R[index]+G[index]+B[index])/3.;

    return gray;

}






/**
 *  @brief compute the R G B components of the output image from its gray level
 *
 *   Given a color image C=(R, G, B), given its gray level
 *
 * @f$ gray= (R+ G+ B)/3 \f$
 *
 * Given a modified gray image gray1
 *
 * This function computes an output color image C1=(R1,G1,B1) which each channel is proportinal
 * to the input channel and whose gray level is gray1,
 *
 * @f$ R1=\frac{gray1}{gray} R    G1=\frac{gray1}{gray} G    B1= \frac{gray1}{gray} B \f$
 *
 * Note that we make a restriction and does not permit a factor  greater than 3
 *
 *
 * @param R  red channel of the input color image
 * @param G  green channel of the input color image
 * @param B  blue channel of the input color image
 * @param gray gray level of the input color image
 * @param gray1 modified gray image
 * @param dim size of the image
 *
 * @return R_out new red channel
 * @return G_out new green channel
 * @return B_out new blue channel
 */

double *compute_color_from_grayscale(double *Rout, double *Gout, double *Bout,
                                     double *R, double *G, double *B, double *gray,
                                     double *gray1, size_t image_size)
{

    int index;
    double  factor, max;

    for(index=0; index<(int)image_size; index++)
    {
        if(gray[index] <= 1.) gray[index]=1.;
        factor=gray1[index]/ gray[index];
        if( factor > 3.) factor=3.;
        Rout[index]=factor* R[index];
        Gout[index]=factor*G[index];
        Bout[index]=factor*B[index];
        if( Rout[index] > 255. || Gout[index] > 255. || Bout[index] > 255.)
        {
            max=R[index];
            if(G[index] > max) max=G[index];
            if( B[index] > max) max=B[index];
            factor= 255. /max;
            Rout[index]=factor* R[index];
            Gout[index]=factor*G[index];
            Bout[index]=factor*B[index];
        }
    }
    return Rout;
    return Gout;
    return Bout;
}

/**
 * @brief  comparison function
 * given x and y  pointers to doubles.
 * Returns  -1 if x < y
 * 0 if x == y
 * +1 if x > y

*/


int myComparisonFunction(const void *x, const void *y)
{


    double dx, dy;

    dx = *(double *)x;
    dy = *(double *)y;

    if (dx < dy)
    {
        return -1;
    }
    else if (dx > dy)
    {
        return +1;
    }
    return 0;
}

/**
 * @brief Simplest color balance
 *
 * Sort the pixels values.
 * Compute the minimium as the sorted array at position dim x s1/100
 * Compute the maximum as the sorted array at position N(1-s2/100)-1.
 * Saturate the pixels according to the computed minimum and maximum.
 * Affine transformation between [minimum, maximum] and [0,255]
 *
 * @param data input array
 * @param s1 the percentage of saturated pixels on the left
 * @param s2 the percentage of saturated pixels on the right
 * @param dim size of the array
 *
 * @return data_out the scaled array.
 */



double *simplest_color_balance(double *data_out, double *data,
                               size_t image_size, float s1,float s2)
{


    int per1, per2;
    int index;
    double min,max, scale;
    double *sortdata;

    sortdata= (double*) malloc(image_size*sizeof(double));

    memcpy(sortdata, data, image_size*sizeof(double));

    qsort(sortdata, image_size, sizeof sortdata[0], &myComparisonFunction);

    per1 = (int) (s1*image_size/100);
    min=sortdata[per1];

    per2=(int)(s2*image_size/100);
    max= sortdata[image_size-1-per2];


    if(max  <= min)
        for(index=0; index<(int)image_size; index++)
            data_out[index]=max;
    else
    {
        scale=255./(max-min);
        for(index=0; index<(int)image_size; index++)
        {
            if(data[index] < min) data_out[index]=0.;
            else if(data[index]> max) data_out[index]=255.;
            else data_out[index]=scale*(data[index]-min);
        }
    }

    return data_out;



}



