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
#include "time.h"


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define PI2  6.283185307179586  /* 2*pi*/
double *  BoxBlur(double* Src, double *Dest, int Width, int Height, int Radius)
{
    int X, Y, Z, MoveIn, MovOut  ;
    double Sum, InvertAmount = 1.0 / (( 2 * Radius + 1) * (2 * Radius + 1));
    double *Temp = (double *)malloc(Height * Width * sizeof(double));
    double *LinePS, *LinePD;

    /* BoxBlur是一种行列可以分离计算的过程，因此先进行行方向的BoxBlur，*/
    /*  得到中间结果，然后再对中间结果进行列方向的BoxBlur，则可以得到最终结果。*/

    for (Y = 0; Y < Height; Y++)            /* 对每行的数据进行BoxBlur，注意这个过程是行方向无关的，即每行的结算结果不会相互影响，因此和适合于并行化*/
    {
        LinePS = Src + Y * Width;           /* 定位到第Y行第一个像素的内存地址*/
        LinePD = Temp + Y * Width;
        Sum = (Radius + 1) * LinePS[0];     /* 计算行第一个像素左右两边 Radius范围内的累加值，超过边缘的用边缘值代替, 注意这里的加1是为了后续的编码方便，请仔细体味这个+1，很有相反的 */
        for (Z = 0; Z < Radius; Z++)        /* 按理说这里的循环内部是要进行判断Z是否大于等于宽度，否则会出现访问非法内存的错误，但是这种情况只在Radius>Width时出现，这个是没有意义的。*/
            Sum += LinePS[Z];               /* 计算第一个点右侧的像素累加值*/
        for (X = 0; X < Width;X++)          
        {   
            MovOut = X - Radius - 1;                        /* 从左侧移出去的那个像素的坐标*/
            if (MovOut < 0) MovOut = 0;                     /* 用边缘值代替*/
            MoveIn = X + Radius;                            /*从右侧移入的那个像素的坐标*/
            if(MoveIn >= Width) MoveIn = Width - 1;         /*用边缘值代替*/
            Sum = Sum - LinePS[MovOut] + LinePS[MoveIn];    /*  新的累加值 = 旧的累加值 - 移出的值 + 移入的值*/
            LinePD[X] = Sum;                                /* 保存到临时内存中*/
        }
    }

    for (X = 0; X < Width; X++)             /* 接下来对临时数据进行列方向的BoxBlur，得到最终结果，编码方式其实一样的，只是列更改为行*/
    {
        LinePS = Temp + X;                  /* 定位到第X列的第一个像素的内存地址*/
        LinePD = Dest + X;
        Sum = (Radius + 1) * LinePS[0];     /* 以下是同样的道理，无需注释了*/
        for (Z = 0; Z < Radius; Z++)
            Sum += LinePS[Z * Width];
        for (Y = 0; Y < Height; Y++)
        {   
            MovOut = Y - Radius - 1;
            if (MovOut < 0) MovOut = 0;
            MoveIn = Y + Radius;
            if(MoveIn >= Height) MoveIn = Height - 1;
            Sum = Sum - LinePS[MovOut * Width] + LinePS[MoveIn * Width];            /* 注意这里的*Width了*/
            LinePD[Y * Width] = Sum * InvertAmount;                 /* 求平均值，浮点的乘法要比除法快多了，这样为什么不用乘法呢。*/
        }
    }
    free(Temp);         /*  注意释放内存*/
    return Dest;
}
/**
 * @Recurisive   Gaussian Blur
 * @param src    the original image
 * @param dst    the blurred image
 * @param width  the width of image
 * @param height the height of image
 * @param sigma  the Coefficient of the Gaussian
 * @param chan   the number of channels
 * @return dst   the output blured gray image
 */
double * IMG_GaussBlur(double* src, double * dst, int width, int height, float sigma, int chan) /*(1)unsigned char*修改为double(2)unsigned char*&修改为double*/
{
	int    i         = 0;
	int    row       = 0;
	int    col       = 0;
	int    pos       = 0;
	int    channel   = 0;
	int    n         = 0;
	int    bufsize   = 0;        
	int size         = 0;
	int rowstride    = 0;
	int itemp0       = 0;
	int itemp1       = 0;    
	float temp       = 0;
	int    channelsize = width*height;

	if (width>height)
	{
		bufsize = width;
	}
	else
	{
		bufsize = height;
	}    

	double* w1    = (double *) malloc (bufsize * sizeof (double));     /*float修改为double*/
	double *w2    = (double *) malloc (bufsize * sizeof (double));     /*float修改为double*/
	double *in    = (double *) malloc (channelsize * sizeof (double));  /*float修改为double*/
	double *out   = (double *) malloc (channelsize * sizeof (double));  /*float修改为double*/

	/****************************计算高斯核***************************/
	double  q= 0;  /*float修改为double*/
	double  q2, q3; /*float修改为double*/  
	double b0;
	double b1;
	double b2;
	double b3;
	double B    = 0;
	int    N    = 3;

	if (sigma >= 2.5)
	{
		q = 0.98711 * sigma - 0.96330;
	}
	else if ((sigma >= 0.5) && (sigma < 2.5))
	{
		q = 3.97156 - 4.14554 * (double) sqrt ((double) 1 - 0.26891 * sigma);
	}
	else
	{
		q = 0.1147705018520355224609375;
	}

	q2 = q * q;
	q3 = q * q2;
	b0 = (1.57825+(2.44413*q)+(1.4281 *q2)+(0.422205*q3));
	b1 = (        (2.44413*q)+(2.85619*q2)+(1.26661 *q3));
	b2 = (                   -((1.4281*q2)+(1.26661 *q3)));
	b3 = (                                 (0.422205*q3));
	B = 1.0-((b1+b2+b3)/b0);

	/*加速方法 减少循环多次/b0*/
	b1 /= b0;
	b2 /= b0;
	b3 /= b0;
	/*计算高斯核结束*/

	/* 处理图像的多个通道*/
	for (channel = 0; channel < chan; channel++)
	{
		/*获取一个通道的所有像素值*/
		for (i = 0, pos = channel; i < channelsize ; i++, pos += chan)
		{
			/* 0-255 => 1-256 */
			in[i] = (double)(src[pos] + 1.0);
		}

		/*纵向处理*/
		for (row=0 ;row < height; row++)
		{
			pos =  row * width;            
			size        = width;
			rowstride    = 1;    
			bufsize        = size;
			size        -= 1;

			temp =  (in + pos)[0]; 
			w1[0] = temp;  /*！！！计算需要相同初值（否则噪声），但是赋值的时候不能用相同的（否则像素重复或缺失）！！！*/
			w1[1] =temp;
			w1[2] =temp;



			for (  n=3; n <= size ; n++)
			{
				w1[n] = (double)(B*(in + pos)[n*rowstride] +    ((b1*w1[n-1] +     b2*w1[n-2] + b3*w1[n-3] )));

			}
			w1[0] =  (in + pos)[0];  /*左边产生噪声（使用不同初始值时产生），右边有3像素黑色不变（out最后3个像素未设定）*/
			w1[1] =(in + pos)[1];
			w1[2] =  (in + pos)[2];

			

			(out + pos)[size]= w2[size]= (double)w1[size];   /*float修改为double*/
			(out + pos)[size-1]=w2[size-1]=(double)(1.0*w1[size-1]) ; /*float修改为double*/
			(out + pos)[size-2]=w2[size-2]= (double)((1.0)*w1[size-2]); /*float修改为double*/
			w2[size]= (double)w1[size-2]; /*float修改为double*/
			w2[size-1]= (double)w1[size-2]; /*float修改为double*/
			w2[size-2]= (double)w1[size-2]; /*float修改为double*/
			for (n = size-3; n>= 0; n--) 
			{
				(out + pos)[n * rowstride] = w2[n] = (double)(B*w1[n] +    ((b1*w2[n+1] +    b2*w2[n+2] + b3*w2[n+3] )));  /*float修改为double*/

			}    

		}    


		/*横向处理*/
		for (col=0; col < width; col++)  /*wbp 在纵向处理的基础上继续横向处理*/
		{                
			size        = height;
			rowstride    = width;    
			bufsize        = size;
			size        -= 1;

			temp  = (out + col)[0*rowstride];  /* wbp 第col列的第一个数据，复制3份，开始前向滤波*/
			w1[0] = temp;
			w1[1] = temp;
			w1[2] = temp;
			for ( n=3; n <= size ; n++)
			{
				w1[n] = (double)(B*(out + col)[n*rowstride] + ((b1*w1[n-1] +    b2*w1[n-2] + b3*w1[n-3] ))); /*float修改为double*/

			}
			w1[0] =  (out + col)[0];
			w1[1] =  (out + col)[rowstride];
			w1[2] =  (out + col)[2*rowstride];


			temp        = w1[size];
			w2[size]    = temp;
			w2[size-1]    = temp;
			w2[size-2]    = temp;
			(in + col)[size * rowstride]=w1[size];
			(in + col)[(size-1) * rowstride]=w1[size-1];
            (in + col)[(size-2) * rowstride]=w1[size-2];

			for (n = size-3; n >= 0; n--)
			{
				(in + col)[n * rowstride] =w2[n]= (double)(B*w1[n] +    ((b1*w2[n+1] +     b2*w2[n+2] + b3*w2[n+3] )));  /*float修改为double*/

			}                
		}
		/*修正偏移的拷贝方法, 但是修正后图像右边及下边会丢失数据？？？// wbp 并且图像左边及上边多出数据*/
                int x,y;
		for(y=0; y<height; y++)
		{
			itemp0 = y*width*chan;
			itemp1 = (y)*width;                                /*+3  数据没丢失，但是是拷贝而不是原数据*/
			for (x=0; x<width; x++)
			{            
				dst[itemp0+x*chan + channel]=in[itemp1+x]-1;    
			}
		}         
	}

	free (w1);
	free (w2);
	free (in);
	free (out);
        return dst;
}
/**
 * @compute the Coefficient of the Gaussian
 * @param winSize the size of Gaussian window
 * @param sigma the param of Gaussian function
 * @return ikern the Coefficient of the Gaussian
 */
int* buildGaussKern(int winSize, float sigma)   /* [1]删除incline[2]修改int sigma为float sigma */
{
    int wincenter, x;
    float   sum = 0.0f;
    wincenter = winSize / 2;
    float *kern = (float*)malloc(winSize*sizeof(float));
    int *ikern = (int*)malloc(winSize*sizeof(int));
    float SQRT_2PI = 2.506628274631f;
    float sigmaMul2PI = 1.0f / (sigma * SQRT_2PI);
    float divSigmaPow2 = 1.0f / (2.0f * sigma * sigma);
    for (x = 0; x < wincenter + 1; x++)
    {
        kern[wincenter - x] = kern[wincenter + x] = exp(-(x * x)* divSigmaPow2) * sigmaMul2PI;
        sum += kern[wincenter - x] + ((x != 0) ? kern[wincenter + x] : 0.0);
    }
    sum = 1.0f / sum;
    for (x = 0; x < winSize; x++)
    {
        kern[x] *= sum;
        ikern[x] = kern[x] * 256.0f;
    }
    free(kern);
    return ikern;
}
/**
 * @The Gaussian blur function of processing the gray image,and the function is simulating the downward function of Convolution
 * @param pixels the input gray image
 * @param pixelsout the output blured gray image
 * @param width  the width of the gray image
 * @param height  the height of the gray image
 * @param sigma  the param of the gaussian function
 * @return pixelsout the output blured gray image
 */
double* GaussBlur1D(double*  pixels,double*  pixelsout, unsigned int  width, unsigned int  height, float sigma)  /*删掉unsigned  int channels,因为是单通道没有用*/
{                                                                                                                             
 /*[1]修改返回值类型void为double类型,函数最后,return pixelsout,这个空间需要开辟和释放*/ 
 /*[2]修改灰度图像数据pixels和模糊后图像pixelsout的数据类型为double类型,因为图像分通道后就是double类型的,后面的tmpBuffer开辟内存需要修改*/ 

    width = 1 * width;  /*3修改为1，因为三个通道变为了1个通道，存储每行数据的宽度变为了原来的1/3.*/
    if ((width % 4) != 0) width += (4 - (width % 4));

    unsigned int  winsize = (1 + (((int)ceil(3 * sigma)) * 2));  /*窗的大小*/
    int *gaussKern = buildGaussKern(winsize, sigma); /*构建高斯核，计算高斯系数*/
    winsize *= 1; /*3改为1，高斯窗的宽度变为原来的1/3*/
    unsigned int  halfsize = winsize / 2;  /*窗的边到中心的距离*/

    double *tmpBuffer = (double*)malloc(width * height* sizeof(double));  /*开辟新的内存存储处理高斯模糊后的数据*/
    int h,w,k;
    for (h = 0; h < height; h++)    /*外层循环，图像的高度*/
    {
        unsigned int  rowWidth = h * width;     /*当前行的宽度为图像的高度乘以每行图像的数据所占的宽度。因为是按行存储的数组。*/

        for (w = 0; w < width; w++) /*w+=channels，可以修改为w++，因为是单通道数据，而不是三通道数据*/
        {
            double rowR = 0.0;  /*存储r分量的数据*/
            int * gaussKernPtr = gaussKern;/*将高斯系数赋值给gaussKernPtr*/
            int whalfsize = w + width - halfsize;
            unsigned int  curPos = rowWidth + w;  /*当前位置*/
            for (k = 1; k < winsize;k++) /* k += channels修改为k++*/
            {
                unsigned int  pos = rowWidth + ((k + whalfsize) % width);
                int fkern = *gaussKernPtr++;
                rowR += (pixels[pos] * fkern);  /*当前像素值乘以高斯系数，rowR这了泛指单通道的当前像素点高斯处理后的数 */
            }

            tmpBuffer[curPos] = rowR/256; /*除以256,((unsigned char)(rowR >> 8))修改为/256*/

        }
    } 
    halfsize = winsize / 2;
    for (w = 0; w < width; w++)
    {
        for (h = 0; h < height; h++)
        {
            double col_all = 0;
            int hhalfsize = h + height - halfsize;
            for (k = 0; k < winsize; k++)
            {
                col_all += tmpBuffer[((k + hhalfsize) % height)* width + w] * gaussKern[k];
            }
            pixelsout[h * width + w] = col_all/256;    /*(unsigned char)(col_all >> 8)修改为/256*/
        }
    }
    free(tmpBuffer);
    free(gaussKern);
    return pixelsout;
}
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
    double *pas15;
    double *pas80;
    double *pas250;
    clock_t start, finish;      
    double duration;   

    image_size=(int)nx* (int)ny;

    pas=(double*) malloc(image_size*sizeof(double));
    pas15=(double*) malloc(image_size*sizeof(double));
    pas80=(double*) malloc(image_size*sizeof(double));
    pas250=(double*) malloc(image_size*sizeof(double));
    /* initialization of the output*/
    
    for(i=0; i<image_size; i++)
        out[i]=0.;

    /* Compute Retinex output*/
    for(n=0; n<nscales; n++)
    {   
       
     if(scale[0]==15)
       {
         BoxBlur(input, pas,nx,ny,9);
         BoxBlur(pas, pas15,nx,ny,9);
         BoxBlur(pas15, pas,nx,ny,11);
       }
         if(scale[1]==80)
       {
         BoxBlur(input, pas,nx,ny,53);
         BoxBlur(pas, pas80,nx,ny,53);
         BoxBlur(pas80, pas,nx,ny,55);
       }
         if(scale[2]==250)
       {
         BoxBlur(input, pas,nx,ny,167);
         BoxBlur(pas, pas250,nx,ny,167);
         BoxBlur(pas250, pas,nx,ny,167);
       }
       start = clock();       
       /*  GaussBlur1D(input,pas,nx,ny,scale[n]/3);*/ 
       /* convolution(input,scale[n],pas,nx,ny);*/
       /*IMG_GaussBlur(input,pas,nx,ny,scale[n]/3,1);*/
       finish = clock();      
       duration = (double)(finish - start) / CLOCKS_PER_SEC;
       printf( "当前尺度%f\n",scale[n]); 
       printf( "模糊部分耗时%f seconds\n", duration );  
       /* printf("尺度：%f\n",scale[n] );*/
        for(i=0; i<(int)image_size; i++)
            out[i]+=w*(log(input[i])-log(pas[i]));
    }

    free(pas);
    free(pas15);
    free(pas80);
    free(pas250);

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
