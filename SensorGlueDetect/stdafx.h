// stdafx.h : 标准系统包含文件的包含文件，
// 或是经常使用但不常更改的
// 特定于项目的包含文件
//

#pragma once
#include "targetver.h"
#include<vector>
#include"ViType.h"
#include <stdio.h>
#include <tchar.h>
#include<ipp.h>
#include"mkl.h"
#include<opencv2\opencv.hpp>



typedef struct hline_t
{
	int p;
	int theta;
}hline;

typedef struct _point_t
{
	int x;
	int y;
}_point;

#define THRESH_LOW    40.f /* Low threshold for edges detection */
#define THRESH_HIGHT  80.f /* Upper threshold for edges detection */
#define BORDER_VAL 0


#ifndef max
#define max(a,b)  (((a) > (b)) ? (a) : (b)) 

#endif
#ifndef min
#define min(a,b)  (((a) > (b)) ? (b) : (a))

#endif


// TODO:  在此处引用程序需要的其他头文件
void Mat2Ipp8u(cv::Mat input, Ipp8u * &dst);//内部申请了内存，记得释放
void imgShow(unsigned char* src, int srcStep, int width, int height);
void Ipp8u2IplimageShow(Ipp8u * &src, int Src_StepBytes, IppiSize roi_size);
int VisBAIP_MorphologyConnectComLabel_FloodFill(IMG_UBBUF pBinary, IMG_UBBUF pMarker, IppiConnectedComp *&pRegin, int &pNum);
void maxDistancePoints(Ipp8u*pSrcImage, int Src_StepBytes, IppiSize roi_size, int relative, IppiRect ROI, const char* filename = NULL);
void LineFitLeastSquares(float *data_x, float *data_y, int data_n, std::vector<float> &vResult);
void hough_line(_point *points, int count, int width, int height, hline* line);
void drawLine(cv::Mat &image, double theta, double rho, cv::Scalar color);
int RegionsSegement(unsigned char *src, int srcStep, IppiSize roiSize, const char *filename);
cv::Mat SetRegionsColor(const char *filename, int markersNum, IppiSize roiSize, Ipp16u *pMarker, int markerStep);
int MorphologyEroAndDiade(unsigned char* src, int srcStep, int srcWidth, int srcHeight, unsigned char *  dst);
void LigitEnhance(Ipp8u * src, int srcStep, IppiSize roiSize);
void loadImageAndsetROI(unsigned char* src, int width, int height, int srcStep, int rel, int thre, IppiRect &ROI);
void LightEvaluate(unsigned char *src, int SrcStep, IppiRect ROI, int num);

