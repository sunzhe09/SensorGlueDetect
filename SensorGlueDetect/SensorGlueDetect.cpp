// SensorGlueDetect.cpp : 定义控制台应用程序的入口点。


#include "stdafx.h"
#include<ipp.h>
#include"string.h"
#include <stdio.h>
#include"mkl.h"
#include<algorithm>
#include<opencv2\opencv.hpp>
#include <opencv2/imgproc/imgproc.hpp>  /*该模块包含图像处理函数*/
#include <opencv2/highgui/highgui.hpp>
#include<vector>
#include<math.h>
#include<io.h>
#include<iostream>


#define THRESH_LOW    40.f /* Low threshold for edges detection */
#define THRESH_HIGHT  80.f /* Upper threshold for edges detection */
#define BORDER_VAL 0


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



#ifndef max
#define max(a,b)  (((a) > (b)) ? (a) : (b)) 

#endif
#ifndef min
#define min(a,b)  (((a) > (b)) ? (b) : (a))

#endif

using namespace std;

void LineFitLeastSquares(float *data_x, float *data_y, int data_n, vector<float> &vResult);
void hough_line(_point *points, int count, int width, int height, hline* line);
void drawLine(cv::Mat &image, double theta, double rho, cv::Scalar color);
int RegionsSegement(unsigned char *src, int srcStep, IppiSize roiSize, const char *filename);
cv::Mat SetRegionsColor(const char *filename, int markersNum, IppiSize roiSize, Ipp16u *pMarker, int markerStep);
int MorphologyEroAndDiade(unsigned char* src, int srcStep, int srcWidth, int srcHeight, unsigned char *  dst);
void LigitEnhance(Ipp8u * src, int srcStep, IppiSize roiSize);
void loadImageAndsetROI(unsigned char* src, int width, int height, int srcStep, int rel, int thre, IppiRect &ROI);
void LightEvaluate(unsigned char *src, int SrcStep, IppiRect ROI);

int main()
{
	
	
	//clock_t start, finish;
	//double totaltime;
	IppStatus status = ippStsNoErr;

	//const char* filename = "C:/Users/bm00133/Desktop/点胶侧面图像/20171121图片/2017-11-21_16_47_35_186.bmp";
	const char* filename = "C:/Users/bm00133/Desktop/标准上下半圈/标准下半圈.bmp";


	const char* path1 = "C:/Users/bm00133/Desktop/新建文件夹/暗图像/*.bmp";
	string exten1 = "C:/Users/bm00133/Desktop/新建文件夹/暗图像/";
	vector<string> filenames;

	struct _finddata_t fileinfo;
	long handle;
	handle = _findfirst(path1, &fileinfo);
	if (handle == -1)
	{
		cout << "fail..." << endl;
		
	}
	else
	{
		
		filenames.push_back(exten1+ fileinfo.name);
	}
	while (!_findnext(handle, &fileinfo))
	{
	
		filenames.push_back(exten1 + fileinfo.name);
	}
	_findclose(handle);


		//声明点
	int *data_x = NULL;
	int *data_y = NULL;
	int data_n = 0;

	for (int num = 0; num < filenames.size(); num++)
	{

		//IplImage的读取方式
		IplImage *pImage = cvLoadImage(filenames.at(num).c_str(), 0);

		IppiRect ROI = { pImage->width / 2-900,pImage->height/2-900,1800, 1800};


		LightEvaluate((unsigned char *)pImage->imageData, pImage->widthStep,  ROI);

		//start = clock();

		//范围roi区域的数据的指针和roi矩形参数
		//loadImageAndsetROI((unsigned char*)pImage->imageData, pImage->width, pImage->height, pImage->widthStep, 40, 180, ROI);
		//cvRectangle(pImage, cvPoint(ROI.x, ROI.y), cvPoint(ROI.width + ROI.x - 1, ROI.height + ROI.y - 1), cvScalar(255), 3, 4, 0);


	/*	int Src_StepBytes = 0;


		int relative = 15;
		IppiSize roi_size = IppiSize();
		roi_size.height = ROI.height;
		roi_size.width = ROI.width;

		Ipp8u *pSrcImage = NULL;*/
		//pSrcImage = ippiMalloc_8u_C1(pImage->width, pImage->height, &Src_StepBytes);
		//ippiCopy_8u_C1R((Ipp8u*)pImage->imageData, pImage->widthStep, pSrcImage, Src_StepBytes, roi_size);

		//pSrcImage = ippiMalloc_8u_C1(ROI.width, ROI.height, &Src_StepBytes);
		//ippiCopy_8u_C1R((Ipp8u*)pImage->imageData + ROI.y*pImage->widthStep + ROI.x, pImage->widthStep, pSrcImage, Src_StepBytes, roi_size);

		//IplImage *sw = cvCreateImageHeader(cvSize(ROI.width, ROI.height), 8u, 1);
		//cvSetData(sw, pSrcImage, Src_StepBytes);

		//LigitEnhance(pSrcImage, Src_StepBytes,roi_size );
	
 		//int flag = RegionsSegement(pSrcImage, Src_StepBytes, roi_size, filename);

	
		/*******提取边缘*******/

		//Ipp8u *pBuffer = NULL;

		//int srcStep = 0, dstStep = 0; int iTmpBufSize = 0;
		//IppiDifferentialKernel filterType = ippFilterScharr;
		//IppiMaskSize mask = ippMskSize3x3;
		//IppiBorderType bodertype = ippBorderRepl;

		//status = ippiCannyBorderGetSize(roi_size, filterType, ippMskSize3x3, ipp8u, &iTmpBufSize);
		//pBuffer = ippsMalloc_8u(iTmpBufSize);
		//status = ippiCannyBorder_8u_C1R(pSrcImage, Src_StepBytes, pSrcImage, Src_StepBytes, roi_size, filterType, ippMskSize3x3, ippBorderRepl, BORDER_VAL, THRESH_LOW, THRESH_HIGHT, ippNormL2, pBuffer);

		//IplImage *show = cvCreateImageHeader(cvSize(roi_size.width, roi_size.height), 8u, 1);
		//cvSetData(show, (uchar*)pSrcImage, Src_StepBytes);

		//cvSaveImage("20.bmp",pImage);


		//计算距离最远的两个点并连接直线
		/*
		//为指针分配内存
		data_x = (int*)malloc(sizeof(int)*roi_size.width*roi_size.height);
		data_y = (int*)malloc(sizeof(int)*roi_size.width*roi_size.height);
		int tempMin = roi_size.width / 2, tempMax = 0;
		int lefttop_index = 0, righttop_index = 0;

		for (int i = 0; i < roi_size.height; i++)
		{
			for (int j = 0; j < roi_size.width; j++)
			{
				if (pSrcImage[i*Src_StepBytes + j] != 0)
				{

					data_x[data_n] = j;
					data_y[data_n] = i;
					data_n++;
					if (max(tempMax, j) > tempMax)
					{
						tempMax = j;
						righttop_index = data_n - 1;
					}

					if (min(tempMin, j) < tempMin)
					{
						tempMin = j;
						lefttop_index = data_n - 1;
					}

				}
			}

		}


		int maxdist = 0;
		int row = 0, col = 0, row1 = 0, col1 = 0;
		int k = max(data_y[righttop_index], data_y[lefttop_index]);

		for (int i = 0; i < data_n; ++i)
		{
			for (int j = i + 1; j < data_n; ++j)
			{
				if (data_y[j]<k &&data_y[i]<k&&abs(data_y[i] - data_y[j]) <= 20 && data_x[i] > tempMin + relative&&data_x[j] > tempMin + relative&&data_x[i] < tempMax - relative&&data_x[j] < tempMax - relative)
				{
					int distance = (data_x[i] - data_x[j])*(data_x[i] - data_x[j]) + (data_y[i] - data_y[j])*(data_y[i] - data_y[j]);
					if (distance > maxdist)
					{
						maxdist = distance;
						row = data_y[i];
						col = data_x[i];
						row1 = data_y[j];
						col1 = data_x[j];
					}
				}
			}

		}


		_point  pt1;
		pt1.x = col;
		pt1.y = row;
		_point  pt2;
		pt2.x = col1;
		pt2.y = row1;

		//finish = clock();
		//totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
		//cout << "\n此程序的运行时间为" << totaltime << "秒！" << endl;
		

		cv::Mat input = cv::imread(filename, 1);
		cv::line(input, cv::Point(pt1.x + ROI.x, pt1.y + ROI.y), cv::Point(pt2.x + ROI.x, pt2.y + ROI.y), cv::Scalar(0, 0, 255));
		cv::line(input, cv::Point(data_x[lefttop_index] + ROI.x, data_y[lefttop_index] + ROI.y), cv::Point(data_x[righttop_index] + ROI.x, data_y[righttop_index] + ROI.y), cv::Scalar(0, 255, 255));
	
	*/

	/*Ipp8u2IplimageShow(pDst, dstStep, roi_size);
	LineFitLeastSquares(data_x, data_y,  data_n, res);
	double rho = res[1] / (sqrt(1 + res[0]*res[0]));
	double thea = atan(res[0])+PI/2 ;
	drawLine(input, thea, rho, cv::Scalar(255,0,0));
	drawLine(input, thea, rho+675, cv::Scalar(0,0,255));*/

	/*ippsFree(pBuffer);*/
  
	//free(data_x);
	//free(data_y);
	//printf("Exit status %d (%s)\n", (int)status, ippGetStatusString(status));

	 }

	cvWaitKey(27);
	return (int)status;
	
	
}


cv::Mat SetRegionsColor(const char *filename,int markersNum,IppiSize roiSize,Ipp16u *pMarker,int markerStep)
{
	cv::Mat res = cv::imread(filename);


	//给连通区域上色
	cv::RNG rng(time(0));
	std::vector<cv::Scalar> color;


	for (int i = 0; i < markersNum; ++i)
	{
		cv::Scalar temp = cv::Scalar(rng.uniform(0, 255), rng.uniform(0, 255), rng.uniform(0, 255));
		color.push_back(temp);
	}

	for (int row = 0; row < roiSize.height;++row)
	{
		for (int col = 0; col < roiSize.width; ++col)
		{
			int label = (int)pMarker[row*markerStep / sizeof(unsigned short) + col];

			if (label != 0)
			{
				//显示部分的操作，采用了opencv的数据结构，后续可以去掉
				res.at<cv::Vec3b>(row, col)[0] = color.at(label - 1).val[0];
				res.at<cv::Vec3b>(row, col)[1] = color.at(label - 1).val[1];
				res.at<cv::Vec3b>(row, col)[2] = color.at(label - 1).val[2];
			}
		}
	}

	return res;

}

int RegionsSegement( unsigned char *pSrcMarker,int srcMarkerStep,IppiSize roiSize,const char *filename)
{

	/*const char *filename = "C:/Users/bm00133/Desktop/点胶侧面图像/20171114条形光方案/2017-11-14_20_53_56_676.bmp";
	IplImage *src = cvLoadImage(filename, 0);*/

	IppStatus status = ippStsNoErr;
	IppiMorphAdvState *pState = NULL;
	
	int minLabel = 1;                      
	int maxLabel = 20000;              
	int markersNum = 0;                   
	Ipp8u* pBuffer = NULL;                
	int bufferSize = 0;

	

	/*Ipp8u *pBinImage = NULL;
	int Bin_StepBytes = 0;*/
	Ipp8u threshold = 0;
	IppCmpOp ippcmpop= ippCmpLessEq;
	//pBinImage = ippiMalloc_8u_C1(roiSize.width, roiSize.height, &Bin_StepBytes);
	ippiComputeThreshold_Otsu_8u_C1R(pSrcMarker, srcMarkerStep, roiSize, &threshold);
	ippiCompareC_8u_C1R(pSrcMarker, srcMarkerStep, threshold, pSrcMarker, srcMarkerStep, roiSize, ippcmpop);
	//ippiThreshold_LTValGTVal_8u_C1R(pSrcMarker, srcMarkerStep, pSrcMarker, srcMarkerStep, roiSize, 10, 255, 10, 0);


	/*显示中间结果*/
	IplImage *show = cvCreateImageHeader(cvSize(roiSize.width, roiSize.height), 8u, 1);
	cvSetData(show, (unsigned char*)pSrcMarker, srcMarkerStep);

	//IppiBorderType
	//对图像进行腐蚀膨胀处理

	MorphologyEroAndDiade(pSrcMarker, srcMarkerStep, roiSize.width, roiSize.height, pSrcMarker);




	//传输进来的应该是二值化图像

	unsigned short *pMarker = NULL;//16u对应的就是两个字节的unsigned short 
	int markerStep = 0;
	pMarker = ippiMalloc_16u_C1(roiSize.width, roiSize.height, &markerStep);

	ippiConvert_8u16u_C1R(pSrcMarker, srcMarkerStep, pMarker, markerStep, roiSize);

	ippiLabelMarkersGetBufferSize_16u_C1R(roiSize, &bufferSize);

	pBuffer = ippsMalloc_8u(bufferSize);//

	ippiLabelMarkers_16u_C1IR(pMarker, markerStep, roiSize, minLabel, maxLabel, ippiNormL1, &markersNum, pBuffer);

	unsigned int *pixelNum = ippsMalloc_32u(markersNum);
	ippsSet_32s(0, (signed int*)pixelNum, markersNum);


	for (int row = 0; row < roiSize.height; ++row)
	{
		for (int col = 0; col < roiSize.width; ++col)
		{
			int  label = (int)pMarker[row*markerStep / sizeof(unsigned short) + col];

			if (label == 0)//顺便筛选掉row坐标不对的区域
			{

				continue;
			}
			else
			{

				pixelNum[label - 1] += 1;	//连通域的像素面积	
								
			}

		}
	}

	
    //	cv::Mat result=SetRegionsColor(filename, markersNum, roiSize, pMarker, markerStep);

	unsigned char *DeleteIndex = ippsMalloc_8u(markersNum);
	ippsSet_8u(0, DeleteIndex, markersNum);


	for (int i = 0; i < markersNum; i++)
	{
		//面积筛选，去除较大或者较小区域
		if (pixelNum[i] <= 15000 )
		{
			DeleteIndex[i] = 1;

		}
	}

	for (int row = 0; row < roiSize.height; row++)
	{
		for (int col = 0; col < roiSize.width; col++)
		{
			int Label = (int)pMarker[col + row*(markerStep / sizeof(unsigned short))];
			if (Label == 0 || 1 == DeleteIndex[Label - 1])
			{
				//删除不合要求区域
				pSrcMarker[col + row*(srcMarkerStep / sizeof(unsigned char))] = 0;
				continue;
			}
			//给满足要求的区域赋值
			pSrcMarker[col + row*(srcMarkerStep / sizeof(unsigned char))] = 255;
		}
	}


	

	
	ippsFree(pBuffer);
	ippFree(pixelNum);
	ippFree(pMarker);
	
	//ippFree(pBinImage);
	//cvReleaseImageHeader(&show);
	ippsFree(DeleteIndex);

	return  status;



}


void LigitEnhance(Ipp8u * src,int srcStep,IppiSize roiSize)
{

	Ipp8u mulFacorValue = 6.5;

	
	ippiMulC_8u_C1RSfs(src,srcStep, mulFacorValue,src,srcStep,roiSize,0);
	

	
	

}

int MorphologyEroAndDiade(unsigned char* src, int srcStep, int srcWidth, int srcHeight, unsigned char *dst)
{
	IppStatus status = ippStsNoErr;

	//Ipp8u pMask[3 * 3] = { 1,1,1,1,0,1,1,1,1 };
	Ipp8u pMask[5 * 5] = { 1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1 };

	//IppiSize maskSize = { 3,3 };
	IppiSize maskSize = { 5,5 };

	int pSpecSize;
	int pBufferSize;

	int dstStep = srcStep;
	IppiSize roiSize = { srcWidth,srcHeight };
	IppiBorderType borderType = ippBorderRepl;
	IppiMorphState *pMorphSpec;
	Ipp8u *pBuffer = NULL;

	Ipp8u *Sub = NULL;
	Sub = ippiMalloc_8u_C1(roiSize.width, roiSize.height, &srcStep);
	memset(Sub, 0, roiSize.width*roiSize.height);

	status = ippiMorphologyBorderGetSize_8u_C1R(roiSize, maskSize, &pSpecSize, &pBufferSize);
	if (status != ippStsNoErr)
	{

		return -1;
	}
	pMorphSpec = (IppiMorphState*)ippsMalloc_8u(pSpecSize);
	pBuffer = (Ipp8u*)ippsMalloc_8u(pBufferSize);

	status = ippiMorphologyBorderInit_8u_C1R(roiSize, pMask, maskSize, pMorphSpec, pBuffer);

	if (status != ippStsNoErr)
	{
		ippsFree(pMorphSpec);
		ippsFree(pBuffer);
		return -1;
	}

	//腐蚀
	status = ippiErodeBorder_8u_C1R(src,srcStep,Sub,srcStep, roiSize, borderType,0,pMorphSpec,pBuffer);
	//膨胀
	status = ippiDilateBorder_8u_C1R(Sub, srcStep, dst, srcStep, roiSize, borderType, 0, pMorphSpec, pBuffer);

	//IplImage *show = cvCreateImageHeader(cvSize(roiSize.width, roiSize.height), 8u, 1);
	//cvSetData(show, (unsigned char*)dst, srcStep);

	ippFree(Sub);
	ippFree(pMorphSpec);
	ippFree(pBuffer);

	return (int)status;

}


std::vector<cv::Point> getPoints(cv::Mat &image, int value)
{
	int nl = image.rows; // number of lines
	int nc = image.cols * image.channels();
	std::vector<cv::Point> points;
	for (int j = 0; j < nl; j++)
	{
		uchar* data = image.ptr<uchar>(j);
		for (int i = 0; i < nc; i++)
		{
			if (data[i] == value)
			{
				points.push_back(cv::Point(i, j));
			}
		}
	}
	return points;
}

void drawLine(cv::Mat &image, double theta, double rho, cv::Scalar color)
{
	if (theta < PI / 4. || theta > 3.*PI / 4.)// ~vertical line
	{
		cv::Point pt1(rho / cos(theta), 0);
		cv::Point pt2((rho - image.rows * sin(theta)) / cos(theta), image.rows);
		cv::line(image, pt1, pt2, cv::Scalar(255), 1);
	}
	else
	{
		cv::Point pt1(0, rho / sin(theta));
		cv::Point pt2(image.cols, (rho - image.cols * cos(theta)) / sin(theta));
		cv::line(image, pt1, pt2, color, 1);
	}
}

void Mat2Ipp8u(cv::Mat input,Ipp8u * &dst)
{
	//dst在外部申请
	IppiSize roi = { input.cols,input.rows };
	int dstStep = 0;
	ippiCopy_8u_C1R((Ipp8u*)input.data, input.step[0], dst, dstStep, roi);
	
}



void hough_line(_point *points, int count, int width, int height, hline* line)
{
	int i, j;
	int p;
	int theta;
	int maxCnt = 0;
	int maxP, maxTheta;
	int maxd;
	int mind;
	int** m;
	int* buf;

	maxd = max(width, height);
	mind = maxd;
	maxd = maxd + (maxd >> 1);
	maxd += mind;
	m = (int**)malloc(181 * sizeof(int*));
	if (!m)
		return;
	buf = (int*)malloc(181 * maxd * sizeof(int));
	if (!buf)
	{
		free(m);
		return;
	}
	memset(buf, 0, 181 * maxd * sizeof(int));
	for (i = 0; i < 181; i++)
	{
		m[i] = &buf[i*maxd];
	}

	for (i = 0; i < count; i++)
	{
		for (theta = 0; theta < 181; theta++)
		{
			p = points[i].x * cos(theta * PI / 180.0) + points[i].y * sin(theta * PI / 180.0) + mind;
			m[theta][p]++;
			assert(p < maxd + mind);
		}
	}
	for (i = 0; i < 181; i++) {
		for (j = 0; j < maxd; j++) {
			if (m[i][j] > maxCnt) {
				maxCnt = m[i][j];
				maxTheta = i;
				maxP = j;
			}
		}
	}
	line->p = maxP - mind;
	line->theta = maxTheta;
	free(buf);
	free(m);
}

/*************************************************************************
最小二乘法拟合直线，y = a*x + b; n组数据; r-相关系数[-1,1],fabs(r)->1,说明x,y之间线性关系好，fabs(r)->0，x,y之间无线性关系，拟合无意义
a = (n*C - B*D) / (n*A - B*B)
b = (A*D - B*C) / (n*A - B*B)
r = E / F
其中：
A = sum(Xi * Xi)
B = sum(Xi)
C = sum(Xi * Yi)
D = sum(Yi)
E = sum((Xi - Xmean)*(Yi - Ymean))
F = sqrt(sum((Xi - Xmean)*(Xi - Xmean))) * sqrt(sum((Yi - Ymean)*(Yi - Ymean)))

**************************************************************************/
void LineFitLeastSquares(float *data_x, float *data_y, int data_n, vector<float> &vResult)
{
	float A = 0.0;
	float B = 0.0;
	float C = 0.0;
	float D = 0.0;
	float E = 0.0;
	float F = 0.0;

	for (int i = 0; i<data_n; i++)
	{
		A += data_x[i] * data_x[i];
		B += data_x[i];
		C += data_x[i] * data_y[i];
		D += data_y[i];
	}

	// 计算斜率a和截距b
	float a, b, temp = 0;
	if (temp = (data_n*A - B*B))// 判断分母不为0
	{
		a = (data_n*C - B*D) / temp;
		b = (A*D - B*C) / temp;
	}
	else
	{
		a = 1;
		b = 0;
	}

	// 计算相关系数r
	float Xmean, Ymean;
	Xmean = B / data_n;
	Ymean = D / data_n;

	float tempSumXX = 0.0, tempSumYY = 0.0;
	for (int i = 0; i<data_n; i++)
	{
		tempSumXX += (data_x[i] - Xmean) * (data_x[i] - Xmean);
		tempSumYY += (data_y[i] - Ymean) * (data_y[i] - Ymean);
		E += (data_x[i] - Xmean) * (data_y[i] - Ymean);
	}
	F = sqrt(tempSumXX) * sqrt(tempSumYY);

	float r;
	r = E / F;

	vResult.push_back(a);
	vResult.push_back(b);
	vResult.push_back(r*r);
}

/******************************************************
roi通过设置矩形的左上角点的坐标和长宽来选择图像中需要计算的区域
在将计算结果映射到源图像中的时候
只需要在相应的横纵坐标上加上矩形起点的横纵坐标即可
dst的步长用src的步长
*******************************************************/

void loadImageAndsetROI(unsigned char* src, int width, int height,int srcStep, int rel,int thre,IppiRect &ROI)
{


	clock_t start, finish;
	double totaltime;


	start = clock();
	//以有图像中心对应的水平垂直两条直线作为搜索线，搜索上下左右四个边界点
	int hasGlue = 0;//判断是否是无胶的图片，搜索不到胶水的roi则报错，并将原图传给dst

	if (!src)
	{
		cout << "传入源图像失败！" << endl;
		return;
	}

	if (thre < 0 || thre>255)
	{
		cout << "阈值参数非法！" << endl;
		return;
	}


	int top=0, bottom=0, left=0, right=0;
	int flag = 0;

	 
	//top

	for ( int i = 0; i < height; i++)
	{
		
		if (src[i*srcStep + width / 2] > thre&&src[(i+20)*srcStep + width / 2]> thre)//这里判断间隔20个像素是和实际胶水的宽度有关系的参数
		{
			top = i;
			break;
		}
		
	}

	//判断右边三分之1位置
	if (top == 0)
	{
		for(int i = 0; i < height; i++)
		{
			if (src[i*srcStep + 2 * width / 3] > thre&&src[(i + 20)*srcStep + 2 * width / 3]> thre)
			{
				top = i;
				break;
			}
		}
		
	}

	//判断左边三分之一
	if (top == 0)
	{
		for (int i = 0; i < height; i++)
		{
			if (src[i*srcStep +   width / 3] > thre&&src[(i + 20)*srcStep +  width / 3]> thre)
			{
				top = i;
				break;
			}
		}

	}



	if (top == 0)
	{
		cout << "搜索roi-top失败！" << endl;
		hasGlue = 1;
	}
	else
	{
		//扫描直到有一行没有白色的区域
		for (int j = top - 1; j > 0; --j)
		{
			flag = 0;
			while (flag<width)
			{
				for (int k = 0; k<width; k++)
				{

					if (src[j*srcStep + k] < thre)
						flag++;
				}

			}

			if (flag == width)
			{
				top = j;
				break;
			}
		}
	}


	


	//bottom
	for (int i = height - 1; i > 0; --i)
	{
		if (src[i*srcStep + width / 2] > thre&& src[(i-20)*srcStep + width / 2] > thre)
		{
			bottom = i;
			break;

		}
	}

	//判断右边三分之1位置
	if (bottom == 0)
	{
		for (int i = height - 1; i > 0; --i)
		{
			if (src[i*srcStep + 2*width / 3] > thre&& src[(i - 20)*srcStep + 2*width / 3] > thre)
			{
				bottom = i;
				break;

			}
		}

	}

	//判断左边三分之一
	if (bottom == 0)
	{
		for (int i = height - 1; i > 0; --i)
		{
			if (src[i*srcStep +  width / 3] > thre&& src[(i - 20)*srcStep +  width / 3] > thre)
			{
				bottom = i;
				break;

			}
		}

	}


	if (bottom == 0)
	{
		cout << "搜索roi-bottom失败！" << endl;
		hasGlue = 1;
	}
	else
	{
		for (int j = bottom + 1; j <height - 1; ++j)
		{
			flag = 0;
			while (flag<width)
			{
				for (int k = 0; k<width; k++)
				{
					if (src[j*srcStep + k] < thre)
						flag++;
				}

			}

			if (flag == width)
			{
				bottom = j;
				break;
			}
		}
	}

	

	//left
	//寻找第一个参考点
	for (int i = 0; i < width; ++i)
	{
		if (src[height / 2 * srcStep + i] > thre&&src[height / 2 * srcStep + i+20] > thre)
		{
			left = i;
			break;
		}
	}

	//上面三分之一
	if (left == 0)
	{
		for (int i = 0; i < width; ++i)
		{
			if (src[height / 3 * srcStep + i] > thre&&src[height / 3 * srcStep + i + 20] > thre)
			{
				left = i;
				break;
			}
		}
	}

	//下面三分之一
	if (left == 0)
	{
		for (int i = 0; i < width; ++i)
		{
			if (src[2*height / 3 * srcStep + i] > thre&&src[2*height / 3 * srcStep + i + 20] > thre)
			{
				left = i;
				break;
			}
		}
	}


	if (left== 0)
	{
		cout << "搜索roi-left失败！" << endl;
		hasGlue = 1;
	}
	else
	{
		//以此边缘点为起点像外边搜索
		for (int j = left - 1; j > 0; j--)
		{
			flag = 0;
			while (flag < height)
			{
				for (int k = 0; k < height; k++)
				{
					if (src[k*srcStep + j] < thre)
						flag++;
				}
			}


			if (flag == height)
			{
				left = j;
				break;
			}
		}
	}

	


	//right
	for (int i = width - 1; i > 0; --i)
	{
		if (src[height / 2 * srcStep + i] > thre&&src[height / 2 * srcStep + i-20] > thre)
		{
			right = i;
			break;
		}
	}

	if (right == 0)
	{
		for (int i = width - 1; i > 0; --i)
		{
			if (src[height / 3 * srcStep + i] > thre&&src[height / 3 * srcStep + i - 20] > thre)
			{
				right = i;
				break;
			}
		}
	}

	if (right == 0)
	{
		for (int i = width - 1; i > 0; --i)
		{
			if (src[2*height / 3 * srcStep + i] > thre&&src[2*height / 3 * srcStep + i - 20] > thre)
			{
				right = i;
				break;
			}
		}
	}



	if (right == 0)
	{
		cout << "搜索roi-right失败！" << endl;
		hasGlue = 1;
	}
	else
	{
		for (int j = right + 1; j < width; j++)
		{
			flag = 0;
			while (flag < height)
			{
				for (int k = 0; k < height; k++)
				{
					if (src[k*srcStep + j] < thre)
						flag++;
				}
			}

			if (flag == height)
			{
				right = j;
				break;
			}
		}
	}

	
	//没有胶水的时候赋值原图并返回
	if (hasGlue == 1)
	{
		ROI.x = 0;
		ROI.y = 0;
		ROI.width = width;
		ROI.height =height;

		
		return;

	}



	//roi向外扩充40个像素
	top -=rel;
	left -= rel;
	right += rel;
	bottom += rel;

	if (top < 0 || left < 0 || right < 0 || bottom < 0)
	{
		cout << "扩充距离偏大，需重新设置"<<endl;
		return ;
	}

	ROI.x = left;
	ROI.y = top;
	ROI.width = right - left+1;
	ROI.height = bottom - top + 1;

	finish = clock();
	totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
	cout << "\n此程序的运行时间为" << totaltime << "秒！" << endl;
	


}


/**********************************************/
// VisBAIP_MorphologyConnectCom_FloodFill
//		FloodFill找连通域，只标记一个连通域
// Input:
// 	 pBinary  输入图像（二值图）
//	 seed  种子点
//newVal  赋的新值
// Output:
//   pRegin   连通域信息
// 	 pMarker  连通域被标记的图像
// Return:
//    int  0-通过 -1-失败 其他-待定
// Author: Jiang He, 11/10/2017
/**********************************************/
int VisBAIP_MorphologyConnectCom_FloodFill(Ipp8u*  pBinary, int Width, int Height, IppiPoint seed, int newVal, IppiConnectedComp &pRegin)
{
	IppStatus status = ippStsNoErr;

	int pBufferSize = 0;
	status = ippiFloodFillGetSize({ Width,Height }, &pBufferSize);
	if (status != ippStsNoErr) return -1;

	int imgStep = Width * sizeof(Ipp8u);
	Ipp8u* pBuffer = (Ipp8u*)malloc(pBufferSize * sizeof(Ipp8u));
	status = ippiFloodFill_8Con_8u_C1IR(pBinary, imgStep, { Width,Height }, seed, newVal, &pRegin, pBuffer);
	free(pBuffer);
	if (status != ippStsNoErr) return -1;

	return 0;
}

/**********************************************/
// VisBAIP_MorphologyConnectComLabel_FloodFill
//		标记所有连通域
// Input:
// 	 pBinary  输入图像（二值图，前景为255，背景为0）
// Output:
//   pRegin   连通域信息
// 	 pMarker  连通域被标记的图像(1~254)
//   pNum     连通域个数
// Return:
//    int  0-通过 -1-失败 其他-待定
// Author: Jiang He, 11/10/2017
/**********************************************/
int VisBAIP_MorphologyConnectComLabel_FloodFill(IMG_UBBUF pBinary, IMG_UBBUF pMarker, IppiConnectedComp *&pRegin, int &pNum)
{


	int status = 0;

	int Width = pBinary.size.width;
	int Height = pBinary.size.height;
	//
	if (pBinary.ptr == NULL || Height < 5 || Width < 5)  return -1;
	if (pMarker.ptr == NULL || pMarker.size.height != pBinary.size.height || pMarker.size.width != pBinary.size.width) return -1;

	memcpy(pMarker.ptr, pBinary.ptr, Width*Height * sizeof(unsigned char));

	Ipp64f pSum = 0;;
	status = ippiSum_8u_C1R(pMarker.ptr, Width * sizeof(Ipp8u), { Width,Height }, &pSum);
	if ((int)pSum % 255 || status != ippStsNoErr)
	{
		return -1;
	}

	//计算前景的像素个数
	int PixNum = (int)(pSum / 255);
	int nowPix = PixNum;

	IMG_INT CCNum = 0;
	int newVal = 1;//从1开始标记，1~254
	vector<IppiConnectedComp> Regin;
	Regin.clear();
	while (nowPix > 0)
	{
		IppiPoint seed = { -1,-1 };
		int flag = 0;
		for (int i = 0; i < Height; i++)
		{
			for (int j = 0; j < Width; j++)
			{
				if (pMarker.ptr[j + i*Width] == 255)
				{
					seed.x = j;
					seed.y = i;
					flag = 1;
					break;
				}
			}
			if (flag)
				break;
		}

		if (seed.x == -1) break;

		IppiConnectedComp ppRegin;
		status = VisBAIP_MorphologyConnectCom_FloodFill(pMarker.ptr, Width, Height, seed, newVal, ppRegin);
		if (status != ippStsNoErr)
		{
			return -1;
		}

		CCNum++;
		Regin.push_back(ppRegin);
		newVal++;
		if (newVal == 255) newVal = 1;
		nowPix = nowPix - ppRegin.area;
	}

	pNum = CCNum;

	pRegin = (IppiConnectedComp*)malloc(CCNum * sizeof(IppiConnectedComp));

	for (int i = 0; i < Regin.size(); i++)
	{
		
	
		if(Regin[i].area>30000)
		{
			pRegin[i] = Regin[i];
			for (int n = Regin[i].rect.y; n < Regin[i].rect.height; n++)
			{
				for (int m = Regin[i].rect.x; m < Regin[i].rect.width; m++)
				{
					
					pMarker.ptr[n*pMarker.linestep + m] == Regin[i].value[0] ?  255 : 0;
						
				}
			}
			
		}

		else
		{
			for (int n = Regin[i].rect.y; n < Regin[i].rect.height; n++)
			{
				for (int m = Regin[i].rect.x; m < Regin[i].rect.width; m++)
				{
					pMarker.ptr[n*pMarker.linestep + m] == Regin[i].value[0] ? 0 : 255;
				}
			}
		}
		
	}

	


	return 0;
}

void LightEvaluate(unsigned char *src,int SrcStep,IppiRect ROI )
{

	Ipp8u* pSrc = NULL;
	int pSrc_step = 0;

	//扣取矩形roi里面的像素
	pSrc = ippiMalloc_8u_C1(ROI.width, ROI.height, &pSrc_step);
	ippiCopy_8u_C1R((Ipp8u*)src + ROI.y*SrcStep + ROI.x, SrcStep, pSrc, pSrc_step, IppiSize{ ROI.width,ROI.height });

	LigitEnhance(pSrc, pSrc_step, IppiSize{ ROI.width,ROI.height });

	Ipp8u threshold = 0;
	IppCmpOp ippcmpop = ippCmpLessEq;

	ippiComputeThreshold_Otsu_8u_C1R(pSrc, pSrc_step, IppiSize{ ROI.width,ROI.height }, &threshold);
	ippiCompareC_8u_C1R(pSrc, pSrc_step, threshold, pSrc, pSrc_step, IppiSize{ ROI.width,ROI.height }, ippcmpop);

	IplImage *show = cvCreateImageHeader(cvSize(ROI.width, ROI.height), 8u, 1);
	cvSetData(show, (uchar*)pSrc, pSrc_step);

	//ippiThreshold_LTValGTVal_8u_C1R(pSrc, pSrc_step, pSrc, pSrc_step, IppiSize{ ROI.width,ROI.height }, 10, 255, 10, 0);


	//连通域分析
	unsigned short *pMarker = NULL; 
	int minLabel = 1;
	int maxLabel = 20000;
	int markersNum = 0;
	Ipp8u* pBuffer = NULL;
	int bufferSize = 0;

	int markerStep = 0;
	pMarker = ippiMalloc_16u_C1(ROI.width, ROI.height, &markerStep);

	ippiConvert_8u16u_C1R(pSrc, pSrc_step, pMarker, markerStep, IppiSize{ ROI.width,ROI.height });

	//计算内存的空间
	ippiLabelMarkersGetBufferSize_16u_C1R(IppiSize{ ROI.width,ROI.height }, &bufferSize);
	//申请指针空间
	pBuffer = ippsMalloc_8u(bufferSize);
	//给连通域标记
	ippiLabelMarkers_16u_C1IR(pMarker, markerStep, IppiSize{ ROI.width,ROI.height }, minLabel, maxLabel, ippiNormL1, &markersNum, pBuffer);

	unsigned int *pixelNum = ippsMalloc_32u(markersNum);
	ippsSet_32s(0, (signed int*)pixelNum, markersNum);



	for (int row = 0; row < ROI.height; ++row)
	{
		for (int col = 0; col < ROI.width; ++col)
		{
			int  label = (int)pMarker[row*markerStep / sizeof(unsigned short) + col];

			if (label == 0)//顺便筛选掉row坐标不对的区域
			{

				continue;
			}
			else
			{

				pixelNum[label - 1] += 1;	//连通域的像素面积	

			}

		}
	}




	unsigned char *DeleteIndex = ippsMalloc_8u(markersNum);
	ippsSet_8u(0, DeleteIndex, markersNum);
	double Area = 0.0;


	for (int i = 0; i < markersNum; i++)
	{
		//面积筛选，去除较大或者较小区域
		if (pixelNum[i] <= 100000|| pixelNum[i]>400000)
		{
			DeleteIndex[i] = 1;

		}
	}

	//删除不合要求区域
	int Lbl = 0,n=0;
	for (int row = 0; row < ROI.height; row++)
	{
		for (int col = 0; col < ROI.width; col++)
		{
			int	 label = (int)pMarker[col + row*(markerStep / sizeof(unsigned short))];
			if (label == 0 || 1 == DeleteIndex[label - 1])
			{	
				src[col + ROI.x + (row + ROI.y)*SrcStep] = 0;
				continue;
			}
			else
			{
				//这里理论上只有一个region，就是内圆中心的部分
				Area += src[col+ROI.x + (row+ROI.y)*SrcStep ];
				src[col+ ROI.x + (row + ROI.y)*SrcStep ] = 255;
				Lbl = label;
				n++;
			}
		
			
			
		}
	}
	

	

	Area /=n;

	cout << "灰度均值为："<<Area << endl;


	

}

void imgShow(unsigned char src,int srcStep,int width,int height)
{
	IplImage *show = cvCreateImageHeader(cvSize(width, height), 8u, 0);
	cvSetData(show, (uchar*)src, srcStep);

	cvWaitKey(0);

	cvReleaseImageHeader(&show);
}




