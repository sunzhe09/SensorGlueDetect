// SensorGlueDetect.cpp : �������̨Ӧ�ó������ڵ㡣
//

#include "stdafx.h"
#include<ipp.h>
#include <stdio.h>
#include"mkl.h"
#include<algorithm>
#include<opencv2\opencv.hpp>
#include<vector>
#include<math.h>

#define THRESH_LOW    40.f /* Low threshold for edges detection */
#define THRESH_HIGHT  60.f /* Upper threshold for edges detection */
#define BORDER_VAL 0
#define PI 3.14159265


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

int main()
{
	
	clock_t start, finish;
	double totaltime;
	IppStatus status = ippStsNoErr;

	const char* filename = "C:/Users/bm00133/Desktop/�㽺����ͼ��/20171121ͼƬ/2017-11-21_16_47_35_186.bmp";
	//mat�Ķ�ȡ��ʽ
	cv::Mat input= cv::imread(filename,1);
	//Ipp8u * dst =nullptr;
    //Mat2Ipp8u(input,* &dst);
	

	 //������
	 int *data_x = NULL;
	 int *data_y = NULL;
	 int data_n = 0;
	 
	/*
	IppiSize roi = { input.cols,input.rows };
	
	//����Ӧ��ֵ
	int srcStep = 0;
	Ipp8u* src = ippiMalloc_8u_C1(roi.width, roi.height, &srcStep);
	ippiCopy_8u_C1R((Ipp8u*)input.data, input.step[0], src, srcStep, roi);

	int threshStep = 0;
	Ipp8u* atsuThresh = ippiMalloc_8u_C1(roi.width, roi.height, &threshStep);
	ippiSet_8u_C1R(0, atsuThresh, threshStep, roi);
	Ipp8u threshValue = 0;
	ippiComputeThreshold_Otsu_8u_C1R(src, srcStep, roi, &threshValue);
	ippiThreshold_LTValGTVal_8u_C1R(src, srcStep, atsuThresh, threshStep, roi, threshValue, 0, threshValue, 255);
	*/
	//IplImage�Ķ�ȡ��ʽ
	IplImage *pImage=cvLoadImage(filename, 0);
	int Src_StepBytes = 0;
	
	start = clock();

	IppiRect ROI = {500,1000 ,1600,600};
	int relative = 15;
	IppiSize roi_size = IppiSize();
	roi_size.height = ROI.height;
	roi_size.width = ROI.width;

	Ipp8u *pSrcImage = NULL;
	//pSrcImage = ippiMalloc_8u_C1(pImage->width, pImage->height, &Src_StepBytes);
	//ippiCopy_8u_C1R((Ipp8u*)pImage->imageData, pImage->widthStep, pSrcImage, Src_StepBytes, roi_size);

	pSrcImage = ippiMalloc_8u_C1(ROI.width, ROI.height, &Src_StepBytes);
	ippiCopy_8u_C1R((Ipp8u*)pImage->imageData+ROI.y*pImage->widthStep+ROI.x, pImage->widthStep, pSrcImage, Src_StepBytes, roi_size);

	LigitEnhance(pSrcImage, Src_StepBytes,roi_size );


	/*IplImage *show1 = cvCreateImageHeader(cvSize(roi_size.width, roi_size.height), 8u, 1);
	cvSetData(show1, (uchar*)pSrcImage, Src_StepBytes);*/
	
	
	/*����ֵ��*/


	Ipp8u *pDstImage = NULL;
	int Dst_StepBytes = 0;
	Ipp8u threshold = 0;
	pDstImage = ippiMalloc_8u_C1(ROI.width, ROI.height, &Dst_StepBytes);
	ippiComputeThreshold_Otsu_8u_C1R(pSrcImage, Src_StepBytes, roi_size,&threshold);
	ippiThreshold_LTValGTVal_8u_C1R(pSrcImage, Src_StepBytes, pDstImage, Dst_StepBytes, roi_size, threshold, 0, threshold, 255);


	//IplImage *show = cvCreateImageHeader(cvSize(roi_size.width, roi_size.height),8u,1);
	//cvSetData(show,(uchar*)pDstImage,Src_StepBytes);

	int flag = RegionsSegement(pDstImage,Dst_StepBytes,roi_size,filename);



	/*******��ȡ��Ե*******/  

	Ipp8u *pBuffer = NULL;


	int srcStep = 0, dstStep = 0; int iTmpBufSize = 0;
	IppiDifferentialKernel filterType = ippFilterScharr;
	IppiMaskSize mask = ippMskSize3x3;
	IppiBorderType bodertype = ippBorderRepl;

	status = ippiCannyBorderGetSize(roi_size, filterType, ippMskSize3x3, ipp8u, &iTmpBufSize);
	pBuffer = ippsMalloc_8u(iTmpBufSize);
	status = ippiCannyBorder_8u_C1R(pDstImage, Dst_StepBytes, pDstImage, Dst_StepBytes, roi_size,filterType, ippMskSize3x3, ippBorderRepl, BORDER_VAL, THRESH_LOW, THRESH_HIGHT, ippNormL2, pBuffer);
	
	IplImage *show = cvCreateImageHeader(cvSize(roi_size.width, roi_size.height), 8u, 1);
	cvSetData(show, (uchar*)pDstImage, Dst_StepBytes);


	//Ϊָ������ڴ�
	data_x = (int*)malloc(sizeof(int)*roi_size.width*roi_size.height);
	data_y = (int*)malloc(sizeof(int)*roi_size.width*roi_size.height);
	int tempMin = roi_size.width / 2, tempMax = 0;
	int lefttop_index = 0, righttop_index = 0;

	for (int i = 0; i<roi_size.height; i++)
	{
		for (int j = 0; j < roi_size.width; j++)
		{
			if (pDstImage[i*Dst_StepBytes + j] != 0)
			{
					
				data_x[data_n] = j;
				data_y[data_n] = i;
				data_n++;
				if (max(tempMax, j) > tempMax)
				{
					tempMax = j;
					righttop_index = data_n - 1;
				}

				if (min(tempMin, j)<tempMin)
				{
					tempMin = j;
					lefttop_index = data_n - 1;
				}

			}
		}

	}


	int maxdist = 0;
	float  row = 0,col = 0, row1 = 0, col1 = 0;

	//status= ippsMinMaxIndx_32f(data_x,data_n,col,&indexmin_x,col1,&indexmax_x);

	for (int i = 0; i < data_n; ++i)
	{
		for (int j = i + 1; j < data_n; ++j)
		{
			int distance = (data_x[i] - data_x[j])*(data_x[i] - data_x[j]) +(data_y[i]-data_y[j])*(data_y[i] - data_y[j]);
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
	

	finish = clock();
	totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
	cout << "\n�˳��������ʱ��Ϊ" << totaltime << "�룡" << endl;
	

	//for (int i = 0; i < boderPoints.size(); ++i)
	//{
	//	for (int j = i + 1; j < boderPoints.size(); ++j)
	//	{	
	//	   distance = abs(boderPoints.at(i).x - boderPoints.at(j).x) +abs(boderPoints.at(i).y - boderPoints.at(j).y);

	//		if (distance > max)
	//		{
	//			max = distance;
	//			row = boderPoints.at(i).y;
	//			col = boderPoints.at(i).x;
	//			row1 = boderPoints.at(j).y;
	//			col1 = boderPoints.at(j).x;
	//		}

	//
	//	}
	//}


	_point  pt1;
	pt1.x = col;
	pt1.y = row;
	_point  pt2;
	pt2.x = col1;
	pt2.y = row1;

	cv::line(input, cv::Point(pt1.x+ROI.x,pt1.y+ROI.y), cv::Point(pt2.x+ROI.x,pt2.y+ROI.y), cv::Scalar(0,0,255));
	
	/*Ipp8u2IplimageShow(pDst, dstStep, roi_size);
	LineFitLeastSquares(data_x, data_y,  data_n, res);
	double rho = res[1] / (sqrt(1 + res[0]*res[0]));
	double thea = atan(res[0])+PI/2 ;
	drawLine(input, thea, rho, cv::Scalar(255,0,0));
	drawLine(input, thea, rho+675, cv::Scalar(0,0,255));*/

	ippsFree(pBuffer);
	
	free(data_x);
	free(data_y);
	printf("Exit status %d (%s)\n", (int)status, ippGetStatusString(status));

	return (int)status;


}


cv::Mat SetRegionsColor(const char *filename,int markersNum,IppiSize roiSize,Ipp16u *pMarker,int markerStep)
{
	cv::Mat res = cv::imread(filename);


	//����ͨ������ɫ
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
				//��ʾ���ֵĲ�����������opencv�����ݽṹ����������ȥ��
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

	/*const char *filename = "C:/Users/bm00133/Desktop/�㽺����ͼ��/20171114���ιⷽ��/2017-11-14_20_53_56_676.bmp";
	IplImage *src = cvLoadImage(filename, 0);*/

	IppStatus status = ippStsNoErr;
	IppiMorphAdvState *pState = NULL;
	
	int minLabel = 1;                      
	int maxLabel = 2000;              
	int markersNum = 0;                   
	Ipp8u* pBuffer = NULL;                
	int bufferSize = 0;


	/*Ipp8u *pBinImage = NULL;
	int Bin_StepBytes = 0;*/
	Ipp8u threshold = 0;

	//pBinImage = ippiMalloc_8u_C1(roiSize.width, roiSize.height, &Bin_StepBytes);
	ippiComputeThreshold_Otsu_8u_C1R(pSrcMarker, srcMarkerStep, roiSize, &threshold);
	ippiThreshold_LTValGTVal_8u_C1R(pSrcMarker, srcMarkerStep, pSrcMarker, srcMarkerStep, roiSize, threshold, 0, threshold, 255);


	/*IplImage *show = cvCreateImageHeader(cvSize(roiSize.width, roiSize.height), 8u, 1);
	cvSetData(show, (unsigned char*)pBinImage, Bin_StepBytes);*/
	//IppiBorderType
	//��ͼ����и�ʴ���ʹ���

	MorphologyEroAndDiade(pSrcMarker, srcMarkerStep, roiSize.width, roiSize.height, pSrcMarker);


	/*��ʾ�м���*/
	//IplImage *show = cvCreateImageHeader(cvSize(roiSize.width, roiSize.height), 8u, 1);
	//cvSetData(show, (unsigned char*)pBinImage, Bin_StepBytes);

	//���������Ӧ���Ƕ�ֵ��ͼ��

	unsigned short *pMarker = NULL;//16u��Ӧ�ľ��������ֽڵ�unsigned short 
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

			if (label == 0)//˳��ɸѡ��row���겻�Ե�����
			{

				continue;
			}
			else
			{

				pixelNum[label - 1] += 1;	//��ͨ����������	
								
			}

		}
	}


    //	cv::Mat result=SetRegionsColor(filename, markersNum, roiSize, pMarker, markerStep);

	unsigned char *DeleteIndex = ippsMalloc_8u(markersNum);
	ippsSet_8u(0, DeleteIndex, markersNum);

	for (int i = 0; i < markersNum; i++)
	{
		//���ɸѡ��ȥ���ϴ���߽�С����
		if (pixelNum[i] <= 25000 )
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
				//ɾ������Ҫ������
				pSrcMarker[col + row*(srcMarkerStep / sizeof(unsigned char))] = 0;
				continue;
			}
			//������Ҫ�������ֵ
			pSrcMarker[col + row*(srcMarkerStep / sizeof(unsigned char))] = 255;
		}
	}

  
	

	
	ippsFree(pBuffer);
	ippFree(pixelNum);
	ippFree(pMarker);
	ippFree(pSrcMarker);
	//ippFree(pBinImage);
	//cvReleaseImageHeader(&show);
	ippsFree(DeleteIndex);

	return  (int)status;



}


void LigitEnhance(Ipp8u * src,int srcStep,IppiSize roiSize)
{

	Ipp8u mulFacorValue = 15;
	
	ippiMulC_8u_C1RSfs(src,srcStep, mulFacorValue,src,srcStep,roiSize,0);

	//IplImage *show = cvCreateImageHeader(cvSize(roiSize.width, roiSize.height), 8u, 0);
	//cvSetData(show, (uchar*)src, srcStep);
	

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

	//��ʴ
	status = ippiErodeBorder_8u_C1R(src,srcStep,Sub,srcStep, roiSize, borderType,0,pMorphSpec,pBuffer);
	//����
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
	IppiSize roi = { input.cols,input.rows };
	int dstStep = 0;
	dst = ippiMalloc_8u_C1(input.cols, input.rows, &dstStep);
	ippiCopy_8u_C1R((Ipp8u*)input.data, input.step[0], dst, dstStep, roi);
	

}

void Ipp8u2IplimageShow(Ipp8u * &src, int Src_StepBytes, IppiSize roi_size)
{

	//int Src_StepBytes = 0;
	/*Ipp8u *pSrcImage = NULL;
	pSrcImage = ippiMalloc_8u_C1(roi_size.width, roi_size.height, &Src_StepBytes);
	ippiCopy_8u_C1R((Ipp8u*)pImage->imageData, pImage->widthStep, pSrcImage, Src_StepBytes, roi_size);*/
	
	IplImage *show = cvCreateImageHeader(cvSize(roi_size.width, roi_size.height), IPL_DEPTH_8U, 1);
	cvSetData(show, (uchar*)src, Src_StepBytes);
	cvNamedWindow("result",0);
	cvShowImage("result",show);
	cvSaveImage("C:/Users/bm00133/Desktop/�㽺����ͼ��/canny.bmp",show);
	cvWaitKey(0);
	cvReleaseImageHeader(&show);
	cvDestroyWindow("result");

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
��С���˷����ֱ�ߣ�y = a*x + b; n������; r-���ϵ��[-1,1],fabs(r)->1,˵��x,y֮�����Թ�ϵ�ã�fabs(r)->0��x,y֮�������Թ�ϵ�����������
a = (n*C - B*D) / (n*A - B*B)
b = (A*D - B*C) / (n*A - B*B)
r = E / F
���У�
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

	// ����б��a�ͽؾ�b
	float a, b, temp = 0;
	if (temp = (data_n*A - B*B))// �жϷ�ĸ��Ϊ0
	{
		a = (data_n*C - B*D) / temp;
		b = (A*D - B*C) / temp;
	}
	else
	{
		a = 1;
		b = 0;
	}

	// �������ϵ��r
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







