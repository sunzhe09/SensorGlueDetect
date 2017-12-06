// stdafx.h : 标准系统包含文件的包含文件，
// 或是经常使用但不常更改的
// 特定于项目的包含文件
//

#pragma once

#include "targetver.h"
#include <stdio.h>
#include <tchar.h>
#include<ipp.h>
#include"mkl.h"
#include<opencv2\opencv.hpp>
#include"ViType.h"



// TODO:  在此处引用程序需要的其他头文件
void Mat2Ipp8u(cv::Mat input, Ipp8u * &dst);//内部申请了内存，记得释放
void imgShow(unsigned char src, int srcStep, int width, int height);
void Ipp8u2IplimageShow(Ipp8u * &src, int Src_StepBytes, IppiSize roi_size);
int VisBAIP_MorphologyConnectComLabel_FloodFill(IMG_UBBUF pBinary, IMG_UBBUF pMarker, IppiConnectedComp *&pRegin, int &pNum);

