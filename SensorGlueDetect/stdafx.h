// stdafx.h : ��׼ϵͳ�����ļ��İ����ļ���
// ���Ǿ���ʹ�õ��������ĵ�
// �ض�����Ŀ�İ����ļ�
//

#pragma once

#include "targetver.h"
#include <stdio.h>
#include <tchar.h>
#include<ipp.h>
#include"mkl.h"
#include<opencv2\opencv.hpp>



// TODO:  �ڴ˴����ó�����Ҫ������ͷ�ļ�
void Mat2Ipp8u(cv::Mat input, Ipp8u * &dst);//�ڲ��������ڴ棬�ǵ��ͷ�
void Ipp8u2IplimageShow(Ipp8u * &src, int Src_StepBytes, IppiSize roi_size);

