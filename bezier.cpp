
#include "opencv2/core.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/highgui.hpp"
#include <stdio.h>
#include <iostream>
#include <vector>
#include<iostream>
#include<fstream>
#include<cstdlib>
#include<cstdio>
#include<cmath>
#include<iomanip>


#define NUM_STEPS 300 //越大，曲线越密，越逼近
using namespace std;
using namespace cv;
class CPoint
{
public:
	float x;
	float y;
	CPoint()
	{
		x=0.0;
		y=0.0;
	}
	CPoint(float a,float b)
	{
		x=a;
		y=b;
	}
 
}; 
 
void curve4(vector<CPoint> &p,  
			double x1, double y1,   //Anchor1  
			double x2, double y2,   //Control1  
			double x3, double y3,   //Control2  
			double x4, double y4)   //Anchor2  
{  
	CPoint tmp0(x1,y1);
	p.push_back(tmp0); 
	double dx1 = x2 - x1;  
	double dy1 = y2 - y1;  
	double dx2 = x3 - x2;  
	double dy2 = y3 - y2;  
	double dx3 = x4 - x3;  
	double dy3 = y4 - y3;  
 
	double subdiv_step  = 1.0 / (NUM_STEPS + 1);  
	double subdiv_step2 = subdiv_step*subdiv_step;  
	double subdiv_step3 = subdiv_step*subdiv_step*subdiv_step;  
 
	double pre1 = 3.0 * subdiv_step;  
	double pre2 = 3.0 * subdiv_step2;  
	double pre4 = 6.0 * subdiv_step2;  
	double pre5 = 6.0 * subdiv_step3;  
 
	double tmp1x = x1 - x2 * 2.0 + x3;  
	double tmp1y = y1 - y2 * 2.0 + y3;  
 
	double tmp2x = (x2 - x3)*3.0 - x1 + x4;  
	double tmp2y = (y2 - y3)*3.0 - y1 + y4;  
 
	double fx = x1;  
	double fy = y1;  
 
	double dfx = (x2 - x1)*pre1 + tmp1x*pre2 + tmp2x*subdiv_step3;  
	double dfy = (y2 - y1)*pre1 + tmp1y*pre2 + tmp2y*subdiv_step3;  
 
	double ddfx = tmp1x*pre4 + tmp2x*pre5;  
	double ddfy = tmp1y*pre4 + tmp2y*pre5;  
 
	double dddfx = tmp2x*pre5;  
	double dddfy = tmp2y*pre5;  
 
	int step = NUM_STEPS;  
 
	while(step--)  
	{  
		fx   += dfx;  
		fy   += dfy;  
		dfx  += ddfx;  
		dfy  += ddfy;  
		ddfx += dddfx;  
		ddfy += dddfy;  
		CPoint tmp1(fx,fy);
		p.push_back(tmp1);  
	}  
	CPoint tmp2(x4,y4);
	p.push_back(tmp2); 
}  
 
int main()
{
	CPoint point[4];
        point[0].x=1.0;
        point[0].y=4.0;
        point[1].x=2.2;
        point[1].y=5;
        point[2].x=6;
        point[2].y=3;
        point[3].x=8;
        point[3].y=9;
        point[1].y=(point[0].y*2+point[3].y)/3;
	vector<CPoint> curvePoint;
	curve4(curvePoint,
			point[0].x,point[0].y,
			point[1].x,point[1].y,
			point[2].x,point[2].y,
			point[3].x,point[3].y
			);

	Mat picin(1000,1000,CV_8UC3);
    for(int i=0;i<curvePoint.size();i++)
		{
		//	cout<<"("<<curvePoint[i].x<<","<<curvePoint[i].y<<")"<<endl;
            int i2= int(curvePoint[i].x*100);
            int j= int(curvePoint[i].y*100);
            for (int r=0;r<3;++r)
                picin.at<Vec3b>(i2,j)[r]=255;
		}
    imshow("bez.jpg",picin);
    waitKey();
    imwrite("bez.jpg",picin);
    return 0;
}
