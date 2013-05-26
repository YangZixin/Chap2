#include <iostream>
#include <math.h>
#include "TGraph.h"
#include "TAxis.h"

//#define phi 1.570796 //当地纬度
//#define g 9.83218 // m/s^2 重力加速度,与纬度有关
#define omega 7.27e-5  // rad/s 地球自转角速度
#define B_2 4e-5 // m^-1  空气阻力常数
#define PI 3.14159265
double *velocity_x = NULL;      
double *velocity_y = NULL;      
double *velocity_z = NULL;      //定义速度分量

double *coordinate_x = NULL;      
double *coordinate_y = NULL;      
double *coordinate_z = NULL;    //定义坐标分量

double *time = NULL;            //计时

using namespace std;

int canon(double PHI, double g)
{
	double phi;
	phi= PHI * PI / 180;
	double vx=500;
	double vy=500;
	double vz=500;
	double x, y, z, t, v_tot, ax, ay, az;
	double dt = 0.01;// s（单位秒）
	int i=0;                   //计数
	int numOfSteps = (2 * vz / (g - 2 * omega * vy)) / dt + 1;//大概估算一下所需时间

	velocity_x = new double[numOfSteps];
	velocity_y = new double[numOfSteps];
	velocity_z = new double[numOfSteps];

	coordinate_x = new double[numOfSteps];
	coordinate_y = new double[numOfSteps];
	coordinate_z = new double[numOfSteps];

	time = new double[numOfSteps];

	while (z >= 0)
	{
		velocity_x[i] = vx;
		velocity_y[i] = vy;
		velocity_z[i] = vz;

		v_tot = sqrt(vx * vx + vy * vy + vz * vz);

		coordinate_x[i] = x;
		coordinate_y[i] = y;
		coordinate_z[i] = z;
//		cout<<coordinate_x[i]<<" "<<coordinate_y[i]<<" "<<coordinate_z[i]<<endl;
		time[i] = t;

		ax = 2 * omega * vy * sin(phi) - B_2 * v_tot * vx;
		ay = -2 * omega *(vz * cos(phi) + vx * sin(phi)) - B_2 * v_tot * vy;
		az = -g + 2 * omega * vy * cos(phi) - B_2 * v_tot * vz;

		x += vx * dt + ax * dt * dt / 2;
		y += vy * dt + ay * dt * dt / 2;
		z += vz * dt + az * dt * dt / 2;

		vx += ax * dt;
		vy += ay * dt;
		vz += az * dt;

		t += dt;
		i++;
	}//炮弹发射==>嘭！

	return i;




}

void plot()
{
	

	TCanvas *c1 = new TCanvas("c1","Canon's Trajectory");

   // create view with axis
   Double_t rx0 = 18500, rx1 = 19200, ry0 = 18500, ry1 = 19200, rz0 = 0, rz1 = 1000;
   Double_t rmin[3], rmax[3];
   rmin[0] = rx0;
   rmin[1] = ry0;
   rmin[2] = rz0;
   rmax[0] = rx1;
   rmax[1] = ry1;
   rmax[2] = rz1;


   TView3D *view = new TView3D(1, rmin, rmax);
   view->ShowAxis();
   TAxis3D *axis = TAxis3D::GetPadAxis(); // Get pointer to axis
   if (axis) {
      axis->SetLabelSize(0.02);    
      axis->SetLabelOffset(-0.02, "z"); 
      axis->SetLabelColor(kBlue);  
      axis->SetAxisColor(kBlue); 

      axis->SetXTitle("East-x");
      axis->SetYTitle("South-y");
      axis->SetZTitle("Altitude-z");
   }
   // draw a box around

	canon(0, 9.78030);
	TPolyLine3D* gr0 = new TPolyLine3D(7461, coordinate_x, coordinate_y, coordinate_z);
	gr0->SetLineWidth(1);
	gr0->SetLineColor(10);
	gr0->Draw();

	canon(10, 9.78186);
	TPolyLine3D* gr1 = new TPolyLine3D(7460, coordinate_x, coordinate_y, coordinate_z);
	gr1->SetLineWidth(1);
	gr1->SetLineColor(1);
	gr1->Draw();

	canon(20, 9.78634);
	TPolyLine3D* gr2 = new TPolyLine3D(7456, coordinate_x, coordinate_y, coordinate_z);
	gr2->SetLineWidth(1);
	gr2->SetLineColor(2);
	gr2->Draw();

	canon(30, 9.79321);
	TPolyLine3D* gr3 = new TPolyLine3D(7450, coordinate_x, coordinate_y, coordinate_z);
	gr3->SetLineWidth(1);
	gr3->SetLineColor(3);
	gr3->Draw();


	canon(40, 9.80166);
	TPolyLine3D* gr4 = new TPolyLine3D(7442, coordinate_x, coordinate_y, coordinate_z);
	gr4->SetLineWidth(1);
	gr4->SetLineColor(4);
	gr4->Draw();

	canon(50, 9.81066);
	TPolyLine3D* gr5 = new TPolyLine3D(7433, coordinate_x, coordinate_y, coordinate_z);
	gr5->SetLineWidth(1);
	gr5->SetLineColor(5);
	gr5->Draw();

	canon(60, 9.81914);
	TPolyLine3D* gr6 = new TPolyLine3D(7424, coordinate_x, coordinate_y, coordinate_z);
	gr6->SetLineWidth(1);
	gr6->SetLineColor(6);
	gr6->Draw();

	canon(70, 9.82606);
	TPolyLine3D* gr7 = new TPolyLine3D(7416, coordinate_x, coordinate_y, coordinate_z);
	gr7->SetLineWidth(1);
	gr7->SetLineColor(7);
	gr7->Draw();

	canon(80, 9.83058);
	TPolyLine3D* gr8 = new TPolyLine3D(7409, coordinate_x, coordinate_y, coordinate_z);
	gr8->SetLineWidth(1);
	gr8->SetLineColor(8);
	gr8->Draw();

	canon(90, 9.83218);
	TPolyLine3D* gr9 = new TPolyLine3D(7403, coordinate_x, coordinate_y, coordinate_z);
	gr9->SetLineWidth(1);
	gr9->SetLineColor(9);
	gr9->Draw();




}
