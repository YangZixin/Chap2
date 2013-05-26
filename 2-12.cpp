#include <iostream>
#include <math.h>
#include "TGraph.h"
#include "TAxis.h"

#define phi 30 //当地纬度
#define g 9.8 // m/s^2 重力加速度,与纬度有关
#define omega 7.27e-5  // rad/s 地球自转角速度
#define B_2 4e-5 // m^-1  空气阻力常数

double *velocity_x = NULL;      
double *velocity_y = NULL;      
double *velocity_z = NULL;      //定义速度分量

double *coordinate_x = NULL;      
double *coordinate_y = NULL;      
double *coordinate_z = NULL;    //定义坐标分量

double *time = NULL;            //计时

using namespace std;

void canon(double vx, double vy, double vz)
{
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
	ofstream OutputFile;
	OutputFile.open("./coordinate1000.dat", ios::out);
	for(int j=0; j<i ; j++)
	{
		OutputFile<<coordinate_x[j]<<"   "<<coordinate_y[j]<<"   "<<coordinate_z[j]<<"   "<<endl;
	}
	OutputFile.close();
	double r = - coordinate_z[i-2] / coordinate_z[i-1];
	double final_x = (coordinate_x[i-2] + r * coordinate_x[i-1]) / (r + 1);
	double final_y = (coordinate_y[i-2] + r * coordinate_y[i-1]) / (r + 1);
	cout<<"The final landing coordinate is ("<<final_x<<","<<final_y<<")"<<endl;
//	return i;




}

void plot()
{
/*	double vx_init, vy_init, vz_init;// m/s
	cout<<"The initial velocities in each directions are: ";
	cout<<"vx(0) = ";
	cin>>vx_init;
	cout<<"vy(0) = ";
	cin>>vy_init;
	cout<<"vz(0) = ";
	cin>>vz_init;*/


	TCanvas *c1 = new TCanvas("c1","Canon's Trajectory");


   // create view with axis
   Double_t rx0 = 0, rx1 = 30000, ry0 = 0, ry1 = 30000, rz0 = 0, rz1 = 17000;
   Double_t rmin[3], rmax[3];
   rmin[0] = rx0;
   rmin[1] = ry0;
   rmin[2] = rz0;
   rmax[0] = rx1;
   rmax[1] = ry1;
   rmax[2] = rz1;

//	int numOfSteps;
 //	cout<<"The times of resonance reaction = ";
//	cin>>numOfSteps;

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
 	
	canon(400, 800, 1200);
	TPolyLine3D* gr1 = new TPolyLine3D(12776, coordinate_x, coordinate_y, coordinate_z);

	gr1->SetLineWidth(2);
	gr1->SetLineColor(8);
	gr1->Draw();


	canon(400, 1200, 800);
	TPolyLine3D* gr11 = new TPolyLine3D(9257, coordinate_x, coordinate_y, coordinate_z);

	gr11->SetLineWidth(2);
	gr11->SetLineColor(2);
	gr11->Draw();

	canon(800, 400, 1200);
	TPolyLine3D* gr2 = new TPolyLine3D(12773, coordinate_x, coordinate_y, coordinate_z);

	gr2->SetLineWidth(2);
	gr2->SetLineColor(9);
	gr2->Draw();

	canon(800, 1200, 400);
	TPolyLine3D* gr21 = new TPolyLine3D(5451, coordinate_x, coordinate_y, coordinate_z);

	gr21->SetLineWidth(2);
	gr21->SetLineColor(6);
	gr21->Draw();

	canon(1200, 800, 400);
	TPolyLine3D* gr3 = new TPolyLine3D(5449, coordinate_x, coordinate_y, coordinate_z);

	gr3->SetLineWidth(2);
	gr3->SetLineColor(7);
	gr3->Draw();

	canon(1200, 400, 800);
	TPolyLine3D* gr31 = new TPolyLine3D(9252, coordinate_x, coordinate_y, coordinate_z);

	gr31->SetLineWidth(2);
	gr31->SetLineColor(17);
	gr31->Draw();


}
