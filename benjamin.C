#include <iostream>
#include <math.h>
// #include "TGraph.h"
// #include "TAxis.h"

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

int canon(double vx, double vy, double vz)
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
	double r = - coordinate_z[i-2] / coordinate_z[i-1];
	double final_x = (coordinate_x[i-2] + r * coordinate_x[i-1]) / (r + 1);
	double final_y = (coordinate_y[i-2] + r * coordinate_y[i-1]) / (r + 1);
	cout<<"The final landing coordinate is ("<<final_x<<","<<final_y<<")"<<endl;
	ofstream OutputFile;
	OutputFile.open("./coordinate.dat", ios::out);
	for(int j=0; j<i ; j++)
	{
		OutputFile<<coordinate_x[j]<<"   "<<coordinate_y[j]<<"   "<<coordinate_z[j]<<"   "<<endl;
	}
	OutputFile.close();

	return i;




}

void plot()
{
	double vx_init, vy_init, vz_init;// m/s
	cout<<"The initial velocities in each directions are: ";
	cout<<"vx(0) = ";
	cin>>vx_init;
	cout<<"vy(0) = ";
	cin>>vy_init;
	cout<<"vz(0) = ";
	cin>>vz_init;

	int i = canon(vx_init, vy_init, vz_init);
	cout<<i;
	int numOfSteps;
 	cout<<"The times of resonance reaction = ";
	cin>>numOfSteps;
	

	TCanvas *c1 = new TCanvas("c1","Canon's Trajectory");

   // create view with axis
   Double_t rx0 = 0, rx1 = 40000, ry0 = 0, ry1 = 30000, rz0 = 0, rz1 = 10000;
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


	TPolyLine3D* gr1 = new TPolyLine3D(numOfSteps, coordinate_x, coordinate_y, coordinate_z);
	gr1->SetLineWidth(1);
	gr1->SetLineColor(42);
	gr1->Draw();

}
