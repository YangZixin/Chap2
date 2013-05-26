#include <iostream>
#include <math.h>

#define phi 1.570796 //当地纬度
#define g 9.83218 // m/s^2 重力加速度,与纬度有关
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

int main()
{
	double vx_init, vy_init, vz_init;// m/s
/*	cout<<"The initial velocities in each directions are: ";
	cout<<"vx(0) = ";
	cin>>vx_init;
	cout<<"vy(0) = ";
	cin>>vy_init;
	cout<<"vz(0) = ";
	cin>>vz_init;
*/
	int i = canon(500, 500, 500);
	cout<<i;
	return 0;

}
