#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <time.h>
#include <math.h>

//all distances are in mm
//all angles are in radians

#define pi 3.1415926535897
#define a_small 1.5
#define R 3000
#define THETA_STEP .001
#define MAX_ANGLE (int)(pi/2/THETA_STEP)
#define PATHS_PER_ANGLE 1000
#define N 2

using namespace std;

///////////////////////////////////////////
//Global Variables
///////////////////////////////////////////



double a_big;
double theta_max;
double L;
double theta_cone;
double theta_cutoff;
double n;
double G;
double b;

//Cone Geometry

//Wall 1 : m = -tan(theta_cone)
//         g = -a_small*tan(theta_cone/2)
//Wall 2 : m =  tan(theta_cone)
//         g = -a_small*tan(theta_cone/2)


///////////////////////////////////////////
// Bounce Calculations
///////////////////////////////////////////


void find_intersection(double m,double g, double m1, double g1, double m2, double g2, double intersection[]){
double x1,x2,y1,y2;

x1 = (g1 - g) / (m - m1);
x2 = (g2 - g) / (m - m2);
y1 = m1*x1 + g1;
y2 = m2*x2 + g2;

if(-g/m < a_small && -g/m > -a_small){  //If the ray intersects the collector (y=0) within the -a_small < x < a_small range
intersection[0] = -g/m;
intersection[1] = 0;
}
else if(y1 < 0 ){
intersection[0] = x1;
intersection[1] = y1;
}
else if (y2 < 0){
intersection[0] = x2;
intersection[1] = y2;
}
else if (y1 > y2){
intersection[0] = x2;
intersection[1] = y2;
}
else{
intersection[0] = x1;
intersection[1] = y1;
}


return;
}


bool bounce(double m, double g){

double m_new;
double g_new;
double theta = pi/2 - atan(m);

double angle_to_surf;

double m1 = -tan(theta_cone);
double m2 = -m1;
double g1 = -a_small*tan(theta_cone/2);
double g2 = g1;

double intersection[2]; // 0 -> x    1 -> y
find_intersection(m, g, m1, g1, m2, g2, intersection);


angle_to_surf = theta + theta_cone/2;
double theta_new = theta + theta_cone/2;

m_new = tan(pi/2 - theta_new);
g_new = intersection[1] - (a_big*m_new);

if(intersection[1] <= 0){
return true;
}

if(intersection[1] >= L){
return false;
}


return bounce (m_new, g_new);
}



///////////////////////////////////////////
//Utility Functions
///////////////////////////////////////////

double radtodeg(double rad){
return rad*180/pi;
}

double degtorad(double deg){
return deg/180*pi;
}



///////////////////////////////////////////
//Main Program
///////////////////////////////////////////


int main(){

ofstream fout;
fout.open("Output.txt");
ofstream csvout;
csvout.open("Output.csv");

fout << "Starting Cone Simulation..."  << endl << endl;
cout << "Starting Cone Simulation..."  << endl << endl;

theta_max = degtorad(25);

///////////////////////////////////////////
//Calculate the Geometry and Parameters
///////////////////////////////////////////

a_big = a_small/(sin(theta_max));
L = (a_small+a_big)/(tan(theta_max));
theta_cone = 2*atan((a_big - a_small)/L);
b = L/(a_big*a_big - a_small*a_small);
n = a_big/a_small;
G = n*n;
theta_cutoff = atan(a_big/R);

cout << "Theta Max: " << radtodeg(theta_max) << endl;
cout << "Small Aperture: " << a_small << endl;
cout << "Big Aperture: " << a_big << endl;
cout << "Distance from Source: " << R << endl;
cout << "Collector Length: " << L << endl;
cout << "Theta Cone: " << radtodeg(theta_cone) << endl;
cout << "Parabola Coefficient (b): " << b << endl;
cout << "Aperture Ratio: " << n << endl;
cout << "Area Ratio: " << G << endl;
cout << "Theta Cutoff: " << radtodeg(theta_cutoff) << endl;

fout << "All Parameters Calculated"  << endl << endl;
cout << "All Parameters Calculated"  << endl << endl;

///////////////////////////////////////////
//Normalize The Power
///////////////////////////////////////////

//i*THETA_STEP is the actual angle

fout << "Normalizing Power Model..."  << endl << endl;
cout << "Normalizing Power Model..."  << endl << endl;

double sum = 0;
double power[MAX_ANGLE];
double weight[MAX_ANGLE];

for(int i = 0; i < MAX_ANGLE; i++){
power[i] = pow(cos(i*THETA_STEP),N);
sum +=pow(cos(i*THETA_STEP),N);
}

for(int i =0; i < MAX_ANGLE; i++){
power[i] /= sum;
}

for(int i = 0; i < MAX_ANGLE ; i++){
    fout << "Angle : " << i*THETA_STEP << " Power: " << power[i] << endl;
    csvout << i*THETA_STEP << ", " << power[i] << endl;
}



///////////////////////////////////////////
//Simulate (Ray-Tracing)
///////////////////////////////////////////

// y = mx + g

// m = tan(pi/2 - theta) * a_big
// g = L - (a_big*m) + path_height

double g,m;
double theta;
int total;

fout << "Starting Ray Tracing..."  << endl << endl;
cout << "Starting Ray Tracing..."  << endl << endl;

for(int i = 0; i < MAX_ANGLE; i++){ //For Each Theta i
theta = i*THETA_STEP; // theta incident
//fout << "Theta: "  << theta << endl;
//cout << "Theta: "  << theta << endl;
total=0;
for(int j = 0; j< PATHS_PER_ANGLE ; j++){  //For Each (of the 1000) paths per Theta i

if(i*THETA_STEP!=0){
m = tan(pi/2 - theta);
g = L - (a_big*m) + j/(2*a_big*m);
//cout << "Theta: "  << theta << "m: " << m <<"g: " << g << endl;

    if(bounce(m,g)){
    total++;
    }
}
else
{
total++;
}

}
weight[i] = ((double)total)/PATHS_PER_ANGLE;
//cout << weight[i] << endl;
}


///////////////////////////////////////////
//Integrate Results (Sum of Multiplication)
///////////////////////////////////////////

fout << "Summing the Multiplication..."  << endl << endl;
cout << "Summing the Multiplication..."  << endl << endl;

double power_captured[2*MAX_ANGLE];
double total_power=0;

for(int i =0; i < MAX_ANGLE; i++){
power_captured[i] = power[i]*weight[i];
total_power += power_captured[i];
}


///////////////////////////////////////////
//Interrupt Results
///////////////////////////////////////////


//power over all transmitted
//power over power without collector
//power over power incident upon collector

cout << "Power Captured: " << total_power << endl;

//These are linear powers for comparison only
//For accurate 3D Power Simulations we need to consider an integral

fout.close();
csvout.close();
return 0;
}
