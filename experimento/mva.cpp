//--------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
//--------------------------------------------------------------------
using namespace std;
//--------------------------------------------------------------------
int main( void ){
 int N = 1000;
 double S0 = 0.0027,	V0 = 2, D0 = S0*V0,
		S1 = 0.002,		V1 = 5, D1 = S1*V1,
		S2 = 0.005,		V2 = 2, D2 = S2*V2,
		S3 = 0.0222,	V3 = 1, D3 = S3*V3,
		S4 = 0.0284,	V4 = 2, D4 = S4*V4,
		Q0, Q1, Q2, Q3, Q4,
		R0, R1, R2, R3, R4,
		U0, U1, U2, U3, U4,
		R, X, Z = 4.0;
 stringstream ss;
 ss << setprecision(3);
 ss << "n;R0;R1;R2;R3;R4;R;X;Q0;Q1;Q2;Q3;Q4;U0;U1;U2;U3;U4;" <<
endl;
 Q0 = Q1 = Q2 = Q3 = Q4 = 0.0; // n = 0;
 for( int n = 1; n <= N; n++ ){
 R0 = S0*(1.0 + Q0);
 R1 = S1*(1.0 + Q1);
 R2 = S2*(1.0 + Q2);
 R3 = S3*(1.0 + Q3);
 R4 = S4*(1.0 + Q4);
 R = R0*V0 + R1*V1 + R2*V2 + R3*V3 + R4*V4;
 X = n/(R+Z);
 Q0 = X*R0*V0;
 Q1 = X*R1*V1;
 Q2 = X*R2*V2;
 Q3 = X*R3*V3;
 Q4 = X*R4*V4;
 U0 = Q0/(1+Q0);
 U1 = Q1/(1+Q1);
 U2 = Q2/(1+Q2);
 U3 = Q3/(1+Q3);
 U4 = Q4/(1+Q4);
 ss << n << ";"
 << R0 << ";" << R1 << ";" << R2 << ";" << R3 << ";" << R4
<< ";"
 << R << ";" << X << ";"
 << Q0 << ";" << Q1 << ";" << Q2 << ";" << Q3 << ";" << Q4
<< ";"
 << U0 << ";" << U1 << ";" << U2 << ";" << U3 << ";" << U4
<< endl;
 }
 string str = ss.str();
 cout << str;
 ofstream fo;
 replace(str.begin(), str.end(),'.',',');
 fo.open("MVA-MS.csv");
 fo << str;
 fo.close();
 return 0;
}