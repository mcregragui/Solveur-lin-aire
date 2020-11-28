#include "FOM.h"
#include <iostream>
#include <math.h>
#include <iostream>
#include <Eigen>
 
using Eigen::MatrixXd;


using namespace Eigen;
using namespace std;

//constructeur
FOM::FOM(MatrixXd A,VectorXd b) :  A_(A),b_(b)
{	
}
//Méthode de arnoldi
void FOM::Arnoldi(VectorXd u, MatrixXd v, MatrixXd H)
{
    int m=u.size();
    VectorXd v1(m),w(m);
    v1=u/(u.norm());
    double h;
    
    for (int j=0; j<m; j++)
    {
        w=A_*v1;
        for (int k=0; k<m; k++)
        {
            v(j,k)=v1[k];
        }
        for (int i=0; i<j; i++)
        {
            h=w.dot(v1);
            H(i,j)=h;
            w=w-h*v1;
        }
        H(j,j)=w.norm();
    }

}

//gradient à pas optimal
VectorXd FOM::Solve()const								
{
	
	MatrixXd A=A_;
	int n = A.rows();
	VectorXd x(n);											
	for (int i = 0; i < n; i++)
	{
		x[i]=0.;
		
	}

	VectorXd r(n),b(b_);
	//cout<<b<<endl;
	r = b -A*x;
	VectorXd rSuivant(n);
	VectorXd xSuivant(n);
	int j = 0;
	double beta=r.norm();
	int nb_iterat_=0;
	while (beta > pow(10,-10))
	{
		beta=r.norm();
		nb_iterat_=nb_iterat_ +1;
		j++;
	}
	cout<<"le nombre d'itération : "<<nb_iterat_<<endl;
    cout<<"----------------FOM------------------------"<<endl;
	return x;
}