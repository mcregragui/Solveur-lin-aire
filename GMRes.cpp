#include "GMRes.h"
#include <iostream>
#include <math.h>
#include <iostream>
#include <Eigen>
 
using Eigen::MatrixXd;


using namespace Eigen;
using namespace std;

//constructeur
GMRes::GMRes(MatrixXd A,VectorXd b) :  A_(A),b_(b)
{	
}
//Méthode de arnoldi
void GMRes::Arnoldi(VectorXd u, MatrixXd& VV, MatrixXd& HH)
{
	//cout<<"started"<<endl;
    int m=u.size();
	//cout<<m<<endl;
    VectorXd v1(m),w(m);
	MatrixXd V(m,m+1),H(m+1,m);
	//cout<<"done"<<endl;
    V.col(0)=u/(u.norm());
    double h;
    //cout<<"loop"<<endl;
    for (int j=0; j<m; j++)
    {
        
        for (int i=0; i<=j; i++)
        {
            H(i,j)=(A_*V.col(j)).dot(V.col(i));//je doute cette expression
        }
		VectorXd sum(m);

		for (int i = 0; i < m; i++)
		{
			sum[i]=0.;
		}

		for (int i=0; i<j; i++)//i<m
		{
			sum=sum+H(i,j)*V.col(i);
		}
		w=A_*V.col(j)-sum;
        H(j+1,j)=w.norm();
		if (H(j+1,j)==0.)
		{
			break;
		}
		V.col(j+1)=w/H(j+1,j);

    }
	VV=V;
	HH=H;

}


//décomp QR
//On utilise la décomposition QR à l'aide de méthode de Given, on trouve d'abord y tel que ||beta*e1-Hm*y||
//soit minimale, Hm=QR tel que Q est Q^(t)*Q=I et R triang sup.
void GMRes::QR(MatrixXd& Q, MatrixXd& R)
{
    int m=A_.rows();
	MatrixXd G(m,m);
	Q=identity(m,m);
	R=A_;
    for (int k=0; k<m; k++)
    {
		double c,s;
        for(j=m-1;j>k;k--)
		{
			/*for(i=0;i<m;i++)
			{
				G(i,i)=1;
			}*/
			G=identity(m,m);
			if (R(k,j)==0){
				c=1;
				s=0;
			}
			else
			{
				if (abs(R(k,j))>abs(R(k-1,j)))
				{
					s=1./(1+pow(R(k,j)/R(k-1,j),2));
					c=s*R(k,j)/R(k-1,j);
				}
				else 
				{
					c=1./(1+pow(R(k,j)/R(k-1,j),2));
					s=c*R(k,j)/R(k-1,j);
				}
			}
			G(j-1,j-1)=c;
			G(j-1,j)=s;
			G(j,j-1)=-s;
			G(j,j)=c;
			R=G.transpos()*R;
			Q=Q*G;
		}
    }


}

//extraction de R //fonctionne
MatrixXd GMRes::ExtractR(MatrixXd C)
{
    int n=C.rows();
    MatrixXd R(n,n);
    for (int i=0; i<n ;i++)
    {
        for (int j=0; j<n ; j++)
        {
            if (j>=i)
            {
                R(i,j)=C(i,j);
            }
            else
            {
                R(i,j)=0.;
            }
            
        }
    }
    return R;
}



//GMRes
VectorXd GMRes::Solve()								
{

	
	
	const int n = A_.rows();
	VectorXd x(n);											
	for (int i = 0; i < n; i++)
	{
		x[i]=0.;
		
	}

	VectorXd r(n);
	//cout<<b<<endl;
	r = b_ -A_*x;
	VectorXd v1(n);
	v1=r/r.norm();
	VectorXd rSuivant(n);
	VectorXd xSuivant(n);
	rSuivant=r;
	xSuivant=x;
	double beta=r.norm();
	int nb_iterat_=0;
	while (beta > pow(10,-10))
	{

		nb_iterat_=nb_iterat_ +1;

		//Arnoldi
		MatrixXd HH(n+1,n),VV(n,n+1);
		this -> Arnoldi(rSuivant,VV,HH);
		MatrixXd H(n,n),V(n,n);
		
		H = HH.block(0,0,n,n);
		V= VV.block(0,0,n,n);
		//cout<<"----------------La matrice H------------------------"<<endl;
		//cout<<H<<endl;
		//cout<<"----------------La matrice V------------------------"<<endl;
		//cout<<V<<endl;
		
		//décomp QR
		HouseholderQR<MatrixXd> qr(n,n);
		qr.compute(H);
        MatrixXd Q(n,n),R(n,n);
        MatrixXd res(n,n);
        Q = qr.householderQ();
        res=MatrixXd::Zero(n,n);
        res=qr.matrixQR();
        R= this ->ExtractR(res);

        //jusqu'à ici extraction et QR fonctionne

		//MatrixXd Q = qr.householderQ();

		//cout<<"----------------La matrice Q------------------------"<<endl;
		//cout<<Q<<endl;


		//itération
		VectorXd e(n);
		e= V.transpose()*v1;
		VectorXd y(n);
		y = H.householderQr().solve(beta*e);
		//cout<<"----------------La solution de Q------------------------"<<endl;
		cout<<y<<endl;
		
		xSuivant = xSuivant+ V*y;
		rSuivant = b_ -A_*xSuivant;
		//rSuivant=(y.dot(e))*VV.col(n);//*VV.col(n+1);
		

		x=xSuivant ;
		r=rSuivant;


		v1=rSuivant/(rSuivant.norm());


		
		beta=rSuivant.norm();
		if (beta > pow(10,-9))
		{
			break;
		}
		
	}

    cout<<"----------------FOM------------------------"<<endl;
	cout<<"le nombre d'itération : "<<nb_iterat_<<endl;
	return x;
}