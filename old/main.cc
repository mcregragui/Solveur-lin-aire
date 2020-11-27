#include<iostream>
#include"Vecteur.h"
#include"Matrix.h"

using namespace std;

int main(){
    Matrix A (2, 2);  
    Matrix B(2*A);                                            
	Vecteur b(2);
	std::cout << "entrez  4 valeurs : "<<std::endl;
	std::cin >> A >> b;
	std::cout << A << b << std::endl;
    A(1,1,3.);
    cout<<A(1,1)<<endl;


};
