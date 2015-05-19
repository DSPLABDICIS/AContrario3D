#include "include/AContrario.hpp"
#include "stdio.h"
//#include "../cluster/include/cluster/cluster.hpp"

using namespace std;

void getxp(boost_incl::mat & xp){
	FILE *fp;
	int i,k,n = 261;
	double Aux;

	int ind[] = {0,1,2,3,6,7,8,9};
	xp.resize(n,10,0);

	fp = fopen("DaniStereo.txt","rt");
	if(fp == NULL)
		throw("Error reading file");

	for(i=0;i<n;i++)                                      // Al rev�s si son enviados de matlab
	{
		for(k=0;k<8;k++)
		{
			fscanf(fp,"%lf",&Aux);
			xp(i,ind[k]) = Aux;
		}
	}

	i = 0;
	cout<<xp(i,ind[0])<<" "<<xp(i,ind[1])<<" "<<xp(i,ind[2])<<" "<<xp(i,ind[3])<<" "<<xp(i,ind[4])<<" "<<xp(i,ind[5])<<endl;
	//       PrintPoint(myV->V[0]);              // mostrar el vector de puntos
	fclose(fp);
}

void printMatrixInfo(boost_incl::mat xp){
	uint i;
	bool _3D;

	_3D = xp.size2() == 8 ? true : false;

	if(_3D)
		cout<<"Cluster\tX\tY\tZ\tVx\tVy\tVz"<<endl;
	else
		cout<<"Cluster\tX\tY\tVx\tVy"<<endl;


	for(i=0; i<xp.size1();i++){
		cout<<xp(i,4)<<"\t";
		cout<<xp(i,0)<<"\t"<<xp(i,1)<<"\t";
		if(_3D)
			cout<<xp(i,6)<<"\t";
		cout<<xp(i,2)<<"\t"<<xp(i,3)<<"\t";
		if(_3D)
			cout<<xp(i,7);
		cout<<endl;
	}
}



int main(void){
	boost_incl::mat xp;
	cluster::AContrario test;
	try{
		cout<<"init"<<endl;
		getxp(xp);
		//test.clusteringTrackingPoints(xp,640,320,0,0.0001); //medidas del test.
		test.clusteringTrackingPoints(xp,5,5,1,3); //medidas del test.
		//printMatrixInfo(xp);
		cout<<"Other method"<<endl<<endl;
		//cluster::clusteringTrackingPoints(xp,640,320); //medidas del test.

	}
	catch(const char* error)
	{
		cout <<error<< endl;
		return -1;
	}
	cout <<"Finished";
	//test.clusteringTrackingPoints()
}


