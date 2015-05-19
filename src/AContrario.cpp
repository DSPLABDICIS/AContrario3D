/*
 * AContrario.cpp
 *
 *  Created on: Aug 5, 2014
 *      Author: crono21
 */

#include "../include/AContrario.hpp"

namespace cluster {

AContrario::AContrario(): FIRST(0), maxx(0), maxy(0), maxz(0) {

}

AContrario::~AContrario() {
}


void AContrario::freeVectorPoints(VectorPoints *vector){


	if(vector != NULL){
		if(vector->V != NULL){
			int i;
			for(i=0;i<vector->nItems;i++)                                      // Al reves si son enviados de matlab
				if(vector->V[i].V!=NULL)
					free(vector->V[i].V);
		}
		free(vector->V);
	}
	free(vector);
}

void AContrario::writeTrackingPointstoFile(mat & Xp, const char* name){
	int size, i, nItems;
	cout <<Xp.size2()<<endl;
	if(Xp.size2() == 6)
		size = 2;					//x,y
	else if(Xp.size2() == 10)
		size = 3;					//x,y,z
	else
		throw("Error on the size of the moving point matrix \"AContrario::writeTrackingPointstoFile\"");

	nItems = Xp.size1();
	FILE *fp = fopen(name,"w+");

	if(fp==NULL)
		throw("Error Opening the file to write \"AContrario::writeTrackingPointstoFile\"");


	if(size == 3)						//3D Data
		for(i=0;i<nItems;i++)
			fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",Xp(i,6),Xp(i,7),Xp(i,8),Xp(i,2),Xp(i,3),Xp(i,9));
	else
		for(i=0;i<nItems;i++)			//2D Data
			fprintf(fp,"%lf\t%lf\t%lf\t%lf\n",Xp(i,0),Xp(i,1),Xp(i,2),Xp(i,3));
}




int AContrario::clusteringTrackingPoints(mat & Xp,
		int          nC,
		int          nR,
		double 		sigma,
		int			maxz,
		double NFA)
{
	VectorPoints* myV = NULL;   //   Primero se declara un vector de puntos para leer los datos del txt
	VectorPoints* myVel = NULL;               // se inicializa verctor de puntos
	VectorNodos* cluster;  //  vector de grupos soluci�n     // NFA que se plantea < 1 para que no se consideren grupos est�ticos
	int solution, nClusters = 0;

	maxx = nC;
	maxy = nR;
	this->maxz = maxz;

	DataConversion(Xp,&myV,&myVel);

	solution = newTecnica(myV,NFA,sigma,myVel,cluster);  // funci�n que entrega grupos , arg de entrada: vector de puntos(posici�n), NFA y vec de puntos(velocidades)
	nClusters = cluster[solution].nCluster;
	saveClusterstoMatrix(cluster[solution],Xp);

	//cout<<"Solution : "<<solution<<endl;
	//printClustertoScreen(cluster[solution]);


	freeVectorPoints(myV);
	freeVectorPoints(myVel);
	freeVectorList(cluster,5);
	free(cluster);
	//printOnTxt(cluster,toPrint);          // Escribir en txt la soluci�n.
	cout<<"Finished clustering Succesfully"<<endl;

	return nClusters;

}

int AContrario::newTecnica(VectorPoints* myV,double NFA, double sigma,VectorPoints* Vel, VectorNodos *& clusters){
	int Li,may1=0,igual1=0,igual0=0,i,M,MS=0,T,r,nMetrics = 5;
	FuncMetric myMetrics[7] = {Index,Bouldin,Dunns,calinski,Silhouette,RMSstd,Rsquared};
	FuncLinkage myLinkage[5] = {singleLinkage,completeLinkage,medianLinkage,meanLinkage,minAreaLinkage};
	FuncEvaluations myEvaluations[4] = {CLVpromedio,CLVgeometrico,CLVarmonica,CLVstd};

	clusters = (VectorNodos*)malloc(5*sizeof(VectorNodos));
	clusters[0] = moduloAgrupamientoVel(myV,myLinkage[0],NFA,sigma,Vel);    // obtengo las cinco soluciones, o 5 conjuntos de grupos solucion
	clusters[1] = moduloAgrupamientoVel(myV,myLinkage[1],NFA,sigma,Vel);
	clusters[2] = moduloAgrupamientoVel(myV,myLinkage[2],NFA,sigma,Vel);
	clusters[3] = moduloAgrupamientoVel(myV,myLinkage[3],NFA,sigma,Vel);
	clusters[4] = moduloAgrupamientoVel(myV,myLinkage[4],NFA,sigma,Vel);

	for(Li=0;Li<5;Li++){
		if(clusters[Li].nCluster > 1)
			may1 ++;
		else if(clusters[Li].nCluster == 1)
			igual1 ++;
	}

	if(may1 >=3)
	{
		Solutions* mySolutions;  //  { Ssl Bsl RMsl RSsl CHsl ; Scl Bcl RMcl RScl CHcl; Smd Bmd RMmd RSmd CHmd ;Smn Bmn RMmn CHmn ; RSmnSmi Bmi RMmi RSmi CHmi}
		int  nVal=0;

		mySolutions = (Solutions*)malloc(may1*sizeof(Solutions));      //  valor obtenidos de las m�tricas
		for(i=0;i<may1;i++)
			mySolutions[i].evals = (double*)malloc(nMetrics*sizeof(double));
		for(T=0;T<5;T++)
			if(clusters[T].nCluster > 1){
				for(M=0;M<nMetrics;M++)
					mySolutions[nVal].evals[M] = myMetrics[M](clusters[T]);   //   min
				mySolutions[nVal].typeTree = T;
				nVal++;
			}
		MS = finallyEvaluation(mySolutions,may1,nMetrics,myEvaluations[2]);
		for(i=0;i<may1;i++)
			free(mySolutions[i].evals);
		free(mySolutions);
		//return clusters[MS];
		return MS;
	}
	if(igual1 >=3)
	{
		Solutions* mySolutions;  //  { Ssl Bsl RMsl RSsl CHsl ; Scl Bcl RMcl RScl CHcl; Smd Bmd RMmd RSmd CHmd ;Smn Bmn RMmn CHmn ; RSmnSmi Bmi RMmi RSmi CHmi}
		mySolutions = (Solutions*)malloc(igual1*sizeof(Solutions));      //  valor obtenidos de las m�tricas
		int nVal=0;

		for(i=0;i<igual1;i++)
			mySolutions[i].evals = (double*)malloc(2*sizeof(double));

		for(T=0;T<5;T++)
			if(clusters[T].nCluster == 1){
				for(M=5;M<6;M++)
					mySolutions[nVal].evals[M-5] = myMetrics[M](clusters[T]);   //   min
				mySolutions[nVal].typeTree = T;
				nVal++;
			}
		MS = finallyEvaluationToOne(mySolutions,igual1,myEvaluations[2]);
		for(i=0;i<may1;i++)
			free(mySolutions[i].evals);
		free(mySolutions);
		//return clusters[MS];
		return MS;
	}
	clusters[0].nCluster=0;
	//return clusters[FIRST];
	return 0;
}

VectorNodos AContrario::moduloAgrupamientoVel(VectorPoints *myV,double (*f)(Nodo N1,Nodo N2),double NFA,double sigma,VectorPoints* Vel)           // Obtener las distancias uno contra uno de cada punto)
{
	VectorTrees* myTree;   											//   Contiene los nodos ates de ser introducidos al �rbol
	VectorNodos myVnodos;

	myTree = InitTreeVelocity(myV,Vel);             				// Inicializo el vector de arboles con el vector de puntos
	InitMinMax(myTree);                   							// Asignar memoria e Iniciar los valores m�nimos y m�ximos de los grupos
	myTree->Dstd=obtenerDesviacionEstandar(*myV);    				// Obtener las desviacioes estandar de los datos ingresados para cada atributo
	myTree->NFA = NFA ;
	CreateBinaryTree(myTree,(*f));        							//crear el �rbol 0-> no normalizada 1 si la quieres normalizada
	//      printTree(myTree.Tree[FIRST]);
	CreateReferences(myTree,maxx, maxy, maxz);    									// se obtienen las referencias
	//     Select_NFA(myTree.Tree[FIRST],myTree);   				//   hace un reccorido en profundida hasta encontrar un grupo significativo por NFA


	getGralSelection(myTree->Tree[FIRST],*myTree,myTree->Tree[FIRST], sigma);	// Algoritmo Gral donde se agrupan respecto a la NFA y la NFAg
	myVnodos = getSolutionGroups(*myTree);
	return myVnodos;
}

void  AContrario::DataConversion(mat & Xp, VectorPoints ** points,  VectorPoints ** velocities)
{
	int nE,size,k,i;

	(*points) = (VectorPoints*)malloc(sizeof(VectorPoints));
	(*points)->nItems = Xp.size1();
	(*points)->V =(Point*)malloc((*points)->nItems*sizeof(Point));

	(*velocities) =(VectorPoints*)malloc(sizeof(VectorPoints));
	(*velocities)->nItems = Xp.size1();
	(*velocities)->V =(Point*)malloc((*velocities)->nItems*sizeof(Point));

	//cout<<Xp.size2()<<endl;
	if(Xp.size2() == 6)
		size = 2;					//x,y
	else if(Xp.size2() == 10)
		size = 3;					//x,y,z
	else
	{
		throw("Error on the size of the moving point matrix \"AContrario::DataConversion\"");
	}

	(*points)->Size = size;
	(*velocities)->Size = size;

	if(size == 3)
	{
		for(i=0;i<(*points)->nItems;i++)                                      // Al reves si son enviados de matlab
		{
			(*points)->V[i].V=(double*)malloc(size*sizeof(double));
			(*velocities)->V[i].V=(double*)malloc(size*sizeof(double));

			(*points)->V[i].V[0] = Xp(i,6);
			(*points)->V[i].V[1] = Xp(i,7);
			(*points)->V[i].V[2] = Xp(i,8);
			(*velocities)->V[i].V[0] = Xp(i,2);
			(*velocities)->V[i].V[1] = Xp(i,3);
			(*velocities)->V[i].V[2] = Xp(i,9);

			(*points)->V[i].Size=size;
			(*velocities)->V[i].Size=size;

//			cout<<"Point : "<<(*points)->V[i].V[0]<<" "<<(*points)->V[i].V[1]<<" "<<(*points)->V[i].V[2]<<" ";
//			cout<<"Velocity "<<(*velocities)->V[i].V[0]<<" "<<(*velocities)->V[i].V[1]<<" "<<(*velocities)->V[i].V[2]<<" "<<endl;
		}
	}
	else{
		for(i=0;i<(*points)->nItems;i++)                                      // Al reves si son enviados de matlab
		{
			(*points)->V[i].V=(double*)malloc(size*sizeof(double));
			(*velocities)->V[i].V=(double*)malloc(size*sizeof(double));

			(*points)->V[i].V[0] = Xp(i,0);
			(*points)->V[i].V[1] = Xp(i,1);
			(*velocities)->V[i].V[0] = Xp(i,2);
			(*velocities)->V[i].V[1] = Xp(i,3);

			(*points)->V[i].Size=size;
			(*velocities)->V[i].Size=size;
		}

	}
	       //PrintPoint(myV->V[0]);              // mostrar el vector de puntos

}

void AContrario::printClustertoScreen(VectorNodos & cluster){
	int i,j,k;

	cout <<"Number of clusters "<<cluster.nCluster<<endl;
	for(i = 0;i<cluster.nCluster; i++){
		cout<<"Items in cluster "<<i<<": "<<cluster.cluster[i].myInfo.myVp.nItems<<endl;
		cout<<"Pos\tx\ty\tvx\tvy"<<endl;
		for(j =0; j<cluster.cluster[i].myInfo.myVp.nItems ;j++){
			cout<<cluster.cluster[i].myInfo.tags[j]<<"\t";
			for(k=0; k<cluster.cluster[i].myInfo.myVp.V[j].Size;k++)
				cout<<cluster.cluster[i].myInfo.myVp.V[j].V[k]<<"\t";
			for(k=0; k<cluster.cluster[i].myInfo.myVpVelocity.V[j].Size;k++)
				cout<<cluster.cluster[i].myInfo.myVpVelocity.V[j].V[k]<<"\t";

			cout<<endl;
		}
		cout<<endl;
	}
}

void AContrario::saveClusterstoMatrix(VectorNodos & cluster, mat & xp){
	int i,j;

	for(i = 0;i<cluster.nCluster; i++)
		for(j =0; j<cluster.cluster[i].myInfo.myVp.nItems ;j++)
			xp(cluster.cluster[i].myInfo.tags[j],4) = i+1;
}

} /* namespace cluster */
