#include "../include/mainSLFunctions.hpp"


namespace cluster{

VectorPoints* CreateVectorPoints(int N,int size)                                       //  Asignar memoria a un  Vector de puntos
{
	VectorPoints* myV;

	myV =(VectorPoints*)malloc(sizeof(VectorPoints));
	myV->nItems = N;
	myV->V =(Point*)malloc(myV->nItems*sizeof(Point));
	return myV;
}

//VectorPoints*  ReadPoints()             // Lee los datos desde un .txt y los guarda como un vector de puntos, Lee por coordenada
//{
//	int nE=0,size=0,j,i;
//
//	VectorPoints *myV;
//	nE = 130;
//	size= 4;
//	myV = CreateVectorPoints(nE,size);
//
//	FILE *fp;
//	fp = fopen("prueba3d.txt","rt");
//	for(j=0;j<size;j++)
//		for(i=0;i<nE;i++)
//			fscanf(fp,"%lf",&(myV->V[i].V[j]));
//
//	fclose(fp);
//	return myV;
//}

double DistEuclideanDstd(Point p1,Point p2)   //  vector para obtener la distancia euclidiana entre dos puntos, considera el caso del �ngulo,                                                      //  es decir toma la distancia m�s cercana entre dos �ngulos.
{
	int i;
	double sum=0;

	for(i=0;i<p1.Size; i++)
		sum += (p1.V[i]-p2.V[i])*(p1.V[i]-p2.V[i]);
	return sqrt(sum);
}

Point FuncAverage(Point r,Point l)                                         // Fucion para encontrar coordenadas del nuevo punto por promedio
{                                                          //  En el caso del �ngulo el punto medio m�s cercano a las componentes
	Point aux;
	int i;

	aux.V=(double*)malloc(r.Size*sizeof(double));
	aux.Size = r.Size;
	for(i=0 ; i < r.Size ; i++)
		aux.V[i] = (r.V[i]+l.V[i])/2.0;
	return aux;
}

void GetVectorDistDstd(VectorTrees* myV,double (*f)(Nodo N1,Nodo N2))           // Obtener las distancias uno contra uno de cada punto
{
	int i,j;

	for(i=0;i< myV->nItems;i++)
	{
		myV->Tree[i]->myInfo.myVdist.Dist[i] = 0;
		for(j=i+1;j<myV->nItems;j++)
		{
			myV->Tree[i]->myInfo.myVdist.Dist[j] = (*f)(*(myV->Tree[i]),*(myV->Tree[j]));
			myV->Tree[j]->myInfo.myVdist.Dist[i] = myV->Tree[i]->myInfo.myVdist.Dist[j];
		}
	}
}

void ProcessDistDstd(VectorTrees* myV,Position myP,double (*f)(Nodo N1,Nodo N2))                                    // Acomodar el triangulo de distancias
{
	int i;

	for(i=0;i<myV->nItems;i++)
	{
		myV->Tree[i]->myInfo.myVdist.Dist[myP.col] =  myV->Tree[i]->myInfo.myVdist.Dist[myV->nItems];
		myV->Tree[i]->myInfo.myVdist.nItems--;
		myV->Tree[myP.row]->myInfo.myVdist.Dist[i] =  (*f)(*(myV->Tree[i]),*(myV->Tree[myP.row]));
		myV->Tree[i]->myInfo.myVdist.Dist[myP.row] =  myV->Tree[myP.row]->myInfo.myVdist.Dist[i];
	}
}
void GetAreas(VectorTrees* myV)           // Obtener las distancias uno contra uno de cada punto
{
	int i,j;

	for(i=0;i< myV->nItems;i++)
	{
		myV->Tree[i]->myInfo.myVdist.Dist[i] = 0;
		for(j=i+1;j<myV->nItems;j++)
		{
			//                       myV->Tree[i]->myInfo.myVdist.Dist[j] = getMinArea(*(myV->Tree[i]),*(myV->Tree[j]));
			myV->Tree[j]->myInfo.myVdist.Dist[i] = myV->Tree[i]->myInfo.myVdist.Dist[j];
		}
	}
}


Position GetMinDist(VectorTrees* myV)                                                     // Ordenar las distancias
{
	int j,i;
	Position myP;

	myP.row = 0;
	myP.col = 1;

	for(i=0;i< myV->nItems;i++)
		for(j=i+1;j<myV->nItems;j++)
			if(myV->Tree[myP.col]->myInfo.myVdist.Dist[myP.row] > myV->Tree[i]->myInfo.myVdist.Dist[j])
			{
				myP.col = j;
				myP.row = i;
			}
	return  myP;
}
//-------------------------------------------------------------------------------------
// criterio de cercania entre nodos

double minAreaLinkage(Nodo N1,Nodo N2)   //  Obtiene el area de la caja que contendr� a todos los elementos del nodod N1 y N2
{
	double val=1.0;
	double auxmin, auxMax;
	int i;

	for(i=0 ; i<N1.myInfo.median.Size ; i++)
	{
		if(N1.myInfo.myVp.Max.V[i] > N2.myInfo.myVp.Max.V[i])
			auxMax = N1.myInfo.myVp.Max.V[i];
		else
			auxMax = N2.myInfo.myVp.Max.V[i];
		if(N1.myInfo.myVp.min.V[i] < N2.myInfo.myVp.min.V[i])
			auxmin = N1.myInfo.myVp.min.V[i];
		else
			auxmin = N2.myInfo.myVp.min.V[i];
		val*=(auxMax-auxmin+1);
	}
	return val;
}


double singleLinkage(Nodo N1,Nodo N2)
{
	int i,j;
	double valAux,Daux;

	valAux = DistEuclideanDstd(N1.myInfo.myVp.V[0],N2.myInfo.myVp.V[0]);

	for(i=0;i<N1.myInfo.myVp.nItems;i++)
		for(j=0;j<N2.myInfo.myVp.nItems;j++)
		{
			Daux = DistEuclideanDstd(N1.myInfo.myVp.V[i],N2.myInfo.myVp.V[j]);
			if(Daux < valAux)
				valAux = Daux;
		}
	return valAux;
}



VectorPoints   join2Vectors(VectorPoints vP1,VectorPoints vP2)
{
	VectorPoints   *newVectorPoints;
	int i,j;

	newVectorPoints = CreateVectorPoints(vP1.nItems+vP2.nItems,vP1.V[0].Size);

	for(i=0;i<vP1.nItems;i++)
		newVectorPoints->V[i] = vP1.V[i];

	for(j=vP1.nItems;j<newVectorPoints->nItems;j++)
		newVectorPoints->V[j] = vP1.V[j-vP1.nItems];

	return *(newVectorPoints);

}

double obtenerMedia(VectorPoints myV,int k)
{
	int i;
	double sum=0;

	for(i=0;i<myV.nItems;i++)
		sum+=myV.V[i].V[k];

	return sum/myV.nItems;
}

double* obtenerDesviacionEstandar(VectorPoints myV)
{
	int i,k;
	double *Dstd;
	double Varianza=0.0,mean=0.0;

	Dstd=(double*)malloc(myV.Size*sizeof(double));

	for(k=0;k<myV.Size;k++)
	{
		mean=obtenerMedia(myV,k);
		//      printf("Mean:%lf\t",mean);
		Varianza=0.0;
		for(i=0;i<myV.nItems;i++)
			Varianza+=(myV.V[i].V[k]-mean)*(myV.V[i].V[k]-mean);
		Dstd[k]=sqrt(Varianza/myV.nItems);
	}
	return Dstd;
}

double dStdLinkage(Nodo N1,Nodo N2)
{
	VectorPoints newVectorPoints;
	Point desvN1,desvN1N2;
	double dStd;

	desvN1.Size = N1.myInfo.myPoint.Size;
	desvN1N2.Size = N1.myInfo.myPoint.Size;
	desvN1.V = obtenerDesviacionEstandar(N1.myInfo.myVp);
	newVectorPoints = join2Vectors(N1.myInfo.myVp,N2.myInfo.myVp);
	desvN1N2.V = obtenerDesviacionEstandar(N1.myInfo.myVp);
	dStd = DistEuclideanDstd(desvN1,desvN1N2);
	return dStd;
}

double completeLinkage(Nodo N1,Nodo N2)
{
	int i,j;
	double valAux,Daux;

	valAux = DistEuclideanDstd(N1.myInfo.myVp.V[0],N2.myInfo.myVp.V[0]);

	for(i=0;i<N1.myInfo.myVp.nItems;i++)
		for(j=0;j<N2.myInfo.myVp.nItems;j++)
		{
			Daux = DistEuclideanDstd(N1.myInfo.myVp.V[i],N2.myInfo.myVp.V[j]);
			if(Daux > valAux)
				valAux = Daux;
		}
	return valAux;
}

double meanLinkage(Nodo N1,Nodo N2)
{
	return DistEuclideanDstd(N1.myInfo.mean,N2.myInfo.mean);
}

double medianLinkage(Nodo N1,Nodo N2)
{
	return DistEuclideanDstd(N1.myInfo.median,N2.myInfo.median);
}
}//end namespace
