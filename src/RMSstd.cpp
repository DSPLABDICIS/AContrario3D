#include "../include/RMSstd.hpp"

namespace cluster{

double NormalHasting(double sigma, double media, double x)
{
	double X=x;
	if(x<0)
		X=-x;
	int i=0;
	double N = 0.0;
	double P=  0.2316419;
	double Z= (X-media)/sigma;
	double T= 1.0/(1.0+P*Z);
	double b[5]= {0.319381530,-0.356563782,1.781477937,-1.821255978,1.330274429};
	double Prod=0.0;

	for(i=0;i<5;i++)
		Prod += b[i]*pow(T,i+1);
	//if(x<0)
		// N =2*(1 - (exp(-Z*Z/2))*Prod/(sqrt(2*PI)));
	// else
	N =2*(exp(-Z*Z/2))*Prod/(sqrt(2*PI));
	// N = (exp(-((x-media)/sigma)*((x-media)/sigma)/2.0))/(sigma*sqrt(2*PI));
	return N;
}

double RMSstdVp(VectorPoints* vN,int n)
{
	int i,j,P,k;

	double suma=0,suma2=0,dist;
	Point myP;

	for(k=0;k < n;k++){
		myP.Size = vN[k].V->Size;
		myP.V =(double*)malloc(myP.Size*sizeof(double));

		myP.V[0]= obtenerMedia(vN[k],0);
		myP.V[1]= obtenerMedia(vN[k],1);
		myP.V[2]= obtenerMedia(vN[k],2);
	}

	for(i = 0 ; i < n ; i++)
	{
		for( j = 0 ; j < vN[i].nItems ; j++)
		{
			dist = DistEuclideanDstd(vN[i].V[j],myP);
			suma += dist*dist;
		}
		suma2 = vN[i].nItems - 1;
	}
	suma=sqrt(suma/(myP.Size*suma2));
	//              printf("suma=%lf\n",suma);
	free(myP.V);

	return suma;
}



double RMSstd( VectorNodos vN)
{
	int i,j,P;

	double suma=0,suma2=0,dist;

	P = vN.cluster[0].myInfo.mean.Size;
	for(i = 0 ; i < vN.nCluster ; i++)
	{
		for( j = 0 ; j < vN.cluster[i].myInfo.myVp.nItems ; j++)
		{
			dist = DistEuclideanDstd(vN.cluster[i].myInfo.myVp.V[j],vN.cluster[i].myInfo.mean);
			suma += dist*dist;
		}
		suma2 += vN.cluster[i].myInfo.myVp.nItems - 1;
	}
	return suma = sqrt(suma/(P*suma2));
}

double Dunns(VectorNodos vN)
{
	double distMin=0,distMax=0;
	double auxMin=0,auxMax=0;
	int i,j,x,y,flag=0,x2;

	for(i = 0 ; i < vN.nCluster ; i++)
	{
		for(x = 0 ; x < vN.cluster[i].myInfo.myVp.nItems ; x++)
		{
			for(x2 = x+1 ; x2<vN.cluster[i].myInfo.myVp.nItems; x2++)
			{
				distMax = DistEuclideanDstd(vN.cluster[i].myInfo.myVp.V[x],vN.cluster[i].myInfo.myVp.V[x2]);
				if(distMax > auxMax)
					auxMax = distMax;
			}
			for(j = i+1 ; j < vN.nCluster ; j++)
				for(y = 0 ; y < vN.cluster[j].myInfo.myVp.nItems ; y++){
					distMin = DistEuclideanDstd(vN.cluster[i].myInfo.myVp.V[x],vN.cluster[j].myInfo.myVp.V[y]);
					if(distMin < auxMin || flag==0)
					{
						auxMin = distMin;
						flag=1;
					}

				}
		}
	}
	//     printf("DunsMin=%lf\tDunsMax=%lf\tDuns=%lf\n",auxMin,auxMax,auxMin/auxMax);
	return auxMin/auxMax;
}

double Silhouette(VectorNodos vN)
{
	int i,j,n,i2,k,N,flag;
	double suma=0,suma2=0,Dint,Dext,ai,bi,auxVal;

	for(n = 0 ; n < vN.nCluster ; n++)
	{
		suma=0;
		for( i = 0 ; i < vN.cluster[n].myInfo.myVp.nItems ; i++)
		{
			N = 0;  Dint = 0;
			for( j = 0 ; j < vN.cluster[n].myInfo.myVp.nItems ; j++)                                          // a(i)
                        		Dint+= DistEuclideanDstd(vN.cluster[n].myInfo.myVp.V[i],vN.cluster[n].myInfo.myVp.V[j]);

			//  calcula sum de distancias del elemneto i del cluster n contra los elementos de los dem�s clusters
			flag=0;
			auxVal=0;
			for(k = 0 ; k < vN.nCluster ; k++)
				if(k!=n)
				{
					Dext = 0;
					for( i2 = 0 ; i2 < vN.cluster[k].myInfo.myVp.nItems ; i2++)
						Dext+= DistEuclideanDstd(vN.cluster[n].myInfo.myVp.V[i],vN.cluster[k].myInfo.myVp.V[i2]);
					Dext /=  vN.cluster[k].myInfo.myVp.nItems;
					if(Dext<auxVal || flag==0)
						auxVal=Dext;                    //  min d(x,y)
					flag=1;
				}

			ai = Dint/(vN.cluster[n].myInfo.myVp.nItems - 1);
			bi = auxVal;
			if(bi > ai)
				suma+= 1-ai/bi;
			if(bi < ai)
				suma+= bi/ai-1;
		}
		suma2 += suma/vN.cluster[n].myInfo.myVp.nItems;
	}

	return   suma2 = suma2/vN.nCluster;
}

double Bouldin(VectorNodos vN)
{
	int i,n,m,k;
	double Z,Zaux,suma2=0,aux;
	double DistMedias[vN.nCluster],Dint;
	//-------------------------------------------------------------------------
	//calcular el promedio de diferencias entre los puntos y la media de su grupo
	for(n = 0 ; n < vN.nCluster ; n++)
	{
		Dint=0;
		for( i = 0 ; i < vN.cluster[n].myInfo.myVp.nItems ; i++)
			Dint+= DistEuclideanDstd(vN.cluster[n].myInfo.myVp.V[i],vN.cluster[n].myInfo.mean); // d(i,mean)
			DistMedias[n] = Dint/(vN.cluster[n].myInfo.myVp.nItems);
	}
	//---------------------------------------------------------------------------------------------------

	//___________________________________________________________________________________
	// obtener la m�xima
	for(m = 0 ; m < vN.nCluster ; m++)
	{
		Z = 0;
		for(k = 0 ;  k < vN.nCluster ; k++)
			if(k!=m)
			{
				aux = DistEuclideanDstd(vN.cluster[m].myInfo.mean,vN.cluster[k].myInfo.mean);
				Zaux = (DistMedias[m]+ DistMedias[k])/aux;
				if(Zaux>Z)
					Z=Zaux;
			}
		suma2 += Z;           // sumatoria de m�ximos
	}
	return   suma2 = suma2/vN.nCluster;   //  DIvisi�n entre numero de clusters
}

void raiseMemory(TreeBinary myV,VectorNodos* myVnodos)
{
	if(myV->flagNFA==1)
	{
		myVnodos->nCluster++;
		myVnodos->cluster = (Nodo*) realloc (myVnodos->cluster,myVnodos->nCluster*sizeof(Nodo));
		myVnodos->cluster[myVnodos->nCluster - 1] = *myV;
	}
}

void  createVectorGroups(TreeBinary myNodo,VectorNodos* myVnodos, void (*G)(TreeBinary B,VectorNodos* Vn))    // REcorrido del �rbol mostrando algo dadoo funci�n G
{
	int count=0;
	TreeBinary aux = NULL;

	//aux = (TreeBinary)malloc(sizeof(Nodo)); //todo: is this needed??
	aux =  myNodo;

	if(aux->lNodo == NULL)
		return;

	do
	{
		(*G)(aux,myVnodos);                          //  funcion para imprimir
		aux = aux->lNodo;
		count++;
		if( aux->myInfo.myVp.nItems == 1)  // aqui va la condicion
		{
			aux = aux->fNodo;
			while( count>0)
			{
				createVectorGroups(aux->rNodo,myVnodos,(*G));
				aux = aux->fNodo;
				count--;
			}
			return;
		}
	}
	while(aux->lNodo != NULL);
	//free(aux);
	return;
}


VectorNodos getSolutionGroups(VectorTrees& myTree)  // Obtener los par�metros de evaluaci�n
{
	VectorNodos myVnodos;    //  Pra guardar Nodos significativos
	myVnodos.cluster=NULL;
	myVnodos.nCluster=0;
	myVnodos.Tree = &myTree;

	createVectorGroups(myTree.Tree[0],&myVnodos,raiseMemory);
	return myVnodos;
}

double calinski(VectorNodos vN)
{
	double Calins=0.0;
	int NC=vN.nCluster,nItemsN=0;

	if(NC==0)
		return Calins;

	int i,k,k2,j,j2;
	Point center;
	double sumSup=0,sumInf,sumFinal=0,Vaux,auxInf;

	int sizeN = vN.cluster[0].myInfo.mean.Size;
	center.Size = sizeN;
	center.V = (double*)calloc(center.Size,sizeof(double));

	for(i=0;i<NC;i++)
	{
		nItemsN+= vN.cluster[i].myInfo.myVp.nItems;
		for(k=0;k<sizeN;k++)
			center.V[k] += vN.cluster[i].myInfo.mean.V[k]*vN.cluster[i].myInfo.myVp.nItems;
	}
	for(k2=0;k2<sizeN;k2++)
		center.V[k2] = center.V[k2]/nItemsN;

	for(j=0;j<NC;j++)
	{
		Vaux = DistEuclideanDstd(vN.cluster[j].myInfo.mean,center);
		sumSup += Vaux*Vaux*vN.cluster[j].myInfo.myVp.nItems;

		sumInf=0;
		for(j2=0;j2<vN.cluster[j].myInfo.myVp.nItems;j2++)
		{
			auxInf =  DistEuclideanDstd(vN.cluster[j].myInfo.myVp.V[j2],vN.cluster[j].myInfo.mean);
			sumInf+= auxInf*auxInf;
		}
		sumFinal+= sumInf;
	}
	//                    printf("sup=%lf\tinf:%lf",sumSup,sumFinal);
	Calins = sumSup*(nItemsN-NC)/(sumFinal*(NC-1));
	free(center.V);
	return Calins;
}

double Rsquared(VectorNodos vN)
{
	double Rs=0;
	int NC=vN.nCluster,nItemsN=0;
	if(NC==0)
		return Rs;

	int i,k,k2,j,j2;
	Point center;
	double sumInf,sumFinal=0,sumT,auxT,auxInf;

	int sizeN = vN.cluster[0].myInfo.mean.Size;

	center.Size = sizeN;
	center.V = (double*)calloc(center.Size,sizeof(double));

	for(i=0;i<NC;i++)
	{
		nItemsN+= vN.cluster[i].myInfo.myVp.nItems;
		for(k=0;k<sizeN;k++)
			center.V[k] += vN.cluster[i].myInfo.mean.V[k]*vN.cluster[i].myInfo.myVp.nItems;
	}
	for(k2=0;k2<sizeN;k2++)
		center.V[k2]= center.V[k2]/nItemsN;
	sumT=0;
	for(j=0;j<NC;j++)
	{
		sumInf=0;
		for(j2=0;j2<vN.cluster[j].myInfo.myVp.nItems;j2++)
		{
			auxInf =  DistEuclideanDstd(vN.cluster[j].myInfo.myVp.V[j2],vN.cluster[j].myInfo.mean);
			sumInf+= auxInf*auxInf;
			auxT =  DistEuclideanDstd(vN.cluster[j].myInfo.myVp.V[j2],center);
			sumT+= auxT*auxT;
		}
		sumFinal+= sumInf;
	}

	Rs = 1 - sumFinal/sumT;
	free(center.V);
	return Rs;
}

double Index(VectorNodos vN)
{
	double I=0;
	int NC=vN.nCluster,nItemsN=0;

	if(NC==0)
		return I;

	int i,k,k2,j,j2,j3;
	Point center;
	double sumInf=0,sumFinal=0,sumT=0,auxT=0,auxInf=0;
	int sizeN = vN.cluster[0].myInfo.mean.Size;

	center.Size = sizeN;
	center.V = (double*)calloc(center.Size,sizeof(double));

	for(i=0;i<NC;i++)
	{
		nItemsN+= vN.cluster[i].myInfo.myVp.nItems;
		for(k=0;k<sizeN;k++)
			center.V[k] += vN.cluster[i].myInfo.mean.V[k]*vN.cluster[i].myInfo.myVp.nItems;
	}
	for(k2=0;k2<sizeN;k2++)
		center.V[k2]= center.V[k2]/nItemsN;

	int flag=0;
	double maxD=0,auxMaxD=0;

	for(j=0;j<NC;j++)
	{
		for(j3=j+1;j3<NC;j3++)
		{
			maxD = DistEuclideanDstd(vN.cluster[j].myInfo.mean,vN.cluster[j3].myInfo.mean);
			if(maxD > auxMaxD || flag==0)
			{
				auxMaxD = maxD;
				flag=1;
			}
		}

		sumInf=0;
		for(j2=0;j2<vN.cluster[j].myInfo.myVp.nItems;j2++)
		{
			auxInf =  DistEuclideanDstd(vN.cluster[j].myInfo.myVp.V[j2],vN.cluster[j].myInfo.mean);
			sumInf+= auxInf;
			auxT =  DistEuclideanDstd(vN.cluster[j].myInfo.myVp.V[j2],center);
			sumT+= auxT;
		}
		sumFinal+= sumInf;
	}

	I = pow((sumT*auxMaxD/(sumFinal*NC)),sizeN);
	free(center.V);
	return I;
}

double IndexVp(VectorPoints* vN,int n)
{
	double I=0;
	int NC,nItemsN=0;
	int i,k,k2,j,j2,j3,j5,sizeN = vN[0].V[0].Size;
	Point center;
	double sumInf=0,sumFinal=0,sumT=0,auxT=0,auxInf=0;
	Point media[n];

	center.Size = sizeN;
	center.V = (double*)calloc(center.Size,sizeof(double));

	NC = n;


	for(i=0;i<NC;i++)
	{
		nItemsN+= vN[i].nItems;
		media[i].Size =  sizeN;
		media[i].V = (double*)calloc(media[i].Size,sizeof(double));
		for(j5=0;j5<vN[i].nItems;j5++)
			for(k=0;k<sizeN;k++)
			{
				media[i].V[k] += vN[i].V[j5].V[k];
				center.V[k] += vN[i].V[j5].V[k];
			}
	}

	for(k2=0;k2<sizeN;k2++)
	{
		center.V[k2]= center.V[k2]/nItemsN;
		for(i=0;i<NC;i++)
			media[i].V[k2]=media[i].V[k2]/vN[i].nItems;
	}

	int flag=0;
	double maxD=0,auxMaxD=0;

	for(j=0;j<NC;j++)
	{
		cout<<"A"<<endl;
		for(j3=j+1;j3<NC;j3++)
		{
			//                      printf("grupo%d=%lf \t contra %d=%lf\n",j,media[j].V[0],j3,media[j3].V[0]);
			//                    printf("grupo%d=%lf \t contra %d=%lf\n",j,media[j].V[1],j3,media[j3].V[1]);

			maxD = DistEuclideanDstd(media[j],media[j3]);
			cout<<"Z : "<<maxD<<endl;
			if(maxD > auxMaxD || flag==0)
			{
				auxMaxD = maxD;
				flag=1;
			}
		}

		sumInf=0;
		for(j2=0;j2<vN[j].nItems;j2++)
		{
			auxInf =  DistEuclideanDstd(vN[j].V[j2],media[j]);
			sumInf+= auxInf;
			auxT =  DistEuclideanDstd(vN[j].V[j2],center);
			sumT+= auxT;
		}
		sumFinal+= sumInf;
	}

	I = pow((sumT*auxMaxD/(sumFinal*NC)),sizeN);

	return I;
}


/*double finallyVal(VectorNodos vN)
{
    VectorPoints myValidations;

    CreateVectorPoints(5,)



}*/

/*void normalizarUp(VectorPoints *P)
{
    double auxMin = P[0];
    double auxMax = P[0];
    int i;

    for(i=0;i<P->nItems;i++)
    {
        if(P->V.V[i]>auxMax)
            auxMax = P->V.V[i];
        if(P->V.V[i]<auxMin)
            auxMin = P->V.V[i];
    }

    for(i=0;i<P->nItems;i++)
        P->V.V[i] = (P->V.V[i] - auxMin)/(auxMax -auxMin);
}*/
}// end namespace
