/*
 * NFAFunctions.cpp
 *
 *  Created on: Jul 29, 2014
 *      Author: dsplab
 */

#include "../include/NFAFunctions.hpp"

namespace cluster {

double getCombinatoriaT(int N, int n1, int n2)
{
	int k,count=1,flag=0;
	double sum=1.0;

	for(k=N ; k>N-n1 ; k--)
	{
		sum *= k;
		if(count > n2 && flag==0)
		{
			count=1;
			flag = 1;
		}
		sum /= count;
		count++;
	}
	return sum;
}

double TrinomialAcumulada(int N,int n1,int n2,double p1,double p2)
{
	int i,j;
	double trinomialAcum=0,cTaux=0,cT=0,Prob=0,difP=0,Ppdif=0,Ppdifaux=0,Ppi=0,Ppj=0;

	difP = 1-p1-p2;
	cTaux = getCombinatoriaT(N,n1,n2);
	Ppi = pow(p1,n1);
	Ppj = pow(p2,n2);
	Ppdifaux = pow(difP,N-n1-n2);

	for(i = n1;i <= N-n2 ; i++)
	{
		cT = cTaux;
		Ppdif = Ppdifaux;
		for(j = n2 ; j <= N-i ; j++)
		{
			Prob = Ppi*Ppj*Ppdif;
			trinomialAcum+=cT*Prob;
			//  printf("%d,%d:Trinomial:%1.30lf\n",i,j,trinomialAcum);
			cT*=1.0*(N-i-j)/(j+1);     //  actualiza cTaux,Ppi,Ppdif respecto a j
			Ppj*=p2;
			Ppdif*=difP;
		}
		cTaux*=1.0*(N-n2-i)/(i+1);        //  actualiza cTaux,Ppi,Ppdif respecto a i
		Ppi*=p1;
		Ppdifaux*=difP;
	}
	//     printf("Trinomial:%1.30lf\n",trinomialAcum);
	return trinomialAcum;
}

double getNFAg(TreeBinary myNodo,TreeBinary nL,TreeBinary nR,VectorTrees myTree)                // GetNFA de grupos
{
	double cte =0;
	int H=20,N;

	N = myTree.Init;
	cte = N*N*pow(H,myTree.size);  // con 4 elem todo: modificaciones para 6 elementos?
	// cte = N*pow(H,myTree.size); // con 2 elem
	//cte*=cte*cte*cte;
	//  printf("N:%d,nL:%d,nR:%d,pl:%1.20lf,pr:%1.20lf\n",N, nL->myInfo.myVp.nItems , nR->myInfo.myVp.nItems , nL->Prob,nR->Prob);
	myNodo->NFAg = cte*TrinomialAcumulada(N, nL->myInfo.myVp.nItems -1 , nR->myInfo.myVp.nItems-1 , nL->Prob,nR->Prob);
	return myNodo->NFAg;
}

int pointsInIntersection(double **box,VectorPoints myVp)    // Devuelve e numero de puntos que se ecuentra dentro de esa regi�n
{
	int i,j,count,P=0;

	for(i=0;i<myVp.nItems;i++)
	{
		count=0 ;
		for(j=0;j<myVp.V[0].Size;j++)
			if(myVp.V[i].V[j] < box[j][1] && myVp.V[i].V[j] > box[j][0])
				count++;
		if(count==4)
			P++;
	}
	return P;
}

void  checkIntersection(TreeBinary nL,TreeBinary nR)
{
	int i,j,count =0,nIr=0,nIl=0;
	double AInte=1.0;
	double **Intersection;
	double ProbAux=0;

	Intersection = (double**)malloc(4*sizeof(double*));

	for(i=0;i<4;i++)
	{
		Intersection[i]=(double*)malloc(3*sizeof(double));
		if(((nL->Ref[i][0] < nR->Ref[i][0]) &&  (nL->Ref[i][1] > nR->Ref[i][0]))|| ((nL->Ref[i][0] < nR->Ref[i][0]) &&  (nL->Ref[i][1] > nR->Ref[i][0])))
			count++;
	}

	if(count==4)             //  Si count=4 entonces existe intersecci�n
	{
		for(j=0;j<4;j++)
		{
			if(nL->Ref[j][0] < nR->Ref[j][0])
				Intersection[j][0]= nR->Ref[j][0];
			else
				Intersection[j][0]= nL->Ref[j][0];

			if(nL->Ref[j][1] > nR->Ref[j][1])
				Intersection[j][1]= nR->Ref[j][1];
			else
				Intersection[j][1]= nL->Ref[j][1];
			AInte *= (Intersection[j][1] - Intersection[j][0]);
		}


		nIl = pointsInIntersection(Intersection,nL->myInfo.myVp);
		nIr = pointsInIntersection(Intersection,nR->myInfo.myVp);

		if(nIl+nIr > 0)
		{
			ProbAux = AInte/(pow(Ancho/Resol,nL->myInfo.myPoint.Size));;           // seg�n la funci�n de probabilidad
			nL->Prob -=  ProbAux;
			nR->Prob -=  ProbAux;
			nR->myInfo.myVp.nItems -= nIr;
			nL->myInfo.myVp.nItems -= nIl;
			return;
		}
		else
		{
			nR->Prob = nR->myInfo.myVp.nItems/((nR->Prob)/(nR->myInfo.myVp.nItems) - AInte);
			nL->Prob = nL->myInfo.myVp.nItems/((nL->Prob)/(nL->myInfo.myVp.nItems) - AInte);
			return;
		}
	}
	else
		return;

}



void GetMinReferenceBarycenter(TreeBinary myNodo,VectorTrees myTree)   // Usando el baricentro obtniene "la caja" o referecias
{
	int i,j,k,m,item=0,w,flag=0;
	double dT,Daux=0;
	double LimS, limi;
	int nRef=20;

	//-------------------------------------------------------------------------------
	// Asignar memoria a la  caja contenedora obtenida de las referencias
	myNodo->Ref = (double**)malloc(myNodo->myInfo.Size*sizeof(double*));
	for(i=0;i<myNodo->myInfo.Size;i++)
		myNodo->Ref[i] = (double*)malloc(3*sizeof(double));   // 0->  lim inferior, 1-> lim superior,  2-> size
	//-------------------------------------------------------------------------------

	//-------------------------------------------------------------------------------
	//  obtener el punto mas cercano al baricentro
	for(m=0; m < myNodo->myInfo.myVp.nItems; m++)
	{
		//               printf("%lf,%lf\t",myNodo->myInfo.myVp.V[m], myNodo->myInfo.median);
		dT = DistEuclideanDstd(myNodo->myInfo.myVp.V[m], myNodo->myInfo.median);
		//             printf("%d,%lf\t",m,dT);
		if( dT < Daux || m==0)
		{
			Daux = dT;
			item = m;
		}
	}
	//-------------------------------------------------------------------------------

	//         LimS =Ancho/Resol;
	//         printf("#de elem:%d\n",myNodo->myInfo.myVp.nItems);
	//PrintPoint(myNodo->myInfo.median);
	//PrintPoint(myNodo->myInfo.myVp.V[item]);


	//     for(h=0;h<myNodo->myInfo.myPoint.Size;h++)
	//       myNodo->Ref[h] = (double*)malloc(3*sizeof(double));
	for(k=0; k<myNodo->myInfo.Size; k++)               // para cada dimension
	{
		Daux=0;
		w=0;

		LimS = myTree.myRef[k][nRef-1];

		for(j=0; j<myNodo->myInfo.myVp.nItems; j++)                   //Para cada punto del vevtor de puntos
			if(j!=item)                                                 // excepto el elemento centro(el + cercano al baricentro)
			{
				dT =(myNodo->myInfo.myVp.V[j].V[k] - myNodo->myInfo.myVp.V[item].V[k]); //Calcula la distancia punto menos Pcentral
				dT *= 2.0;                                                               // Tama�o de la caja 2veces la distancia

				//-------------------------------------------------------------------------------
				//evaluar que las dimensiones no caigan fuerade la ventana total
				if(dT/2.0 >0 && (myNodo->myInfo.myVp.V[item].V[k] - dT/2.0)<0)
					dT = myNodo->myInfo.myVp.V[j].V[k];

				if(dT/2.0 < 0 && (myNodo->myInfo.myVp.V[item].V[k] - dT/2.0) > LimS)
					dT = LimS - myNodo->myInfo.myVp.V[j].V[k];
				//-------------------------------------------------------------------------------
				if(dT<0)
					dT*=-1.0;

				if( dT > Daux)
				{
					Daux = dT;
					flag=j;
				}

			}
		//           printf("Daux:%lf\t",Daux);
		w=-1;
		if(Daux>LimS)
			myNodo->Ref[k][2]=LimS;
		else
			do
			{
				w++;
				myNodo->Ref[k][2] = myTree.myRef[k][w];
			}
			while(Daux > myTree.myRef[k][w] && w <  nRef);

		//          printf("D:%lf\t",myNodo->Ref[k][2]);


		dT =(myNodo->myInfo.myVp.V[flag].V[k] - myNodo->myInfo.myVp.V[item].V[k]);
		// printf("dt:%lf",dT);
		if(dT>0)
		{
			if((myNodo->myInfo.myVp.V[item].V[k] - dT)>0)
			{
				if(myNodo->myInfo.myVp.V[item].V[k] + myNodo->Ref[k][2]/2.0 >LimS)
				{
					myNodo->Ref[k][0] = LimS-myNodo->Ref[k][2];
					myNodo->Ref[k][1] = LimS-0.00001;
				}
				else{
					myNodo->Ref[k][0] = myNodo->myInfo.myVp.V[item].V[k] - myNodo->Ref[k][2]/2.0;
					myNodo->Ref[k][1] = myNodo->myInfo.myVp.V[item].V[k] + myNodo->Ref[k][2]/2.0;
				}
			}
			else
			{
				myNodo->Ref[k][0] = 0;
				myNodo->Ref[k][1] = myNodo->Ref[k][2];
			}
		}
		else
		{
			if((myNodo->myInfo.myVp.V[item].V[k] - dT)< LimS)
			{
				if(myNodo->myInfo.myVp.V[item].V[k] - myNodo->Ref[k][2]/2.0 < 0)
				{
					myNodo->Ref[k][0] = 0.0;
					myNodo->Ref[k][1] = myNodo->Ref[k][2];
				}
				else
				{
					myNodo->Ref[k][1] = myNodo->myInfo.myVp.V[item].V[k] + myNodo->Ref[k][2]/2.0;
					myNodo->Ref[k][0] = myNodo->myInfo.myVp.V[item].V[k] - myNodo->Ref[k][2]/2.0;
				}
			}
			else
			{
				myNodo->Ref[k][0] = LimS-myNodo->Ref[k][2];;
				myNodo->Ref[k][1] = LimS-0.00001;
			}
		}
		//                             printf("Min:%lf\t",myNodo->Ref[k][0]);
		//                                       printf("max:%lf\n",myNodo->Ref[k][1]);

	}
	//               printf("%lf\n",myNodo->Ref.V[k]);
}

double Combinatoria(int N,int n)
{
	int i;
	double sum=1.0;

	for(i=n+1; i<N+1; i++)
		sum *= ((1.0*i)/(i-n));
	return sum;
}

double binomialAcum(int n,int x,double Prob)
{
	double bAcum=0,c,difProb,probP,probN;
	int i;

	difProb=1-Prob;
	probP=pow(Prob,x);
	probN=pow(difProb,n-x);
	c = Combinatoria(n,x);
	for(i=x;i<=n;i++)
	{
		bAcum += c*probP*probN;
		c*=1.0*(n-i)/(i+1);
		probP*=Prob;
		probN*=(1.0/difProb);

	}
	//        printf("\nbAT:%1.25lf",bAcum);
	return bAcum;
}

double getNFA(TreeBinary myNodo, VectorTrees myTree, double sigma)                 // Obtiene el valor de la NFA del Nodo de entrada
{
	double Binomial,cte,Atotal=1.0,tempp;
	int h,nReg=20;

	if(myNodo->ready==1)
		return myNodo->NFA;

	//cout<<myNodo->name<<endl;
	myNodo->Prob =1;                                   // Iniciar probabilidad
	GetMinReferenceBarycenter(myNodo,myTree);         //  Obtener las referencias todo: verificar si negativos influencian el resultado
	myNodo->Area=1.0;                                 // Iniciar �rea

	for(h=0;h<myTree.size;h++)                        //  Obtener Area de Universo  y area de la caja
	{
		//printf("A:%lf\n",myNodo->Ref[h][2]);
		myNodo->Area*= myNodo->Ref[h][2];

		Atotal*= myTree.myRef[h][nReg-1];

	}
	//printf("At=%lf\tAr=%lf\n",Atotal,myNodo->Area);
	myNodo->Prob = (myNodo->Area/(Atotal));
	tempp = myNodo->Prob;

	//   Para la velocidad
	//---------------------------------------------------------------------------------------------------------------

	double P,Vmean;
	double media[]={0.0,0.0,0.0}; //Modelo de camara Movil
	double *S;

	//printf("Vx%lf,Vx%lf\n",myNodo->myInfo.myVpVelocity.V[FIRST].V[0],myNodo->myInfo.myVpVelocity.V[FIRST].V[1]);


	//VXmean=obtenerMedia(myNodo->myInfo.myVpVelocity,0);
	//VYmean=obtenerMedia(myNodo->myInfo.myVpVelocity,1);
	//Px = NormalHasting(sigma,0.0,-mediaX+VXmean);
	//Py = NormalHasting(sigma,0.0,-mediaY+VYmean);
	//myNodo->Prob *= Px*Py;

	for(h=0;h<myTree.size;h++)                        //  Obtener media de velocidades
	{
		Vmean = obtenerMedia(myNodo->myInfo.myVpVelocity,h);
		P = NormalHasting(sigma,0.0,-media[h]+Vmean);
		myNodo->Prob *= P;
	}

	//S=obtenerDesviacionEstandar(myNodo->myInfo.myVpVelocity);    //Todo: VErificar si se requiere calcular

	// printf("MediaX:%lf\tMediaY:%lf\tSY:%lf\n",VXmean,VYmean,S[1]) ;

	//   printf("gpo:%d\tprob NormalX:%1.30lf\tprob NormalY:%1.30lf\n",myNodo->name,Px,Py) ;

	//       myNodo->Prob = myNodo->myInfo.myVp.nItems/myNodo->Area; // Doctora Dora


	//-------------------------------------------------------------------------------------------------------------


	//       myNodo->Prob = (myNodo->Area/(Atotal));///(myNodo->myInfo.myVp.nItems*myNodo->myInfo.myVp.nItems);                 //  probabilidada igual a la raz�n de areas
	//       myNodo->Prob = (myNodo->Area/(Atotal));                 //  probabilidada igual a la raz�n de areas
	if(myNodo->Prob==1)
		Binomial=1;
	else
		Binomial=binomialAcum(myTree.Init,myNodo->myInfo.myVp.nItems-1,myNodo->Prob);

	//       cte=pow(20,myTree.size)*myTree.Init*myTree.Init;  / con 4 elem
	cte=pow(nReg,myTree.size)*myTree.Init*myTree.Init*myTree.Init;  // con 2 elem todo: cambio de la función para tres elemntos?
	myNodo->NFA = cte * Binomial;
	//              printf("NFA:%1.50lf\n",myNodo->NFA);
	//printf("Prob:%lf\t,Bin:%1.40lf \t,cte:%lf \tNFA:%1.40lf\n",myNodo->Prob,Binomial,cte,myNodo->NFA);
	myNodo->ready=1;

	//cout <<std::scientific<<setw(5)<<"Ppos  "<<myNodo->Prob << "\tPtotal: "<< myNodo->Prob <<"\tN: "<< myNodo->myInfo.myVp.nItems <<"\tBinA: "<< Binomial <<"\tNFAO: "<< myNodo->NFA <<"   "<<myNodo->name<<endl;

	return myNodo->NFA;
}

int getNpoints(TreeBinary myNodo)
{
	int m = myNodo->myInfo.myVp.nItems;
	return    m;
}

int Conditions(TreeBinary myNodo,VectorTrees myTree,double NFAc, double sigma)   //  compara la NFA con la NFA de grupos; 0 si NFAg < NFA  y 1 si no
{
	double H;
	double NFAl=0.0,NFAr=0.0;

	NFAl= getNFA(myNodo->lNodo,myTree, sigma);
	NFAr= getNFA(myNodo->rNodo,myTree, sigma);
	getNFAg(myNodo,myNodo->lNodo,myNodo->rNodo,myTree);
	//      printf("NFAg:%lf\t",myNodo->NFAg );
	if(myNodo->NFAg < NFAc)
		return 1;
	return 0;
}

void getGralSelection(TreeBinary myNodo,VectorTrees& myTree,TreeBinary nodoEval, double sigma)  // Proceso de selecci�n de nodos con NFA <= 1
{           // Evaluation 0= no ha sido evaluado por lo tanto debe compararse con una nueva NFA, =1  si ya fue evaluado compararse con la NFA recibida
	int count=0,flag;
	TreeBinary aux = NULL;
	double NFAc=0;

	aux =  myNodo;
	//printf("Estoy evaluando los hijos del nodo: %d contra el nodo: %d\n",myNodo->name,nodoEval->name);
	NFAc = nodoEval->NFA;
	if(aux==nodoEval)
	{
		NFAc=getNFA(aux,myTree, sigma);
		if(aux->NFA > myTree.NFA)
		{
			//printf("Mi NFA de: %d es mayor a e\n",aux->name);
			if(getNpoints(aux->lNodo) > 1)
				getGralSelection(aux->lNodo,myTree,aux->lNodo, sigma);
			if(getNpoints(aux->rNodo) > 1)
				getGralSelection(aux->rNodo,myTree,aux->rNodo, sigma);
			return;
		}
		//printf("Mi NFA es de: %lf es menor que el umbral \n",aux->NFA);
		myNodo->flagNFA=1;
	}
	do
	{
		if(getNpoints(aux->lNodo) > 1  && getNpoints(aux->rNodo) > 1)   //  Verifica si los hijos son aptos para
		{
			flag = Conditions(aux,myTree,NFAc, sigma);    //  me devuelve 1 si la NFAg < NFA  ; 0 si es mayor o igual
			if(flag == 1)
			{
				//                        printf(" los hijos de: %d tienen una NFAg menor que la NFA de: %d \n",aux->name,nodoEval->name);
				nodoEval->flagNFA=0;
				getGralSelection(aux->lNodo,myTree,aux->lNodo,sigma);
				getGralSelection(aux->rNodo,myTree,aux->rNodo,sigma);
				//}
				aux = aux->fNodo;
				while( count>0)
				{
					getGralSelection(aux->rNodo,myTree,nodoEval,sigma);
					aux = aux->fNodo;
					count--;
				}
				return;
			}
			// printf(" los hijos de: %d tienen una NFAg mayor que la NFA de: %d=%1.10lf \n",aux->name,nodoEval->name,nodoEval->NFA);
		}
		else
		{
			//                printf(" los hijos de: %d no cumplen condici�n de hijos\n",aux->name);
			if(getNpoints(aux->lNodo) > 1)
				getGralSelection(aux->lNodo,myTree,nodoEval,sigma);
			if(getNpoints(aux->rNodo) > 1)
				getGralSelection(aux->rNodo,myTree,nodoEval,sigma);

			aux = aux->fNodo;
			while( count>0)
			{
				getGralSelection(aux->rNodo,myTree,nodoEval,sigma);
				aux = aux->fNodo;
				count--;
			}
			return;
		}
		aux = aux->lNodo;
		count++;
	}
	while(getNpoints(aux->lNodo) > 1);
	return;
}

void   Select_NFA(TreeBinary myNodo,VectorTrees myTree,double sigma)  // Proceso de selecci�n de nodos con NFA <= 1
{           // Evaluation 0= no ha sido evaluado por lo tanto debe compararse con una nueva NFA, =1  si ya fue evaluado compararse con la NFA recibida
	int count=0,flag;
	TreeBinary aux;
	double NFAc=0;

	aux = (TreeBinary)malloc(sizeof(Nodo));
	aux =  myNodo;

	getNFA(aux,myTree,sigma);
	if(aux->NFA < myTree.NFA)
	{
		myNodo->flagNFA=1;
		return;
	}
	else
	{
		if(getNpoints(aux->lNodo) > 1)
			Select_NFA(aux->lNodo,myTree,sigma);
		if(getNpoints(aux->rNodo) > 1)
			Select_NFA(aux->rNodo,myTree,sigma);
		return;
	}
}

void   SelectNodos2(TreeBinary myNodo,VectorTrees myTree)    // REcorrido del �rbol mostrando algo dadoo funci�n G
{
	int count=0;
	TreeBinary aux;

	aux = (TreeBinary)malloc(sizeof(Nodo));
	aux =  myNodo;

	if(aux->lNodo == NULL)
		return;

	if((aux->myInfo.myVp.Max.V[2]-aux->myInfo.myVp.min.V[2])< myTree.NFA*myTree.Dstd[2]  && (aux->myInfo.myVp.Max.V[3]-aux->myInfo.myVp.min.V[3]) < myTree.NFA*myTree.Dstd[3])
	{
		aux->flagNFA=1;
		return;
	}

	do
	{
		aux = aux->lNodo;
		count++;
		if((aux->myInfo.myVp.Max.V[2]-aux->myInfo.myVp.min.V[2])< myTree.NFA*myTree.Dstd[2]  && (aux->myInfo.myVp.Max.V[3]-aux->myInfo.myVp.min.V[3]) < myTree.NFA*myTree.Dstd[3])
		{
			aux->flagNFA=1;
			aux = aux->fNodo;
			while( count>0)
			{
				SelectNodos2(aux->rNodo,myTree);
				aux = aux->fNodo;
				count--;
			}
			return;
		}
		if( aux->myInfo.myVp.nItems == 1)  // aqui va la condicion
		{
			aux = aux->fNodo;
			while( count>0)
			{
				SelectNodos2(aux->rNodo,myTree);
				aux = aux->fNodo;
				count--;
			}
			return;
		}
	}
	while(aux->lNodo != NULL);
	free(aux);
	return;
}

void CreateReferences(VectorTrees* myTree,int maxx,int maxy, int maxz)                    //  Obtiene las referencias de cada variable(Normalizadas o no Normalizadas)
{
	int j,k,val[myTree->size];
	//    double tol=0.0001;
	int nReg=20;

	//    val[0]=100;      //  Tama�o m�ximo en X              %640
	//  val[1]=100;      //  Tama�o m�ximo en y             %480
	/*    val[2]=100;                //  Tama�o m�ximo en Vx   %400
    val[3]=400;              //  Tama�o m�ximo en Vy    %400*/

	//      Para  im�genes reales
	int D=1;
	float _3dval = 1.0;
	if(myTree->size==3){
		_3dval = 0.1f;
		D = 10;}

	val[0]=D*maxx;      //  Tama�o m�ximo en X              %640 todo: Measures??
	val[1]=D*maxy;      //  Tama�o m�ximo en y             %480*/

	if(myTree->size == 3)
		val[2]=D*maxz;

	//     val = Ancho/Resol;
	//                       val=6*myTree->Dstd[k];

	myTree->myRef=(double**)malloc(myTree->size*sizeof(double*));
	for(k=0 ; k<myTree->size;k++)
	{
		myTree->myRef[k]=(double*)malloc(nReg*sizeof(double));
		//    val=6*myTree->Dstd[k];
		//                   printf("val:%d\n",val);
		for(j=0 ; j<nReg;j++)
		{
			myTree->myRef[k][j] = pow(exp((log(val[k]))/(nReg-1)),j)*_3dval;
			//                          if(k==0)
			//                       printf("%lf\t",myTree->myRef[k][j]);
		}
		//             printf("\n\n");
	}
}

} /* namespace cluster */
