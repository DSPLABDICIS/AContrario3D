/*
 * BinaryTreeFunctions.cpp
 *
 *  Created on: Jul 29, 2014
 *      Author: dsplab
 */

#include "../include/BinaryTreeFunctions.hpp"

namespace cluster {

void freeTreeBinary(TreeBinary tree){
	static int count =0;

	// Recursive running until the last children
	if(tree->rNodo != NULL)
		freeTreeBinary(tree->rNodo);
	if(tree->lNodo != NULL)
		freeTreeBinary(tree->lNodo);

	// delete the memory when it is a last level node
	int i;
	free(tree->myInfo.myVp.V);
	free(tree->myInfo.myVpVelocity.V);
	free(tree->myInfo.tags);
	free(tree->myInfo.myVdist.Dist);
	free(tree->myInfo.myVp.min.V);
	free(tree->myInfo.myVp.Max.V);
	if(tree->Ref!=NULL){
		for(i=0;i<tree->myInfo.Size;i++)
			free(tree->Ref[i]);
		free(tree->Ref);
	}


	// If it is a group node extra memory to delete
	if(tree->myInfo.level>1){
		free(tree->myInfo.myPoint.V);
		free(tree->myInfo.mean.V);
		free(tree->myInfo.median.V);
	}
	//free the Tree itself
	free(tree);
	count++;
}

void freeVectorNode(VectorNodos* Ntree){

	//Free Vector Trees
	if(Ntree->Tree!=NULL){
	if(Ntree->Tree->Tree[0] !=NULL)
		freeTreeBinary(Ntree->Tree->Tree[0]);
	free(Ntree->Tree->Tree);
	for(int i=0; i<Ntree->Tree->size;i++)
		free(Ntree->Tree->myRef[i]);
	free(Ntree->Tree->myRef);
	free(Ntree->Tree->Dstd);}
}


void freeVectorList(VectorNodos* List,int size){

	for(int i=0;i<size;i++){
			freeVectorNode(List+i);
			free(List[i].cluster);
			free(List[i].Tree);
	}
}



VectorDist CreateVectorDist(int N)                                                //  Vector para guardar las distancias entre los elementos
{
	VectorDist myV;
	int i;

	myV.nItems = N;
	myV.Dist =(double*)malloc(myV.nItems*sizeof(double));
	for(i=0;i<myV.nItems;i++ )
		myV.Dist[i]=0;
	return myV;
}

void CreateVectorDist2(VectorDist &myV, int N)                                                //  Vector para guardar las distancias entre los elementos
{
	int i;

	myV.nItems 	= N;
	myV.Dist 	=(double*)malloc(myV.nItems*sizeof(double));
	for(i=0;i<myV.nItems;i++)
		myV.Dist[i]=0;
}

TreeBinary InitNodo(int Npuntos, int Nnodos, int size)
{
	TreeBinary nNodo;
	int static Name=1;

	nNodo  = (TreeBinary)malloc(sizeof(Nodo));                                  // memoria para nodo


	nNodo->myInfo.level						= Npuntos;
	nNodo->myInfo.Size						= size;
	nNodo->myInfo.myVp.nItems 				= Npuntos;
	nNodo->myInfo.myVp.Size					= size;
	nNodo->myInfo.myVpVelocity.nItems 		= Npuntos;  							// con velocidad
	nNodo->myInfo.myVpVelocity.Size			= size;
	nNodo->myInfo.myVp.V 					= (Point*)malloc(nNodo->myInfo.myVp.nItems*sizeof(Point));
	nNodo->myInfo.myVpVelocity.V 			= (Point*)malloc(nNodo->myInfo.myVpVelocity.nItems*sizeof(Point));   // Con velocidad
	nNodo->myInfo.tags 						= (int*)malloc(nNodo->myInfo.myVpVelocity.nItems*sizeof(int));

	//Promptly going to change
	nNodo->myInfo.median.Size 				= size;
	nNodo->myInfo.mean.Size 				= size;
	//****************************

	CreateVectorDist2(nNodo->myInfo.myVdist,Nnodos);


	//----------------------------------------------------
	nNodo->NFA 		= 1000.0;
	nNodo->NFAg 	= 0.0;
	nNodo->flagNFA 	= 0;
	nNodo->ready	= 0;
	nNodo->name		= Name++;
	nNodo->Ref      = NULL;

	nNodo->rNodo = NULL;
	nNodo->lNodo = NULL;
	nNodo->fNodo = NULL;

	return nNodo;
}

void CalculateBarycenterMean(Nodo* nNodo,Nodo* A,Nodo* B)                          // Calcular Baricentro = a la media de sus puntos mNodo=(mA*nA+mB*nB)/nNodo
{
	int i;

	nNodo->myInfo.mean.V=(double*)malloc(nNodo->myInfo.myPoint.Size*sizeof(double));
	for(i=0;i<nNodo->myInfo.myPoint.Size;i++)
		nNodo->myInfo.mean.V[i] = ((A->myInfo.mean.V[i]) * (A->myInfo.myVp.nItems) + (B->myInfo.mean.V[i]) * (B->myInfo.myVp.nItems))/(nNodo->myInfo.myVp.nItems);
}

void CalculateBarycenterMedian(Nodo* nNodo,Nodo* A,Nodo* B)                         // Punto max y min de su nube de puntos
{
	int k;

	nNodo->myInfo.myVp.min.V=(double*)malloc(nNodo->myInfo.myPoint.Size*sizeof(double));
	nNodo->myInfo.myVp.Max.V=(double*)malloc(nNodo->myInfo.myPoint.Size*sizeof(double));
	nNodo->myInfo.median.V=(double*)malloc(nNodo->myInfo.myPoint.Size*sizeof(double));

	for(k=0;k<nNodo->myInfo.myPoint.Size;k++)
	{
		if(A->myInfo.myVp.min.V[k] < B->myInfo.myVp.min.V[k])
			nNodo->myInfo.myVp.min.V[k] = A->myInfo.myVp.min.V[k];
		else
			nNodo->myInfo.myVp.min.V[k] = B->myInfo.myVp.min.V[k];

		if(A->myInfo.myVp.Max.V[k] > B->myInfo.myVp.Max.V[k])
			nNodo->myInfo.myVp.Max.V[k] = A->myInfo.myVp.Max.V[k];
		else
			nNodo->myInfo.myVp.Max.V[k] = B->myInfo.myVp.Max.V[k];
		nNodo->myInfo.median.V[k] = (nNodo->myInfo.myVp.Max.V[k] + nNodo->myInfo.myVp.min.V[k])/2.0;   // Obtener el punto medio del conjunto (max-min)/2
	}
}

void LoadInfo(Nodo* nNodo,Nodo* A,Nodo* B)                                        // llenar la nueva informacion a partir de los dos nuevos nodos
{
	int i,j,ind;

	nNodo->myInfo.level = A->myInfo.level > B->myInfo.level ? A->myInfo.level +1 : B->myInfo.level +1;
	CalculateBarycenterMean(nNodo,A,B);
	CalculateBarycenterMedian(nNodo,A,B);

	for(i=0,ind =0;i<A->myInfo.myVp.nItems;i++)
	{
		nNodo->myInfo.myVp.V[i]= A->myInfo.myVp.V[i];
		nNodo->myInfo.myVpVelocity.V[i]= A->myInfo.myVpVelocity.V[i];
		nNodo->myInfo.tags[ind++] = A->myInfo.tags[i];
	}
	for(j=A->myInfo.myVp.nItems;j < (B->myInfo.myVp.nItems + A->myInfo.myVp.nItems);j++)
	{
		nNodo->myInfo.myVp.V[j]= B->myInfo.myVp.V[j-A->myInfo.myVp.nItems];
		nNodo->myInfo.myVpVelocity.V[j]= B->myInfo.myVpVelocity.V[j-A->myInfo.myVpVelocity.nItems];
		nNodo->myInfo.tags[ind++] = B->myInfo.tags[j-A->myInfo.myVpVelocity.nItems];
	}
}

TreeBinary getNewTree(Nodo* R, Nodo* L,Point (*f)(Point r,Point l)) // Obtener nuevo Nodo a partir de otros dos mediante funci�n f*
{
	int i, ind;
	TreeBinary nTree;

	nTree = InitNodo(R->myInfo.myVp.nItems + L->myInfo.myVp.nItems,L->myInfo.myVdist.nItems,L->myInfo.myPoint.Size);
	nTree->myInfo.myPoint = (*f)(R->myInfo.myPoint,L->myInfo.myPoint);

	LoadInfo(nTree,R,L);
	R->Dir = 1;           //  define la variable Dir 0 para nodo derecho,  1 para nodo izquierdo
	L->Dir = 0;
	nTree->rNodo = R;
	nTree->lNodo = L;
	R->fNodo = nTree;
	L->fNodo = nTree;

	return nTree;
}

void CreateBinaryTree(VectorTrees *myTree,double (*f)(Nodo N1,Nodo N2))//
{
	Position P;
	int countName=myTree->Init,i;

	GetVectorDistDstd(myTree,(*f));   //  se obtienen las distancias de todo el conjunto
	do
	{
		P = GetMinDist(myTree);

		myTree->Tree[P.row] = getNewTree(myTree->Tree[P.col], myTree->Tree[P.row], FuncAverage); //con la informacion de los dos nodos hijos se inicializa el nuevo nodo y se va a la posici�n del primero de ellos
		myTree->Tree[P.row]->name =  countName;

		if(P.col != (myTree->nItems - 1))                                                               //  En caso de ser necesario, cambiar el segundo punto por el ultimo de la lista
			myTree->Tree[P.col] = myTree->Tree[myTree->nItems - 1];

		myTree->nItems--;
		countName++;
		if(myTree->nItems > 1)                                 // se realiza el proceso de calculo de las nuevas distancias para obtener los nuevos dos puntos
			ProcessDistDstd(myTree,P,(*f));
	}
	while(myTree->nItems > 1);
	return;
}

void InitMinMax(VectorTrees* myTree)
{
	Position Pos;
	int i,Tam,dim,element;

	for(element=0;element<myTree->nItems;element++)
	{
		myTree->Tree[element]->myInfo.myVp.min.V=(double*)malloc(myTree->size*sizeof(double));
		myTree->Tree[element]->myInfo.myVp.Max.V=(double*)malloc(myTree->size*sizeof(double));

		for(dim=0;dim<myTree->size;dim++)
		{
			myTree->Tree[element]->myInfo.myVp.Max.V[dim] = myTree->Tree[element]->myInfo.myPoint.V[dim];
			myTree->Tree[element]->myInfo.myVp.min.V[dim] = myTree->Tree[element]->myInfo.myPoint.V[dim];
		}
	}
}

VectorTrees InitTree(VectorPoints* Vp)   //  Iniciar el �rbol
{

	VectorTrees myTrees;
	int element,dim, eInit = 1;


	myTrees.Tree = (TreeBinary*)malloc(Vp->nItems*sizeof(TreeBinary));
	myTrees.nItems = Vp->nItems;
	myTrees.Init = Vp->nItems;
	myTrees.size = Vp->V[0].Size;


	for(element=0 ; element<myTrees.nItems ; element++)
	{
		myTrees.Tree[element] = InitNodo(eInit,myTrees.nItems,Vp->V[element].Size);  // Inicializar cada Nodo
		myTrees.Tree[element]->myInfo.myPoint = Vp->V[element];                    // Paso los elementos
		myTrees.Tree[element]->myInfo.level=eInit;                               // Nivel de cada nodo =1
		myTrees.Tree[element]->myInfo.myVp.nItems=eInit;                          // Numero de elementos=1
		myTrees.Tree[element]->myInfo.median.V=(double*)malloc(Vp->V[element].Size*sizeof(double));
		myTrees.Tree[element]->myInfo.median = myTrees.Tree[element]->myInfo.myPoint; // Baricentro es el mismo punto
		myTrees.Tree[element]->myInfo.mean.V=(double*)malloc(Vp->V[element].Size*sizeof(double));
		myTrees.Tree[element]->myInfo.mean = myTrees.Tree[element]->myInfo.myPoint; // Baricentro es el mismo punto
		myTrees.Tree[element]->myInfo.myVp.V[0] = myTrees.Tree[element]->myInfo.myPoint; // Mi vector de puntos solo contien el mismo punto.
	}
	return myTrees;
}

VectorTrees* InitTreeVelocity(VectorPoints* Vp,VectorPoints* Vveloc)   //  Iniciar el �rbol
{

	VectorTrees *myTrees = (VectorTrees*) malloc(sizeof(VectorTrees));
	int element,eInit = 1;

	myTrees->Tree 	= (TreeBinary*)malloc(Vp->nItems*sizeof(TreeBinary));
	myTrees->nItems = Vp->nItems;
	myTrees->Init 	= Vp->nItems;
	myTrees->size 	= Vp->Size;

	//cout<<myTrees.nItems<<endl;
	for(element=0 ; element<myTrees->nItems ; element++)
	{
		myTrees->Tree[element] 								= InitNodo(eInit,myTrees->nItems,Vp->Size);// Inicializar cada Nodo
		myTrees->Tree[element]->myInfo.myPoint 				= Vp->V[element];							// Paso los elementos
		myTrees->Tree[element]->myInfo.median 				= Vp->V[element]; 							// Baricentro es el mismo punto
		myTrees->Tree[element]->myInfo.mean					= Vp->V[element]; 							// Baricentro es el mismo punto
		myTrees->Tree[element]->myInfo.tags[0] 				= element; 									// Tag del punto en la matriz original
		myTrees->Tree[element]->myInfo.myVp.V[0] 			= Vp->V[element]; 							// Mi vector de puntos solo contien el mismo punto.
		myTrees->Tree[element]->myInfo.myVpVelocity.V[0] 	= Vveloc->V[element]; 						// Mi vector de puntos solo contien el mismo punto.
	}
	return myTrees;
}

void ExplorerTree3(TreeBinary myNodo,FILE *fp,void (*G)(TreeBinary B,FILE *A))    // REcorrido del �rbol escribir inf en archivo .txt
{
	int count=0;
	TreeBinary aux;

	aux = (TreeBinary)malloc(sizeof(Nodo));
	aux =  myNodo;

	if(aux->lNodo == NULL)
		return;

	do
	{
		(*G)(aux,fp);                          //  funcion para imprimir
		aux = aux->lNodo;
		count++;
		if( aux->myInfo.myVp.nItems == 1)  // aqui va la condicion
		{
			aux = aux->fNodo;
			while( count>0)
			{
				ExplorerTree3(aux->rNodo,fp,(*G));
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

void ExplorerTree2(TreeBinary myNodo,int *val,void (*G)(TreeBinary B,int *A))    // REcorrido del �rbol mostrando algo dadoo funci�n G
{
	int count=0;
	TreeBinary aux;

	aux = (TreeBinary)malloc(sizeof(Nodo));
	aux =  myNodo;

	if(aux->lNodo == NULL)
		return;
	do
	{
		(*G)(aux,val);                          //  funcion para imprimir
		aux = aux->lNodo;
		count++;
		if( aux->myInfo.myVp.nItems == 1)  // aqui va la condicion
		{
			aux = aux->fNodo;
			while( count>0)
			{
				ExplorerTree2(aux->rNodo,val,(*G));
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





} /* namespace cluster */
