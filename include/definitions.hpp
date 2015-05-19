/*
 * definitions.h
 *
 *  Created on: Aug 5, 2014
 *      Author: crono21
 */

#ifndef DEFINITIONS_H_
#define DEFINITIONS_H_

#define PI 3.14159265
#define Ancho 640
#define Resol 1.5

#include "stdio.h"
#include "../../KLT/include/boost_incl.hpp"

namespace cluster{

using namespace std;
using namespace boost_incl;

typedef struct SolutionsTag
{
    int typeTree;  //  0=SL  1=CL  2 Mean  3=Mdian  4=minArea
    double*  evals;  //  Ya sean 5 � 2
}Solutions;

typedef struct PositionTag
{
   int row;
   int col;
}Position;


typedef struct PointTag
{
	short int Size;
	double* V;
}Point;

typedef struct VectorPointsTag
{
	int nItems;
	short int Size;
	Point* V;
	Point min;
	Point Max;
}VectorPoints;

typedef struct VectorDistTag
{
   int nItems;
   double *Dist;
}VectorDist;

typedef struct InfoTag
{
  int level;                 // Nivel en el arbol orden ascendente
  short int Size;			 // Tamaño de los datos
  Point median;              // Mediana
  Point mean;                // Media
  Point myPoint;             // punto(coordenadas que representan al nodo)
  VectorPoints myVp;         // puntos inicales que contiene el grupo
  VectorPoints myVpVelocity; // puntos inicales que contiene el grupo
  int * tags;			     // Valor de posicion en la matriz inicial de los puntos en myVp y myVpVelocity empieza en 0
  VectorDist myVdist;        // Vector de las distanciasde desde su punto  respecto al punto de los demas nodos
}Info;

typedef struct NodoTag
{
	Info myInfo;                  //  Informacion del nodo
	struct NodoTag* rNodo;            //   nodo derecho
	struct NodoTag* lNodo;            //   nodo izquierdo
	struct NodoTag* fNodo;            //   nodo padre
	double  **Ref;             //   posiciones de las cajas envolventes [coordenada][min, median, max]
	int flagNFA;                   //   Bandera que indica si es o no significativo
	int ready;                     //   Bandera que indica si el nodo ya ha sido evaludo
	int name;                      //   Etiqueta n�mero de elemento
	double NFA;                    //   NFA
	double NFAg;                   //   NFAg
	double Area;                   //   Area del rect�ngulo
	double Prob;                   //   Probabilidad
	int Dir;                       //   Nodo izquierdo = 0 , nodo derecho = 0
	double* Dstd;                  //   Desviaci�nes est�ndar del grupo
	//    VectorPoints velocity;
}Nodo, * TreeBinary;

typedef struct  VectorTreesTag
{
    int nItems;                   // N�mero de elementos de mi vector de nodos
    int Init;                     // numero inicial de arboles
    TreeBinary* Tree;             // Puntero de arboles
    Position *Pos;                // Posiciones para ordenar las distancias entre los nodos
    int flagNorm;                 //  Bandera de la normalizaci�n
//     double Norms[4];                //  puntero de valores de inicializacion:  para 4 dim:  width,height,vel,tetha
    int dimConfig;                //  posisci�n del elemento variable  en este caso la velocidad  = 2
    double minV,maxV;             //   Velocidad min y velocidad maxima
    double   NFA;                //  Velocidad del traker
    int size;                     //  tama�o de cada punto
    double **myRef;          //  vector de referencias  (Ndim  20)
    double *Dstd;            //  Desviaci�n est�ndar del conjunto de puntos de entrada para cada atributo
}VectorTrees;

typedef struct VectorNodosTag
{
	int   nCluster;
	VectorTrees *Tree;
	TreeBinary cluster;

}VectorNodos;

}

#endif /* DEFINITIONS_H_ */
