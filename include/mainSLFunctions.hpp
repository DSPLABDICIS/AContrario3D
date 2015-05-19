/*
 * mainSLFunctions.hpp
 *
 *  Created on: Jul 29, 2014
 *      Author: dsplab
 */

#ifndef MAINSLFUNCTIONS_HPP_
#define MAINSLFUNCTIONS_HPP_

#include "definitions.hpp"

namespace cluster{

typedef double (*FuncLinkage)(Nodo N1,Nodo N2);

double singleLinkage(Nodo N1,Nodo N2);
double completeLinkage(Nodo N1,Nodo N2);
double medianLinkage(Nodo N1,Nodo N2);
double meanLinkage(Nodo N1,Nodo N2);
double minAreaLinkage(Nodo N1,Nodo N2);

void GetVectorDistDstd(VectorTrees* myV,double (*f)(Nodo N1,Nodo N2));
Position GetMinDist(VectorTrees* myV);
double obtenerMedia(VectorPoints myV,int k);
double* obtenerDesviacionEstandar(VectorPoints myV);
Point FuncAverage(Point r,Point l);
void ProcessDistDstd(VectorTrees* myV,Position myP,double (*f)(Nodo N1,Nodo N2));
double DistEuclideanDstd(Point p1,Point p2);


}


#endif /* MAINSLFUNCTIONS_HPP_ */
