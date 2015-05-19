/*
 * RMSstd.h
 *
 *  Created on: Jul 29, 2014
 *      Author: dsplab
 */



#ifndef RMSSTD_H_
#define RMSSTD_H_

#include "definitions.hpp"
#include "mainSLFunctions.hpp"
#include <iostream>
#include <iomanip>

namespace cluster{

typedef double (*FuncMetric)(VectorNodos vN);

double Index(VectorNodos vN);
double Bouldin(VectorNodos vN);
double Dunns(VectorNodos vN);
double Silhouette(VectorNodos vN);
double calinski(VectorNodos vN);
double RMSstd( VectorNodos vN);
double Rsquared(VectorNodos vN);

VectorNodos getSolutionGroups(VectorTrees& myTree);
double obtenerMedia(VectorPoints myV,int k);
double* obtenerDesviacionEstandar(VectorPoints myV);
double NormalHasting(double sigma, double media, double x);

} //end namespace
#endif /* RMSSTD_H_ */
