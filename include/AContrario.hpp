/*
 * AContrario.h
 *
 *  Created on: Aug 5, 2014
 *      Author: crono21
 */

#ifndef ACONTRARIO_H_
#define ACONTRARIO_H_

#include "definitions.hpp"
#include "evaluationMethods.hpp"
#include "RMSstd.hpp"
#include "BinaryTreeFunctions.hpp"
#include "NFAFunctions.hpp"
#include <iostream>

namespace cluster {

class AContrario {
public:
	AContrario();
	virtual ~AContrario();

	int clusteringTrackingPoints(mat & Xp, int nC, int nR,double sigma=0.2, int maxz = 0, double NFA = 0.001f); //Sigma for gaussian density probability in NFA
	void writeTrackingPointstoFile(mat & Xp, const char* name);


private:
	void freeVectorPoints(VectorPoints *vector);
	void printClustertoScreen(VectorNodos & cluster);
	void saveClusterstoMatrix(VectorNodos & cluster, mat & xp);
	VectorNodos moduloAgrupamientoVel(VectorPoints *myV,double (*f)(Nodo N1,Nodo N2),double NFA,double sigma,VectorPoints* Vel);
	int newTecnica(VectorPoints* myV,double NFA,double sigma,VectorPoints* Vel, VectorNodos *& clusters);
	void  DataConversion(mat & Xp, VectorPoints ** points,  VectorPoints ** velocities);

	int FIRST;
	int maxx;									//Reference boxes sizes
	int maxy;
	int maxz;
};

} /* namespace cluster */

#endif /* ACONTRARIO_H_ */
