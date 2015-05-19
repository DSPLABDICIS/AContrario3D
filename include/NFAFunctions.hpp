/*
 * NFAFunctions.hpp
 *
 *  Created on: Jul 29, 2014
 *      Author: dsplab
 */

#ifndef NFAFUNCTIONS_HPP_
#define NFAFUNCTIONS_HPP_

#include "definitions.hpp"
#include "mainSLFunctions.hpp"
#include "RMSstd.hpp"

namespace cluster {

void CreateReferences(VectorTrees* myTree,int maxx, int maxy, int maxz); //Optional depth value if 3D
void getGralSelection(TreeBinary myNodo,VectorTrees& myTree,TreeBinary nodoEval,double sigma);

} /* namespace cluster */

#endif /* NFAFUNCTIONS_HPP_ */
