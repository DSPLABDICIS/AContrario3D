/*
 * BinaryTreeFunctions.hpp
 *
 *  Created on: Jul 29, 2014
 *      Author: dsplab
 */

#ifndef BINARYTREEFUNCTIONS_HPP_
#define BINARYTREEFUNCTIONS_HPP_

#include "definitions.hpp"
#include "mainSLFunctions.hpp"

namespace cluster {

void CreateBinaryTree(VectorTrees *myTree,double (*f)(Nodo N1,Nodo N2));
VectorTrees* InitTreeVelocity(VectorPoints* Vp,VectorPoints* Vveloc);
void freeVectorList(VectorNodos* List,int size);
void freeVectorNode(VectorNodos* Ntree);
void freeTreeBinary(TreeBinary tree);

void InitMinMax(VectorTrees* myTree);




} /* namespace cluster */

#endif /* BINARYTREEFUNCTIONS_HPP_ */
