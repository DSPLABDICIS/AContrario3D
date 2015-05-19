/*
 * evaluationMethods.h
 *
 *  Created on: Jul 29, 2014
 *      Author: dsplab
 */

#ifndef EVALUATIONMETHODS_H_
#define EVALUATIONMETHODS_H_

#include "definitions.hpp"

namespace cluster {

typedef int (*FuncEvaluations)(Solutions* P,int nLinkages,int nMetrics);

int finallyEvaluation(Solutions* P,int nLinkage,int nmetrics,FuncEvaluations Feval);
int finallyEvaluationToOne(Solutions* P,int nLinkage,FuncEvaluations Feval);

int CLVpromedio(Solutions* P,int nLinkages,int nMetrics);
int CLVgeometrico(Solutions* P,int nLinkages,int nMetrics);
int CLVarmonica(Solutions* P,int nLinkages,int nMetrics);
int CLVstd(Solutions* P,int nLinkages,int nMetrics);





} /* namespace cluster */

#endif /* EVALUATIONMETHODS_H_ */
