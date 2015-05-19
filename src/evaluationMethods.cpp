/*
 * evaluationMethods.cpp
 *
 *  Created on: Jul 29, 2014
 *      Author: dsplab
 */

#include "../include/evaluationMethods.hpp"

namespace cluster {

int finallyEvaluation(Solutions* P,int nLinkage, int nMetrics,FuncEvaluations Feval)
{
    double auxMin[5] ;
    double auxMax[5] ;
    int i,j,k,sol=0;

       for(i=0 ; i< nMetrics ; i++)   //   //  para cada metrica de validacion  //  obtener m�nimos y m�ximos para la normalizaci�n.
         {
             auxMax[i] = P[0].evals[i];           //  se le asigna  el valor del primero de ellos
             auxMin[i] = P[0].evals[i];           //  se le asigna  el valor del primero de ellos
             for(j=0 ; j< nLinkage ; j++)   //   para cada linkage
                  {
                     if(P[j].evals[i] > auxMax[i])
                           auxMax[i] = P[j].evals[i];
                     if(P[j].evals[i] < auxMin[i])
                            auxMin[i] = P[j].evals[i];
                   }
         }
    for(k=0 ; k<nLinkage ; k++)   // Normalizaci�n
        {
            P[k].evals[0] = (P[k].evals[0] - auxMin[0])*0.9/(auxMax[0]- auxMin[0]) + 0.1;    //  Index Max
            P[k].evals[1] = 1 - (P[k].evals[1] - auxMin[1])*0.9/(auxMax[1]- auxMin[1]);    //  Bouldin Min
            P[k].evals[2] = (P[k].evals[2] - auxMin[2])*0.9/(auxMax[2]- auxMin[2]) + 0.1;    //  Dunns  Max
            P[k].evals[3] = (P[k].evals[3] - auxMin[3])*0.9/(auxMax[3] - auxMin[3])+0.1; //silhouete   Max
            P[k].evals[4] = (P[k].evals[4] - auxMin[4])*0.9/(auxMax[4] - auxMin[4])+0.1;  // calinski  Max
        }
    sol = Feval(P,nLinkage,5);
    return P[sol].typeTree;
}

int finallyEvaluationToOne(Solutions* P,int nLinkage,FuncEvaluations Feval)
{
    double auxMin[5] ;
    double auxMax[5] ;
    int i,j,k,nValidation=1,sol=0;


       for(i=0 ; i< nValidation ; i++)   //  para cada metrica de validacion  //  obtener m�nimos y m�ximos para la normalizaci�n.
         {
             auxMax[i] = P[0].evals[i];           //  se le asigna  el valor del primero de ellos
             auxMin[i] = P[0].evals[i];           //  se le asigna  el valor del primero de ellos
              for(j=0 ; j< nLinkage ; j++)   //   para cada linkage
                  {
                     if(P[j].evals[i] > auxMax[i])
                           auxMax[i] = P[j].evals[i];
                     if(P[j].evals[i] < auxMin[i])
                            auxMin[i] = P[j].evals[i];
                   }
         }
    //  Para minimos: RMS, Rsou Bouldin
    for(k=0 ; k<nLinkage ; k++)
            P[k].evals[0] = 1 - (P[k].evals[0] - auxMin[0])*0.9/(auxMax[0]- auxMin[0]);    //  RMSstd  Min

    sol = Feval(P,nLinkage,nValidation);

return P[sol].typeTree;
}

int CLVpromedio(Solutions* P,int nLinkages,int nMetrics)
 {
    int i,j,sol;
    double prom,auxProm;

    for(j=0;j<nLinkages;j++)
      {
         prom=0;
         for(i=0;i<nMetrics;i++)
            prom+=P[j].evals[i];
            if(prom > auxProm || j==0)
              {
                 auxProm = prom;
                 sol=j;
             }
       }
    return sol;
 }

int CLVgeometrico(Solutions* P,int nLinkages,int nMetrics)
{
    int i,j,sol;
    double geom,auxGeom;

    for(j=0;j<nLinkages;j++)
      {
         geom=1;
         for(i=0;i<nMetrics;i++)
            geom*=P[j].evals[i];
         if(geom > auxGeom || j==0)
            {
                auxGeom = geom;
                sol=j;
             }
       }
    return sol;
}

int CLVarmonica(Solutions* P,int nLinkages,int nMetrics)
{
    int i,j,sol,elem=nMetrics,nM;  //  valor cualquiera;
    double armo,auxArmo=0.0,valMin;

    for(j=0;j<nLinkages;j++)
      {
         armo=0.0;
         if(nMetrics==5)
           {
            elem = 0;
            valMin = P[j].evals[0];
            for(nM=1;nM<nMetrics;nM++)
                if(P[j].evals[nM] < valMin)
                    {
                        valMin = P[j].evals[nM];
                        elem = nM;
                     }
             }
         for(i=0;i<nMetrics;i++)
             if(i!=elem)
                 armo+=1.0/P[j].evals[i];

         armo = 1.0/armo;
         if(armo > auxArmo || j==0)
             {
                 auxArmo = armo;
                 sol=j;
             }
       }
    return sol;
}

int CLVstd(Solutions* P,int nLinkages,int nMetrics)
{
    int i,j,i2,sol=0;
    double mean,vari,Dstd,auxDstd;

   for(j=0;j<nLinkages;j++)
      {
         mean=0;
         vari=0;
         for(i=0;i<nMetrics;i++)
             mean += P[j].evals[i];
             mean/=nMetrics;
         for(i2=0;i2<nMetrics;i2++)
             vari += (P[j].evals[i2]-mean)*(P[j].evals[i2]-mean);
             Dstd = sqrt(vari);
             if(Dstd < auxDstd || j==0)
             {
                 auxDstd = Dstd;
                 sol=j;
             }
         cout<<"Linkage("<<j<<"):"<<Dstd<<endl;
       }
   cout<<endl;
    return sol;
}


} /* namespace cluster */
