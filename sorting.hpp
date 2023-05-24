#ifndef __SORTING_H
#define __SORTING_H

#include "iostream"
#include "list"
#include "Eigen/Eigen"
#include "map"

using namespace std;
using namespace Eigen;

namespace SortLibrary {

  template<typename T>
  void Merge(vector<T>& v,
             const unsigned int& sx,
             const unsigned int& cx,
             const unsigned int& dx){
      unsigned int i,j,k;
      i=sx; j=cx+1;  k=0;
      vector<T> b;
      b=v;
      while((i<=cx)&& (j<=dx))
      {
          if(v[i]<=v[j])
          {

              b[k]=v[i];
              i=i+1;
          }
          else
              {
               b[k]=v[j];
               j=j+1;
              }
       k=k+1;
       }
      for (; i<=cx; i=i+1, k=k+1)
          b[k]=v[i];
      for (; j<=dx; j=j+1, k=k+1)
          b[k]=v[j];
      for ( i=sx; i<=dx;i=i+1)
          v[i]=b[i-sx];

  }

  template<typename T>
  void MergeSort(vector<T>& v,
                 const unsigned int& sx,
                 const unsigned int& dx){
      unsigned int cx;
      if(sx<dx)
      {
          cx=(sx+dx)/2;
          MergeSort( v, sx, cx);
          MergeSort( v, cx +1, dx);
          Merge(v,sx,cx,dx);

      }

  }

}

#endif // __SORTING_H
