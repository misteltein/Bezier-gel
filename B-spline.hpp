#ifndef BS_H
#define BS_H
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cfloat>
//#include "Vector3Dedit.hpp"
/*
A B-spline is a generalization of the Besier curve.
Let a vector known as the knot vector be defined
   T = (t_0, t_1, ..., t_m),
where T is a nondecreasing sequence with
   t_i in [0,1],
and define control points
   P = (P_0, P_1, ..., P_n).
Define the degree as p := m-n-1.
The knots t_(p+1),...,t_(m-p-1) are called internal knots.
Define the basis function as

   N_i,0(t) = [ 1 if t_i <= t < t_(i+1) and t_i<t_(i+1),
              [ 0 otherwise

                  t - t_i                       t_(i+j+1) - t
   N_i,j(t) = --------------- N_i,(j-1)(t) + ---------------------N_(i+1),(j-1)(t) ,
              t_(i+j) - t_i                  t_(i+j+1) - t_(i+1)

where j = 1, 2, ..., p.
Then the curve defined by
   C(t) = Sum_(i=0)^(n) P_i N_i,p (t)

*/

template <typename VECTOR>
class B_spline
{
   public:
      B_spline(){}
      B_spline(const std::vector<VECTOR>& ps, const int& degree=3);
      void echo(int color=0)const;
      VECTOR C(const double& t)const;
   private:
      std::vector<VECTOR> P;//points
      std::vector<double> T;//knots
      double N(const int& i, const int& j, const double& t)const;
      int n;
      int m;
      int p;//degree
};

template <typename VECTOR>
B_spline<VECTOR>::B_spline(const std::vector<VECTOR>& ps, const int& degree)
{
   P = ps;
   n = ((int)P.size())-1;
   p = (degree>(n))?(n):(degree);//if p=n, Bézier curve.
   m = p+n+1;
   {//make T
      for(int i=0;i<=p;++i)
      {
         T.push_back(0.0);
      }
      const double delta_t  = 1.0/((m-p)-p);
      double       t_       = delta_t;
      for(int i=p+1;i<=m-p-1;++i)
      {
         T.push_back(t_);t_+=delta_t;
      }
      for(int i=m-p;i<=m;++i)
      {
         T.push_back(1.0);
      }
   }
}
template <typename VECTOR>
double B_spline<VECTOR>::N(const int& i, const int& j, const double& t)const
{
   if(0==j)
   {
      return ((T[i]<=t)&&(t<T[i+1]))?1.0:0.0;
   }
   if(j>p){std::cout<<"fail"<<std::endl;}
   double unit_A = ((t-T[i])/(T[i+j]-T[i]))*N(i,j-1,t);
   double unit_B = ((T[i+j+1]-t)/(T[i+j+1]-T[i+1]))*N(i+1,j-1,t);
   if(isnan(unit_A)){unit_A=0.0;}
   if(isnan(unit_B)){unit_B=0.0;}
   return unit_A+unit_B;
}

template <typename VECTOR>
void B_spline<VECTOR>::echo(int color)const
{
   constexpr double delta = 0.001;
   int max = 1.0/delta;
   for(int i=0;i<=max;++i)
   std::cout<<C(i*delta)<<" "<<color<<std::endl; 
}

template <typename VECTOR>
VECTOR B_spline<VECTOR>::C(const double& t)const
{
   if(std::abs(t)<10*DBL_EPSILON){return P[0];}
   if(std::abs(t-1.0)<10*DBL_EPSILON){return P[n];}

   VECTOR res;
   for(int i=0;i<=n;++i)
   {
      res += P[i]*N(i,p,t);
   }
   return res;
}

//int main()
//{
//   std::vector<Vector3D> pts;
//   for(int i=0;i<=10;++i)
//   {
//      pts.push_back(Vector3D(i,i*i,i));
//   }
//   std::for_each(pts.begin(),pts.end(),[](auto d){std::cout<<d<<std::endl;});
//   B_spline<Vector3D> bs(pts);
//   bs.echo();
//}

#endif
