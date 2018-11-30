#ifndef CONSTEXPR_MATH_H
#define CONSTEXPR_MATH_H

namespace cexpr_math{
template <typename T>
constexpr T sqr(T s)
{
   return s*s;
}

template <typename T>
constexpr T sqrt(T s)
{
   T x = s / 2.0;
   T prev = 0.0;
   while(x!=prev)
   {
      prev=x;
      x=(x+s/x)/2.0;
   }
   return x;
}

template<typename T>
constexpr T sqrt_imp(T a)
{
   T x,h,g;
   x=DBL_EPSILON;
   g=1.0;
   h=sqr(1.0-a*x*x);
   while(h<sqr(g))
   {
      x+=x*(h*(8.0+h*(6.0+5.0*h))/16.0);
      g=h;
      h=sqr(1.0-a*x*x);
   }
   return x;
}

template <typename T>
constexpr double max(const T& v1, const T& v2){return (v1<v2)?v2:v1;}

template <typename T>
constexpr T pow(const T& a, const int& b)
{
   T result=a;
   for(int i=1;i<b;++i){result*=a;}
   return result;
}

}

#endif
