#include <iostream>
#include <vector>
#include <list>
#include <random>
#include <float.h>
#include "constexpr_math.hpp"
#include "Vector3Dedit.hpp"
#include "B-spline.hpp"

//parameter
//random dispacement is restricted on sphere of radius r=[0.5*LIMITS_A:LIMITS_A]
constexpr double LIMITS_A = 1.0;
//in making a random branch, number of the joint point is restricted [0.5*LIMITS_JOINT_POINT_NUM:LIMITS_JOINT_POINT_NUM]
constexpr int    LIMITS_JOINT_POINT_NUM = 20;
//the number of branch is restricted [LIMITS_BRANCH_NUM]
constexpr int    LIMITS_BRANCH_NUM = 50;
//in making smooth branch, the number of split is SPLIT_FACTOR * (num. of branch's jointpoint)
constexpr double SPLIT_FACTOR = 1.0;

std::mt19937 thread_mt;
Vector3D random_d()
{
   std::uniform_real_distribution<double>    dist_d(-LIMITS_A,LIMITS_A);
   Vector3D d;
   do{
      d.x=dist_d(thread_mt);
      d.y=dist_d(thread_mt);
      d.z=dist_d(thread_mt);
   }while
   (
      ((d.x)*(d.x)+(d.y)*(d.y)+(d.z)*(d.z))>(cexpr_math::sqr(LIMITS_A))||
      ((d.x)*(d.x)+(d.y)*(d.y)+(d.z)*(d.z))<0.5*(cexpr_math::sqr(LIMITS_A))
   );
   return d;
}

int main()
{
   std::vector<Vector3D> points;
   std::vector<std::list<int>>        branches;
   std::uniform_int_distribution<int> dist_sp(LIMITS_JOINT_POINT_NUM/2,LIMITS_JOINT_POINT_NUM);
   points.push_back(Vector3D (0.0,0.0,0.0));
   int s=0;
   while(s<LIMITS_BRANCH_NUM)
   {
      std::list<int> bch;
      for(int i=0,sp_n=dist_sp(thread_mt);i<sp_n;++i)
      {
         if(i==0)
         {
            std::uniform_int_distribution<int> dist_t(0,int(points.size())-1);
            const int picked = [&]()
            {
               bool end_f=false;
               int result;
               do{
                  result = dist_t(thread_mt);
                  end_f=true;
                  for(int b=0,sbze=branches.size();b<sbze;++b)
                  {
                     const std::list<int>& bch_ = branches[b];
                     if(result==bch_.front()||result==bch_.back())
                     {
                        end_f=false;break;
                     }
                  }
               }while(!end_f);
               return result;
            }();
            const Vector3D p = points[picked]+random_d();
            bch.push_back(picked);
            points.push_back(p);
            for(int c=0,scze=branches.size();c<scze;++c)
            {
               bool exist_f = false;
               std::list<int>& bch = branches[c];
               std::list<int> fst;
               std::list<int> snd;
               for(auto it=bch.begin();it!=bch.end();++it)
               {
                  if(!exist_f){fst.push_back(*it);};
                  if(picked==*it)
                  {
                     exist_f=true;          
                  }
                  if(exist_f){snd.push_back(*it);};
               }
               if(exist_f)
               {
                  bch=fst;
                  branches.push_back(snd);
                  break;
               }
            }
         }
         else
         {
            const Vector3D p = points.back()+random_d();
            points.push_back(p);
         }
         bch.push_back(points.size()-1);
      }
      branches.push_back(bch);
      ++s;
   }
   
   for(int i=0,size=branches.size();i<size;++i)
   {
      const std::list<int>& bch = branches[i];
      std::vector<Vector3D> control_points;
      for(auto it=bch.begin();it!=bch.end();++it)
      {
         control_points.push_back(points[(*it)]);
      }
      //smooth it!
      B_spline<Vector3D> tmp(control_points,(int)control_points.size());//<- means bezier

      for(int t = 0,tmax=(int(SPLIT_FACTOR*((int)control_points.size())));t<=tmax;++t)
      {
         std::cout<<tmp.C(t*(1.0/tmax))<<" "<<i<<std::endl;
      }
      std::cout<<std::endl<<std::endl;
   }
   return EXIT_SUCCESS;
}
