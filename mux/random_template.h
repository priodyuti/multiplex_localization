#include<iostream>
#include<fstream>
#include<cstdlib>
#include<random>
#include<ctime>
//#include<algorithm>
//#include<iterator>
#include<math.h>
//#include<ctime>
#include<stdlib.h>
#include<iomanip>
#include<sstream>
#include<string>

//Paper: Sarika Jalan and Priodyuti Pradhan, Localization of multilayer networks by optimized single-layer rewiring, Phys. Rev. E 97, 042314 (2018)

// Written: Priodyuti Pradhan, Complex System Lab,IIT Indore
// contact: priodyutipradhan@gmail.com  

class edges{
public:
 int x;
 int y;
};

template<class X>unsigned long int ER_Random(unsigned int N, X **g, double p);
double ran_number();
template <class X>int Edge_Rewire(X **a,X **a_New,unsigned int N, unsigned long int Edges, int EDGE[]);
unsigned long int Random_Index( unsigned long int k);


unsigned long int Random_Index(unsigned long int k)
{
 std::random_device rd;                                                  
 std::mt19937 gen(rd());
 std::uniform_int_distribution<> dis(0, k);
 unsigned long int Random_Index = dis(gen);
 return Random_Index;
}

double ran_number()
{
  std::random_device rd;                                                  
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0, 1);
  
  return dis(gen);
}

template<class X>unsigned long int ER_Random(unsigned int N, X **g, double p) 
{ 
 unsigned int nu;
 double x;
 while(1){
   unsigned long edges=0;
   for(int i=0;i<N-1;i++){
     for(int j=i+1;j<N;j++){
        x=ran_number();  
        //cout<<x<<endl; 
        if(x<p){
          g[i][j]=1;
          g[j][i]=1;
          edges++;
        }
        else{
         g[i][j]=0;
         g[j][i]=0;
        }
     }
     g[i][i]=0;
   }
   nu = dfs(g,N);// Check for connectedness using DFS
   if(nu==N)
     return edges;
  }
}

template <class X>int Edge_Rewire(X **a,X **a_New,unsigned int N,unsigned long int Edges, int EDGE[])
{
  unsigned long int T=((N*(N-1))/2)-Edges;
  unsigned long int k=0,l=0;
  edges *e1, *e2;
  e1 = new edges[Edges];   //Edges present
  e2 = new edges[T];       // Edges not present
  for(int i=0;i<N-1;i++){
    for(int j=i+1;j<N;j++){
      if(a[i][j]==1){
        e1[k].x=i;
        e1[k].y=j;
        k=k+1;
      }
      else{
        e2[l].x=i;
        e2[l].y=j;
        l=l+1;
       }
     }                               
  }

 //CHOOSE INTEGER FROM A RANGE WITH UNIFORM DISTRIBUTION
 while(1)
 {
   Matrix_Copy(a,a_New,N);
   unsigned long int Random_Index1 = Random_Index(k-1);
   a_New[e1[Random_Index1].x][e1[Random_Index1].y] = 0;   //Choose only one edge and remove it randomly
   a_New[e1[Random_Index1].y][e1[Random_Index1].x] = 0;
   EDGE[0] = e1[Random_Index1].x;
   EDGE[1] = e1[Random_Index1].y;
 
   unsigned long int Random_Index0 = Random_Index(l-1);
   a_New[e2[Random_Index0].x][e2[Random_Index0].y] = 1; //randomly select a pair of vertices and add an edge between them 
   a_New[e2[Random_Index0].y][e2[Random_Index0].x] = 1;
   EDGE[2] = e2[Random_Index0].x;
   EDGE[3] = e2[Random_Index0].y;

   int nu = dfs(a_New,N);// Check for connectedness using DFS
   if(nu == N)
   {
     delete [] e1; delete [] e2;
     return 1; 
   }
 }
}

