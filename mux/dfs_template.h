#include<iostream>
#include<fstream>
#include<cstdlib>
#include<random>
#include<ctime>
#include<algorithm>
#include<iterator>
#include<math.h>
#include<ctime>
#include<stdlib.h>
#include<iomanip>
#include<sstream>
#include<string>
#include<new>
#include <sys/stat.h>
#include <dirent.h>

#define WHITE 1
#define GRAY 2
#define NIL 0
#define BLACK 3


//Paper: Sarika Jalan and Priodyuti Pradhan, Localization of multilayer networks by optimized single-layer rewiring, Phys. Rev. E 97, 042314 (2018)
// Written: Priodyuti Pradhan, Complex System Lab,IIT Indore
// contact: priodyutipradhan@gmail.com  


using namespace std;

template<class X> void Display_Network(unsigned long int n, X **g);
template<class X> unsigned int dfs(X **g,unsigned int n);
template <class X> int DFS_Visit(X **g,unsigned int u, unsigned int n, unsigned int *color, unsigned int *pi);
template<class X> unsigned int Max_Degree(X **g,unsigned long int N);
template<class X> unsigned int Max_Degree_Index(X **g, unsigned int N);
template<class X>void Degree_sequence(unsigned int N, X **G, unsigned int Deg[]);
void Normalize_Deg_Vector(unsigned int *Deg, double *Ndeg, unsigned int N);
template<class X>void Matrix_Copy(X **a,X **a_New,unsigned long int N);
template<class X> int Allocate_2Darray(X **g,unsigned long int N, unsigned long int M);
template<class X>void Free_2Darray(X **g, unsigned long int N);
template<class X> int Network_File(X **g, unsigned int N, string path, unsigned int File_Num);
template <class X> unsigned long int Access_Network(X **A, string path);
int create_directory(string path);
int file_copy(string src, string dest);
int plot_data(string src_path);

template<class X> int Allocate_2Darray(X **g,unsigned long int N,unsigned long int M)
{
 for(unsigned long int i = 0; i < N; i++)
 { 
    try{   
       g[i] = new X[M];        //Create N x M dimensional matrix
    } catch (bad_alloc xa) {
   cout << "Allocation Failure\n";
   return 1;
  }  
  }
//::::::::::: Initialize to zeros::::::::::::::
 for(unsigned long int i=0;i<N;i++)
   for(unsigned long int j=0;j<M;j++)
      g[i][j] = 0; 
}

template<class X>void Free_2Darray(X **g, unsigned long int N)
{
  for(unsigned long int i = 0; i < N; i++)
     delete [] g[i];
  delete [] g;
}

template<class X> void Display_Network(unsigned long int n,X **g)
{
  for(unsigned long int i=0;i<n;i++)
  {
    for(unsigned long int j=0;j<n;j++)
      cout<<g[i][j]<<' ';
   cout<<endl;
  }
}

int create_directory(string path)
{
    DIR *pDir;
    string a = "mkdir -p ";
    string in = a + path;
  
    pDir = opendir(path.c_str()); 
    if (pDir == NULL) {
      std::cout<<"Creating Directory..."<<std::endl;
      //const int dir= system("mkdir -p //home//priodyuti//Desktop//myfolder1");
      const int dir = system(in.c_str());
      if (dir< 0)
          return 0;
   }                                                                                                
   else{ 
     std::cout<<"Directory already exists..."<<std::endl;
     return 1;
   } 
   closedir (pDir);
}


 
int file_copy(string src, string dest)
{  
  cout << "File copying..."<<endl;
  string a = "cp ";
  src += " ";                                                             
  string path = a+src+dest;
  //const int dir= system("cp //home//priodyuti//Desktop//test.txt //home//priodyuti//Desktop//t.txt");
  const int dir = system(path.c_str());

  return 1;
}   

int plot_data(string src_path)
{  
  string a = "gnuplot -p -e \"plot '";
  string c = "'\"";
  string gnu_path = a+src_path+c;
  //system("gnuplot -p -e \"plot '/home/priodyuti/Desktop/data.txt'\"");
  system(gnu_path.c_str());

  return 1;
}

template<class X>void Degree_sequence(unsigned int N, X **G, unsigned int Deg[])
{
   for(int i=0;i<N;i++){
    int counter=0;
    for(int j=0;j<N;j++){
      if(G[i][j]==1) 
        counter++;
    }
    Deg[i] = counter;
   }
}

void Normalize_Deg_Vector(unsigned int *Deg,double *Ndeg, unsigned int N)
{ // Normalize in norm 2 
  double sum = 0.0;
  for(int i=0;i<N;i++)
  {
    sum = sum+(Deg[i]*Deg[i]);
  }
  sum = sqrt(sum);
  for(int i=0;i<N;i++)
  { 
   Ndeg[i] = Deg[i]/sum;
  }
}

template<class X> unsigned int Max_Degree(X **g, unsigned long int N)
{
 unsigned long int Max_degree=0, max_deg_index;
 
 for(unsigned long int i=0;i<N;i++)
 {
   unsigned long int Degree=0;
   for(unsigned long int j=0;j<N;j++)
   {
     if(g[i][j]==1)
       Degree=Degree+1;
   }
   if(Degree>Max_degree)
   {
     Max_degree = Degree;
     max_deg_index = i+1;
   } 
 }
 return Max_degree;
}

template<class X> unsigned int Max_Degree_Index(X **g, unsigned int N)
{
 int Max_degree=0,max_deg_index;
 
 for(int i=0;i<N;i++)
 {
  unsigned int Degree=0;
  for(int j=0;j<N;j++)
  {
    if(g[i][j]==1)
      Degree = Degree+1;
    if(Degree > Max_degree)
    {
       Max_degree = Degree;
       max_deg_index = i+1;
    }   
  }
}
return max_deg_index;
}

template <class X> unsigned long int Access_Network(X **A, string path)
{
  unsigned long int counter=0;
  ifstream infile(path);
  if(!infile)
  {
      cout<<"File not opened to Read"<<endl;
      return 1;
  }
  unsigned long int a, b;
  while (infile >> a >> b)
  {
    A[a-1][b-1]=1;
    A[b-1][a-1]=1;
    counter++;
  }
  infile.close();
  return counter;
}

template<class X> int Network_File(X **g, unsigned int N, string path, unsigned int File_Num)
{
  string filename = path + to_string(File_Num) +".txt";
  ofstream outFile(filename);
  if(!outFile)
  {
   cout<<"File not opened for Writing"<<endl;
   return 1; 
  }
  for(int i=0; i<N-1; i++)
   for(int j=i+1; j<N; j++)
     if(g[i][j] == 1)
       outFile<<i+1<<' '<<j+1<<endl;
  outFile.close();
}

template<class X>void Matrix_Copy(X **a_source, X **a_target, unsigned long int N)
{
 for(unsigned long int i=0; i<N; i++)
   for(unsigned long int j=0; j<N; j++)
       a_target[i][j] = a_source[i][j];  
}

//::::::::::::::::;;This below DFS code will check for connected ness::::::::::::::::::
template<class X> int DFS_Visit(X **g, unsigned int u, unsigned int n, unsigned int *color, unsigned int *pi)
{
  unsigned int *Adj;
  try{
      Adj = new unsigned int[n];
  } catch (bad_alloc xa) {
   cout << "Allocation Failure\n";
   return 1;
  }
  unsigned int i,j,v;
  color[u] = GRAY;
  for(i=0,j=0;i<n;i++)
  {
   if(g[u][i]==1)
    Adj[j++]=i;
  }
  for(i=0;i<j;i++)
  {
    v = Adj[i];
    
    if(color[v]==WHITE) 
    {
      pi[v]=u;
      DFS_Visit(g,v,n,color,pi);
    }
  } 
  color[u]=BLACK;
  delete [] Adj;
}

template<class X> unsigned int dfs(X **g, unsigned int n)
{
  unsigned int *color,*pi;
  color = new unsigned int[n];
  pi = new unsigned int[n];

  for(unsigned int i=0;i<n;i++)
  {
     color[i] = WHITE;
     pi[i] = NIL;
  }   
  cout<<"Checking for connectedness..."<<endl;
  DFS_Visit(g,0,n,color,pi);
  for(unsigned int i=0;i<n;i++)
   if(color[i] != BLACK) 
   {
     delete [] color;
     delete [] pi;
     return 0;
   }
 delete [] color;
 delete [] pi; 
 return n; 
}

