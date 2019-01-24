#include<iostream>
#include<fstream>
#include<cstdlib>
#include"dfs_template.h"
#include"Eigen_template.h"
#include"random_template.h"
#include<math.h>
#include<ctime>
#include<stdlib.h>
#include<iomanip>
#include<sstream>
#include<string>
//#include<algorithm>
//#include<iterator>
#include<chrono>
#include <new>

//Paper: Sarika Jalan and Priodyuti Pradhan, Localization of multilayer networks by optimized single-layer rewiring, Phys. Rev. E 97, 042314 (2018)

// Written: Priodyuti Pradhan, Complex System Lab,IIT Indore
// contact: priodyutipradhan@gmail.com  

using namespace std;
using namespace std::chrono;

class Mux;
class Monoplex;

int Rewiring(Mux ob_old, Mux ob_new, unsigned int t);
void layer_rewiring(bool** a_New, Monoplex *A_New[], Mux ob_old, Mux ob_new, unsigned long int Layer_Index);


class Monoplex{
protected:
  unsigned int n, k;
  bool** A;
  double p; 
public:
  Monoplex(unsigned int a, unsigned int b){ n = a; k = b; p = k/double(n);}
  Monoplex(unsigned int a){ n = a;}
  unsigned long int edges;
  void give_mat();
  bool** give_mat(unsigned int a,unsigned int b);
  void mat_init_random();
  void mat_display();
  bool access_layers(unsigned int p, unsigned int q);
  void find_ipr();
  void free_mat();
  void extract_mux_layer(Monoplex *X[],unsigned int l);
  int rewire_edges(bool** A_New, unsigned long int Edges); 
};

class Mux: public Monoplex { 
  unsigned int N;
public:
  unsigned int l;
  unsigned long int *Edges; 
  Mux(unsigned int a, unsigned int b, unsigned int c, unsigned int d): Monoplex(c, d)
  { 
     N = a; 
     l = b;
     Edges = new unsigned long int[l]; 
  }
  void init_mux();
  void display_mux();
  int construct_mux(Monoplex *X[]);
  void conn_inter_layer(unsigned int k);
  void conn_intra_layer(Monoplex *X[], unsigned int k);
  void insert_intra_layer(bool** a_new, unsigned int k);
  double find_ipr_mux();
  void free_mux();
  void copy_mux(Mux ob_old, Mux ob_new);
  void store_network(string path, unsigned int iter);
  friend int Rewiring(Mux ob_old, Mux ob_new, unsigned int t, unsigned int N1);
  friend void layer_rewiring(bool** a_New, Monoplex *A_New[], Mux ob_old, Mux ob_new, unsigned long int Layer_Index);
};

void Monoplex::give_mat() { A = new bool*[n]; Allocate_2Darray(A, n, n); }

void Monoplex::mat_init_random() { edges = ER_Random(n, A, p); }

bool Monoplex::access_layers(unsigned int p, unsigned int q) { return (A[p][q]); }

void Monoplex::free_mat() { Free_2Darray(A, n); }

void Monoplex::mat_display() { Display_Network(n, A); cout<<endl; }

int Monoplex::rewire_edges(bool** A_New, unsigned long int Edges)
{
  int *EDGE;
  try{  
  EDGE = new int[4]; // allocate space for an int
  } catch (bad_alloc xa) {
   cout << "Allocation Failure\n";
   return 1;
  }
  
  Edge_Rewire(A, A_New, n, Edges, EDGE);
  delete [] EDGE;
}

void Monoplex:: find_ipr()
{
  double *Ev = new double[n];
  double max_eig_val = Eigen11(A,n,Ev);
  cout<<"Max Eigen value of A: "<<max_eig_val<<endl;
  double ipr = IPR(Ev, n);
  cout<<"IPR of A: "<<ipr<<endl;
  unsigned int max_deg = Max_Degree(A, n);
  cout<<"Max Degree of A: "<<max_deg<<endl;
  delete [] Ev;
}

void Monoplex:: extract_mux_layer(Monoplex *X[],unsigned int l)
{    // k is the index of the matrix
 for(int k=1; k<=l; k++)
   for(int i = n*(k-1), p = 0; i<k*n; i++, p++)
     for(int j = n*(k-1), q = 0; j<k*n ; j++,q++)
          X[k-1]->A[p][q]  = A[i][j];
}

void Mux::init_mux() { A = new bool*[N]; Allocate_2Darray(A, N, N); }

void Mux::free_mux() { Free_2Darray(A, N); }

void Mux::display_mux() { Display_Network(N, A); cout<<endl; }

int Mux::construct_mux(Monoplex *X[])
{
  for(int k=1;k<l;k++)
    conn_inter_layer(k);

  for(int k=1;k<=l;k++)
    conn_intra_layer(X, k);

  for(int i=0;i<l;i++)
    Edges[i]= X[i]->edges;
     
  return 0;
}

double Mux:: find_ipr_mux()
{
  double *Ev = new double[N];
  double max_eig_val = Eigen11(A, N, Ev);
  double c = IPR(Ev, N);
  unsigned int max_deg = Max_Degree(A, N);
  //cout<<"IPR of mux, max_deg, max_eig_val: "<<c<<' '<<max_deg<<' '<<max_eig_val<<endl;
  delete [] Ev;
  return c; 
} 

void Mux::conn_inter_layer(unsigned int k)
{
  for(int i=(k-1)*n,j=k*n; i<k*n; i++,j++)
  {
     A[i][j] = 1;
     A[j][i] = 1;
  }
}

void Mux:: conn_intra_layer(Monoplex *X[], unsigned int k)
{    // k is the index of the matrix
 for(int i = n*(k-1), p = 0; i<k*n; i++, p++)
   for(int j = n*(k-1), q = 0; j<k*n ; j++,q++)
      A[i][j] = X[k-1]->access_layers(p,q);
}

void Mux:: insert_intra_layer(bool** a_new, unsigned int k)
{    // k is the index of the matrix
  for(int i = n*(k-1), p = 0; i<k*n; i++, p++)
    for(int j = n*(k-1), q = 0; j<k*n ; j++,q++)
       A[i][j] = a_new[p][q];
}

void Mux::copy_mux(Mux ob_old, Mux ob_new)
{
   for(int i=0;i<N;i++)
     for(int j=0;j<N;j++)
        ob_new.A[i][j] = ob_old.A[i][j];
}

void Mux:: store_network(string path, unsigned int iter)
{
  Network_File(A, N, path, iter);
}

void layer_rewiring(bool** a_New, Monoplex *A_New[], Mux ob_old, Mux ob_new, unsigned long int Layer_Index)
{       
   A_New[Layer_Index]->rewire_edges(a_New,ob_old.Edges[Layer_Index]);
   ob_new.copy_mux(ob_old, ob_new);
   ob_new.insert_intra_layer(a_New, Layer_Index+1);      
}

int Rewiring(Mux ob_old, Mux ob_new, unsigned int t, unsigned int N1)
{  
   int *EDGE = new int[4];
   unsigned int l =  ob_old.l;
   Monoplex *A_New[l];
   for(int i=0;i<l;i++)
   {
     A_New[i] = new Monoplex(N1);
     A_New[i]->give_mat();
   }  
   ob_old.extract_mux_layer(A_New,l);
   
   bool **a_New = new bool*[N1];
   Allocate_2Darray(a_New, N1, N1);    

   if(t==1)
   { 
     layer_rewiring(a_New, A_New, ob_old, ob_new, 0);
   }
   if(t==2)
   {
     unsigned long int Layer_Index = Random_Index(1);
     if(Layer_Index == 0)   
        layer_rewiring(a_New, A_New, ob_old, ob_new, Layer_Index);
     else
      layer_rewiring(a_New, A_New, ob_old, ob_new, Layer_Index);
   }
   if(t==3) 
   {
     unsigned long int Layer_Index = Random_Index(2);
     if (Layer_Index == 0)
       layer_rewiring(a_New, A_New, ob_old, ob_new, Layer_Index);
     else if(Layer_Index == 1)
       layer_rewiring(a_New, A_New, ob_old, ob_new, Layer_Index); 
     else
       layer_rewiring(a_New, A_New, ob_old, ob_new, Layer_Index);
   }
   if(t==4) 
   {
     unsigned long int Layer_Index = Random_Index(3);
     if (Layer_Index == 0)
       layer_rewiring(a_New, A_New, ob_old, ob_new, Layer_Index);
     else if(Layer_Index == 1)
       layer_rewiring(a_New, A_New, ob_old, ob_new, Layer_Index); 
     else if(Layer_Index == 2)
       layer_rewiring(a_New, A_New, ob_old, ob_new, Layer_Index);     
     else
       layer_rewiring(a_New, A_New, ob_old, ob_new, Layer_Index);
   }
   if(t==5) 
   {
     unsigned long int Layer_Index = Random_Index(4);
     if (Layer_Index == 0)
       layer_rewiring(a_New, A_New, ob_old, ob_new, Layer_Index);
     else if(Layer_Index == 1)
       layer_rewiring(a_New, A_New, ob_old, ob_new, Layer_Index); 
     else if(Layer_Index == 2)
       layer_rewiring(a_New, A_New, ob_old, ob_new, Layer_Index);     
     else if(Layer_Index == 3)
       layer_rewiring(a_New, A_New, ob_old, ob_new, Layer_Index);
     else 
       layer_rewiring(a_New, A_New, ob_old, ob_new, Layer_Index);
   }
   if(t==6) 
   {
     unsigned long int Layer_Index = Random_Index(5);
     if(Layer_Index == 0)
       layer_rewiring(a_New, A_New, ob_old, ob_new, Layer_Index);
     else if(Layer_Index == 1)
       layer_rewiring(a_New, A_New, ob_old, ob_new, Layer_Index); 
     else if(Layer_Index == 2)
       layer_rewiring(a_New, A_New, ob_old, ob_new, Layer_Index);     
     else if(Layer_Index == 3)
       layer_rewiring(a_New, A_New, ob_old, ob_new, Layer_Index);
     else if(Layer_Index == 4)
       layer_rewiring(a_New, A_New, ob_old, ob_new, Layer_Index);
     else
       layer_rewiring(a_New, A_New, ob_old, ob_new, Layer_Index); 
   } 
   for(int i=0;i<l;i++)
      A_New[i]->free_mat(); 
   Free_2Darray(a_New, N1);
}

