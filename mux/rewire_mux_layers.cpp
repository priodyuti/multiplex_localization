#include<iostream>
#include<fstream>
#include<cstdlib>
#include"rewire_mux_layers.h"
#include<math.h>
#include<ctime>
#include<stdlib.h>
#include<iomanip>
#include<sstream>
#include<string>
#include<chrono>

//Paper: Sarika Jalan and Priodyuti Pradhan, Localization of multilayer networks by optimized single-layer rewiring, Phys. Rev. E 97, 042314 (2018)

// Written: Priodyuti Pradhan, Complex System Lab,IIT Indore
// contact: priodyutipradhan@gmail.com  

using namespace std;
using namespace std::chrono;

int main()
{
  high_resolution_clock::time_point t_start = high_resolution_clock::now();

  unsigned int k, N1, l, t;
  unsigned long int Iter;
  double p;

  cout<<" Number of layers: ";  cin >> l;
  cout<<" Number of layers to be rewired: "; cin >> t;  
  cout<<" Number of nodes in a layer: ";  cin >> N1;
  cout<<" Average degree: ";  cin>>k;
  cout<<" Number of Iterations: "; cin>>Iter;
  p = k/double(N1);
  cout<<" Connection probability: "<<p<<endl;
  
  string path = "//home//priodyuti//Desktop//Mux-4-1//";
  for(unsigned int i=1;i<=l;i++)
    path = path + "ER-";
  path = path + to_string(t)+"-"+to_string(k)+"-"+to_string(N1)+"//";
 
  create_directory(path);

  string ipr_path = path + "ipr.txt";
  ofstream fout(ipr_path);
  if(!fout) { cout<<"File not Opened"<<endl; return 1; }

  string path1 = path + "Networks//";
  create_directory(path1);

  unsigned int N = l*N1;
  
  Mux ob(N,l,N1,k);
  Monoplex* A[l];
  for(int i=0;i<l;i++)
  {
     A[i] = new Monoplex(N1,k);
     A[i]->give_mat();
     A[i]->mat_init_random();
     //A[i]->mat_display();
     A[i]->find_ipr();
  }
  ob.init_mux();
  
  ob.construct_mux(A);  
  //ob.display_mux();
  double c1 = ob.find_ipr_mux();
  //ob.store_network(path1, 1);
   
  double Tem = 0.9;   //Temparature parameture for simulated annealing
  
  for(unsigned long int iter=2; iter<=Iter; iter++)
  {
    Mux ob_New(N,l,N1,k);
    ob_New.init_mux();

    int Flag = 0;
    Rewiring(ob, ob_New, t, N1); 
    //ob_New.display_mux();
    double c2 = ob_New.find_ipr_mux();
    
    if(c2>c1)
    {
      ob_New.copy_mux(ob_New, ob);
      //Matrix_Copy(g_New, g, N);
      c1 = c2; 
      Flag = 1;
    }
    else
    { 
      double Ex = exp((c2-c1)/(100*Tem));
      double Ran = ran_number();
      if(Ran<Ex)
      {
        ob_New.copy_mux(ob_New, ob);
        //Matrix_Copy(g_New, g, N);
        c1 = c2;
      } 
    }
    fout<<Flag<<' '<<c1<<' '<<c2<<endl;
    Tem = Tem*0.998;
    //if(iter % 5000 == 1)
    //   ob.store_network(path1, iter);
    ob_New.free_mux();
  }

   ob.store_network(path, 1);
   ob.free_mux();
   for(int i=0;i<l;i++)
        A[i]->free_mat(); 

 cout<<"\nCompletes... "<<endl;
 //:::::::::::::::::::Measuring time taken by the abobe call :::::::::::::::
 high_resolution_clock::time_point t_end = high_resolution_clock::now();
 auto duration1 = duration_cast<seconds>( t_end - t_start ).count();
 cout <<"Full Program Time Duration: " <<duration1<<" Seconds"<<endl;
 //:::::::::::::::::::::: END  :::::::::::::::::::::::::::::::::::::::::::::
 return 0;
}


