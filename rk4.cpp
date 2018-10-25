#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
using namespace std;


struct RK
{
  float b= 8.0/3.0, r= 28, s=10;
  float x,y,z,dh,h=0.02,to;
  float k1,k2,k3,k4,m1,m2,m3,m4,p1,p2,p3,p4;
  int i,j,contador=0;

  void Paso ()
  {
    cout << "\nIngresa el tamano de paso h: "<< endl; //Pide la dimension tanto de la matriz, como de los vectores
    cin>>dh;
  }

  void line()//Iteraciones
  {
    for(int i=1;i<41;i++)
    cout<<"--";
    cout<<"\n";
  }

  float F1 (float x, float y) //dx/dt = s(y-x) = k
  {
  return (s*(y - x));
  }

  float F2 (float x, float y, float z) //dy/dt = x*(r-z)-y = m
  {
    return ((x*(r-z)) - y);
  }

  float F3 (float x, float y, float z) // dz/dt = x*y - b*z = p
  {
    return ((x*y) - (b*z));
  }

  void RK4 ()
  {
    //condiciones iniciales
    x=1,y=1,z=1,to=0;
    int contador=0;

    cout <<setw(65)<<"Tabulacion de valores de t, x, y, z usando RK4:\n";
    line();
    cout<<setw(15)<<"t"<<setw(15)<<"X"<<setw(15)<<"Y"<<setw(15)<<"Z\n";
    line();


    //for(h=dh;h<=5.0;h=h+dh)
    while(to<=30)
    {
      cout<<setw(15)<<to<<setw(15)<<x<<setw(15)<<y<<setw(15)<<z<<endl;

      k1= F1(x,y);
      m1= F2(x,y,z);
      p1= F3(x,y,z);

      k2= F1(x+0.5*k1*h  , y+0.5*m1*h  );
      m2= F2(x+0.5*k1*h  , y+0.5*m1*h , z+0.5*p1*h );
      p2= F3(x+0.5*k1*h  , y+0.5*m1*h  , z+0.5*p1*h );

      k3= F1(x+0.5*k2*h  , y+0.5*m2*h );
      m3= F2(x+0.5*k2*h  , y+0.5*m2*h , z+0.5*p2*h  );
      p3= F3(x+0.5*k2*h  , y+0.5*m2*h , z+0.5*p2*h );

      k4= F1(x+k3*h , y+m3*h);
      m4= F2(x+k3*h , y+m3*h , z+p3*h);
      p4= F3(x+k3*h  , y+m3*h , z+p3*h) ;

      x+= h*(k1+2*k2+2*k3+k4)/6;
      y+= h*(m1+2*m2+2*m3+m4)/6;
      z+= h*(p1+2*p2+2*p3+p4)/6;

    /*  cout << k1 << endl;
      cout << k2 << endl;
      cout << k3 << endl;
      cout << k4 << endl;*/

      //cout<<setw(15)<<h<<setw(15)<<x<<setw(15)<<y<<setw(15)<<z<<endl;
      to += h;
      contador ++;

    }
    cout << "Numero de iteraciones: " << contador*15 <<endl;
  }

  void Euler ()
  {
    x=1,y=1,z=1,to=0;
    int contador=0;

    cout <<setw(65)<<"Tabulacion de valores de t, x, y, z usando Euler:\n";
    line();
    cout<<setw(15)<<"t"<<setw(15)<<"X"<<setw(15)<<"Y"<<setw(15)<<"Z\n";
    line();

    while(to<=30)
    {
    cout<<setw(15)<<to<<setw(15)<<x<<setw(15)<<y<<setw(15)<<z<<endl;

    k1= F1(x,y);
    m1= F2(x,y,z);
    p1= F3(x,y,z);

    x= x + h*k1;
    y= y + h*m1;
    z= z + h*p1;

    to += h;
    contador ++;
    }
    cout << "Numero de iteraciones: " << contador*6 <<endl;

  }

  void RK2 ()
  {
    x=1,y=1,z=1,to=0;
    int contador=0;

    cout <<setw(65)<<"Tabulacion de valores de t, x, y, z usando RK2:\n";
    line();
    cout<<setw(15)<<"t"<<setw(15)<<"X"<<setw(15)<<"Y"<<setw(15)<<"Z\n";
    line();


    while(to<=30)
    {
      cout<<setw(15)<<to<<setw(15)<<x<<setw(15)<<y<<setw(15)<<z<<endl;

      k1= F1(x,y);
      m1= F2(x,y,z);
      p1= F3(x,y,z);

      k2= F1(x+0.5*k1*h  , y+0.5*m1*h );
      m2= F2(x+0.5*k1*h  , y+0.5*m1*h , z+0.5*p1*h );
      p2= F3(x+0.5*k1*h  , y+0.5*m1*h  , z+0.5*p1*h );

      x= x + h*(k1+k2);
      y= y + h*(m1+m2);
      z= z + h*(p1+m2);

      to += h;
      contador ++;
      }
        cout << "Numero de iteraciones: " << contador*9 <<endl;
  }

};



int main()
{
  RK rk;

  rk.Euler();
  rk.RK2();
  rk.RK4();

  return 0;
}


// https://people.sc.fsu.edu/~jburkardt/py_src/lorenz_ode/lorenz_ode.html
//https://stackoverflow.com/questions/41101960/runge-kutta-4-in-c-for-lorenz-system
// https://github.com/oliviercotte/Linear-differential-equation/blob/master/rk4.m
// https://github.com/MarioAriasGa/lorenz/blob/master/solver.cpp
// https://www.mathworks.com/matlabcentral/answers/360751-solving-lorenz-attractor-equations-using-runge-kutta-rk4-method
// https://stackoverflow.com/questions/9338804/runge-kutta-in-c-for-lorenz-equation
