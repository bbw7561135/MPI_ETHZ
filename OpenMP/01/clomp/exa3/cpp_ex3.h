#include <stdio.h>
#include <vector>
#include <math.h>
#include <omp.h>
using namespace std;

class Point 
{
   public:
   double x,y, z;

   Point(double a, double b, double c);
   Point(Point &point);

   void move(Point &p);
   void print();
};

class Triangle
{
public:
   Point *p1, *p2, *p3;
   Triangle();
   Triangle(Point &p1,Point &p2,Point &p3);
   Triangle(Triangle &t);
   ~Triangle();

   void move(Point &p);
   void print();
};

class Rectangle
{
public:
   Point *p1, *p2, *p3, *p4;
   Rectangle();
   Rectangle(Point &p1, Point &p2, Point &p3, Point &p4);
   Rectangle(Rectangle &r);
   ~Rectangle();

   void move(Point &p);
   void print();
};

class House
{
public:
    Rectangle * frontB, * backB, * sideR, * sideL, *roofR, *roofL, *bottom;
    Triangle * frontT, *backT;

   House();
   House(House &h);
   ~House();
   void move(Point &p);
   void setup();
   void print();
};
