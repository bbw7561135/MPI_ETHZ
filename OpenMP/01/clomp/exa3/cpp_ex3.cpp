#include "cpp_ex3.h"

#define PI 3.141592654

#define STEPS 1000.0

int main()
{
   House house;
   int arrived;
   std::vector<House*> * movie = new std::vector<House*>;
   house.setup();
   house.print();


   #pragma omp parallel
   {
       int i;
       Point start(-30.0, -14.0, -30.0 ); 
       Point end(30.0, 14.0, 30.0); 
       House *h;
       int n_threads = omp_get_num_threads();

       double x_step, y_step, z_step;

       x_step = (end.x - start.x) / STEPS;
       y_step = (end.y - start.y) / STEPS;
       z_step = (end.z - start.z) / STEPS;

       Point steps(x_step, y_step, z_step);

       #pragma omp single
       {
	   arrived = 0;
	   h = new House(house);
	   movie->push_back(h);
       }

       /* move the house along a line and back again */
       for (i=0; i<STEPS; i++)
       {
	   house.move(start); 

	   #pragma omp critical
	   {
	       arrived++;
	       if (arrived == n_threads)
	       {
		   h = new House(house);
		   movie->push_back(h);
		   arrived = 0;
	       }
	   }

	   #pragma omp barrier

	   start.move(steps);
       }

       house.move(start); 
       #pragma omp critical
       {
	   arrived++;
	   if (arrived == n_threads)
	   {
	       h = new House(house);
	       movie->push_back(h);
	   }
       }
   }
   printf("\nThe following should be the same\n");
   house.print();
}

/**************************************************************
 * Setup
 **************************************************************/
void House::setup()
{
   Point a(0.0,0.0,0.0);
   Point b(0.0,1.0,0.0);
   Point c(1.0,1.0,0.0);
   Point d(1.0,0.0,0.0);
   Point e(0.0,0.0,1.0);
   Point f(0.0,1.0,1.0);
   Point g(1.0,1.0,1.0);
   Point h(1.0,0.0,1.0);
   Point i(0.5,1.5,0.0);
   Point j(0.5,1.5,1.0);

   frontB = new Rectangle(a, b, c, d);
   backB = new Rectangle(e,f,g,h);
   sideR = new Rectangle(c,d,h,g);
   sideL = new Rectangle(a,b,f,e);
   roofR = new Rectangle(c,i,j,g);
   roofL = new Rectangle(b,i,j,e);
   bottom = new Rectangle(a,e,h,d);
   frontT = new Triangle(b,c,i);
   backT = new Triangle(f,g,j);
}

Point::Point(double a, double b, double c)
{
   x = a;
   y = b;
   z = c;
}

Point::Point(Point &point)
{
   x = point.x;
   y = point.y;
   z = point.z;
}

Triangle::Triangle()
{
   p1 = p2 = p3 = NULL;
}

Triangle::Triangle(Triangle &t)
{
   p1 = new Point(*t.p1);
   p2 = new Point(*t.p2);
   p3 = new Point(*t.p3);
}

Triangle::Triangle(Point &p1,Point &p2,Point &p3)
{
   this->p1 = new Point(p1);
   this->p2 = new Point(p2);
   this->p3 = new Point(p3);
}

Rectangle::Rectangle()
{
   p1 = p2 = p3 = p4 = NULL;
}

Rectangle::Rectangle(Rectangle &r)
{
   p1 = new Point(*r.p1);
   p2 = new Point(*r.p2);
   p3 = new Point(*r.p3);
   p4 = new Point(*r.p4);
}

Rectangle::Rectangle(Point &p1, Point &p2, Point &p3, Point &p4)
{
   this->p1 = new Point(p1);
   this->p2 = new Point(p2);
   this->p3 = new Point(p3);
   this->p4 = new Point(p4);
}

House::House()
{
    frontB = backB = sideR = sideL = roofR = roofL = NULL;
    frontT = backT = NULL;
}

House::House(House &h)
{
   frontB = new Rectangle(*h.frontB);
   backB = new Rectangle(*h.backB);
   sideR = new Rectangle(*h.sideR);
   sideL = new Rectangle(*h.sideL);
   roofR = new Rectangle(*h.roofR);
   roofL = new Rectangle(*h.roofL);
   frontT = new Triangle(*h.frontT);
   backT = new Triangle(*h.backT);
}

/**************************************************************
 * Routines for moving objects.  Objects are moved in space by 
 * the amount specified in Point.
 **************************************************************/
void Point::move(Point &p)
{
   x += p.x;
   y += p.y;
   z += p.z;
}

void Triangle::move(Point &p)
{
   p1->move(p); 
   p2->move(p); 
   p3->move(p); 
}

void Rectangle::move(Point &p)
{
   p1->move(p); 
   p2->move(p); 
   p3->move(p); 
   p4->move(p); 
}

void House::move(Point &p)
{
   #pragma omp single nowait
   {
      frontB->move(p);
   }
   #pragma omp single nowait
   {
      backB->move(p);
   }
   #pragma omp single nowait
   {
      sideR->move(p);
   }
   #pragma omp single nowait
   {
      sideL->move(p);
   }
   #pragma omp single nowait
   {
      roofR->move(p);
   }
   #pragma omp single nowait
   {
      roofL->move(p);
   }
   #pragma omp single nowait
   {
      frontT->move(p);
   }
   #pragma omp single nowait
   {
      backT->move(p);
   }
   #pragma omp single nowait
   {
      bottom->move(p);
   }

}


/**************************************************************
 * Destructors 
 **************************************************************/
Triangle::~Triangle()
{
   if (p1)
   {
      delete p1;
      delete p2;
      delete p3;
   }
}
Rectangle::~Rectangle()
{
   if (p1)
   {
      delete p1;
      delete p2;
      delete p3;
      delete p4;
   }
}

House::~House()
{
   if (frontB)
   {
	delete frontB;
	delete backB;
	delete sideR;
	delete sideL;
	delete roofR;
	delete roofL;
	delete frontT;
	delete backT;
   }
}

/**************************************************************
 * Routines for printing points 
 **************************************************************/
#define ROUND(x) ((((10.0*x) - (unsigned long long)(10.0*x)) > 0.5) ? (double)((unsigned long long)(10.0*x) + 1)/10.0 : (double)((unsigned long long)(10.0*x))/10.0)
void Point::print()
{

   printf("(%lf,%lf,%lf)", ROUND(x), ROUND(y),ROUND(z));
}

void Triangle::print()
{
   printf ("[");
   p1->print();
   printf (",");
   p2->print();
   printf (",");
   p3->print();
   printf ("]");
}

void Rectangle::print()
{
   printf ("[");
   p1->print();
   printf (",");
   p2->print();
   printf (",");
   p3->print();
   printf (",");
   p4->print();
   printf ("]");
}

void House::print()
{
   printf("Front: ");
   frontB->print();
   printf (",");
   frontT->print();
   printf("\nBack:  ");
   backB->print();
   printf (",");
   backT->print();
   printf("\nLeft:  ");
   sideL->print();
   printf (",");
   roofL->print();
   printf("\nRight: ");
   sideR->print();
   printf (",");
   roofR->print();
   printf ("\n\n");

}

