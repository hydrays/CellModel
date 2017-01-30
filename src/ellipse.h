#ifndef ELLIPSE
#define ELLIPSE

#include <math.h>

class Ellipse
{
public:
    int ellipse_id;
    double c1;
    double c2;
    double v1, v2;
    double u1, u2;
    double a, b;
    double area;

public:
    //Ellipse(){};
    //~Ellipse(){};
    
    Ellipse(int id, double c1, double c2, double v1, double v2, double u1, double u2)
	: ellipse_id(id), c1(c1), c2(c2), v1(v1), v2(v2), u1(u1), u2(u2)
	{
	    if (v1*u1 + v2*u2 != 0)
	    {
		std::cout << "Elliplse error: eigenvectors not orthogonal \n";
		getchar();
	    }
	    a = sqrt(v1*v1 + v2*v2);
	    b = sqrt(u1*u1 + u2*u2);
	}
    
    bool isInside(const double x, const double y)
	{
	    double newx = ((x-c1)*v1 + (y-c2)*v2)/a;
	    double partx = newx*newx/(a*a);
	    if ( partx > 1 )
	    {
		return false;
	    }
	    double newy = ((x-c1)*u1 + (y-c2)*u2)/b;
	    double party = newy*newy/(b*b);
	    if ( party > 1 )
	    {
		return false;
	    }
	    return (partx + party) <= 1;
	}
};

#endif //ELLIPSE
