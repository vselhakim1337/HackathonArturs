#include <iostream>
#include <math.h>
#include <limits>
#include <float.h>

#define DEBUG 1
#define EPS 0.0000001

struct Point
{
    double x;
    double y;
    Point(double x, double y) : x(x), y(y) {}
};
std::ostream& operator << ( std::ostream& outs, const Point &p)
{
    return outs << "(" << p.x << "," << p.y << ")";
}

inline double dsin(const double angle){return sin(angle * M_PI/180.0);}
inline double dcos(const double angle){return cos(angle * M_PI/180.0);}
inline double datan(const double fract){ return atan(fract) * 180.0/M_PI;}


class Solver 
{
protected:
    double angle = 0;
public:
    virtual Point next(const double signX, const double signY, 
        const Point &prev) = 0;
    virtual void reduceStep() = 0;
    virtual void updateAngle(const double signY) = 0;
    double getAngle()
    {
        return angle;
    }
};
class RotationMatrix : public Solver
{
private:
    double step = 45;
public:
    Point next(const double signX, const double signY,
        const Point &prev)
    {
        const double angle = signX*signY*step;
        double nextX = dcos(angle) * prev.x + (-dsin(angle) * prev.y);
        double nextY = dsin(angle) * prev.x + dcos(angle) * prev.y;
        return Point(nextX, nextY);
    }
    void reduceStep()
    {
        step = step / 2.0;
    }
    void updateAngle(const double signY)
    {
        angle += signY * step;
    }
};
class Cordic : public Solver
{
private:
    int iter = 0;
    double powOf2fact = 1;
    //const double factorK = 0.6072529;
    const double factorK = 1.0;
    const double fact00 = 1.0;
    double fact01 = 0.0;
    double fact10 = 0.0;
    const double fact11 = 1.0;

public:
    Point next(const double signX, const double signY,
        const Point &prev)
    {
        fact01 = -(signX * signY * powOf2fact);
        fact10 = signX * signY * powOf2fact;
        double nextX = fact00 * prev.x + fact01 * prev.y;
        double nextY = fact10 * prev.x + fact11 * prev.y;

#ifdef DEBUG
        printf("Iteration: %d\n", iter);
        printf("Angle: %d\n", iter);
        printf("(%f,%f)\n", nextX, nextY);
#endif
        nextX*=factorK;
        nextY*=factorK;
        // cosine(theta) multiplier
        //nextX*=dcos(datan(powOf2fact));
        //nextY*=dcos(datan(powOf2fact));

#ifdef DEBUG        
        printf("(%f,%f)\n", nextX, nextY);
#endif
        return Point(nextX, nextY);
    }
    void reduceStep()
    {
        ++iter;
        powOf2fact = 1.0/pow(2, iter);  
    }
    void updateAngle(const double signY)
    {
#ifdef DEBUG
        printf("Angle: %f\n\n", angle);
#endif
        angle = angle - signY * datan(powOf2fact);
    }

};


double approx(const Point &bPoint, Solver *solver)
{
    Point nPoint(bPoint.x, bPoint.y);

#ifdef DEBUG
    int iter = 50;
#endif

    double signX = bPoint.x > 0 ? 1 : -1.0;
    while (std::abs(nPoint.y) > EPS)
    {
        double signY = nPoint.y > 0 ? -1.0 : 1.0;
        nPoint = solver->next(signX, signY, nPoint);
        solver->updateAngle(signY);
        solver->reduceStep();
#ifdef DEBUG
        if (--iter <= 0)
            break;
#endif
    }
    return solver->getAngle();
}

int main()
{
    double x = 0.5;
    double y = 0.866025;

    double x0 = 0.866;
    double y0 = 0.5;
    
    Point p(x,y);
    Point p0(x0,y0);

    RotationMatrix rmSolver;
    Cordic cordicSolver;

    //std::cout << approx(p, &rmSolver) << '\n';
    std::cout << approx(p0, &cordicSolver) << '\n';

    return 0;
}