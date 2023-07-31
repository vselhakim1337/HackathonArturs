#include <iostream>
#include <math.h>
#include <limits>
#include <float.h>

#define DEBUGn 1
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
    virtual void reset() = 0;
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
    void reset(){}
};
class Cordic : public Solver
{
private:
    int iter = 0;
    double powOf2fact = 1;
    const double factorK = 0.6072529;
    //const double factorK = 1.0;
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
        printf("Angle: %f\n\n", (-signY * datan(powOf2fact)));
#endif
        angle = angle - signY * datan(powOf2fact);
    }
    void reset()
    {
        angle = 0.0;
        iter = 0;
        powOf2fact = 1;
        fact01 = 0.0;
        fact10 = 0.0;
    }
};


double approx(const Point &bPoint, Solver *solver)
{
    solver->reset();
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

#define MAX_ITER 2

#include <cstdint>
std::int32_t tan_table[] = {
    572840,
 267745,
 131759,
  65621,
  32779,
  16385,
   8192,
   4096,
   2048,
   1024
};

#define FLOAT_TO_FP(x, p) ((std::int32_t)(x * (double)(1 << (p))))
#define FP_TO_FLOAT(x, p) ((double)(x) / (double)(1 << (p)))

double float_cordic_approx_angle(double x, double y) {
    double xp = x;
    double yp = y;
    double zp = 0.0;

    double xn = 0.0;
    double yn = 0.0;
    double zn = 0.0;
    for (std::size_t i = 0; i < MAX_ITER; ++i) {
        double k = 1.0/pow(2.0, i+1);
        double sign = yp < 0.0 ? 1.0 : -1.0;
        xn = xp - sign*k*yp;
        yn = yp + sign*k*xp;
        zn = zp - sign*tan(k);

        xp = xn;
        yp = yn;
        zp = zn;
    }

    return zn;
}


double fixed_point_cordic_approx_angle(double x, double y) {
    std::int32_t xp = FLOAT_TO_FP(x, 20);
    std::int32_t yp = FLOAT_TO_FP(y, 20);
    std::int32_t zp = 0.0;

    std::int32_t xn = 0.0;
    std::int32_t yn = 0.0;
    std::int32_t zn = 0.0;
    for (std::size_t i = 0; i < MAX_ITER; ++i) {
        std::int32_t sign = yp < 0 ? 1 : -1;
        xn = xp - sign*(yp >> 1);
        yn = yp + sign*(xp >> 1);
        zn = zp - sign*tan_table[i];

        xp = xn;
        yp = yn;
        zp = zn;
    }

    return FP_TO_FLOAT(zn, 20);
}

int main()
{
    double x = 0.5;
    double y = 0.866025;

    double x0 = 0.866;
    double y0 = 0.5;

    double x1 = 0.951057;
    double y1 = 0.309017;


    Point p0(x0,y0);
    Point p1(x1,y1);

    RotationMatrix rmSolver;
    Cordic cordicSolver;

    //std::cout << approx(p, &rmSolver) << '\n';

    // std::cout << approx(p0, &cordicSolver) << '\n';
    // std::cout << "====================" << std::endl;
    // std::cout << approx(p0, &cordicSolver) << '\n';

    std::cout << "=====================" << std::endl;
    std::cout << float_cordic_approx_angle(x0,y0) / 3.14159265358979323846 * 180 << '\n';
    std::cout << "=====================" << std::endl;
    std::cout << fixed_point_cordic_approx_angle(x0,y0) / 3.14159265358979323846 * 180 << '\n';


    return 0;
}