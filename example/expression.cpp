#include <iostream>
#include <numeric>
#include <map>
#include <fstream>
#include <chrono>
#include "geometry/point.h"
#include "geometry/shapes.h"

int main()
{
    cp::geometry::point<double> A(1,2),B(0,4);
    cp::geometry::line_segment<double> AB(A,B);
    auto L=bisector(AB);

}