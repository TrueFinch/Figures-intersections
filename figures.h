#ifndef INTERSECTIONS_FIGURES_H
#define INTERSECTIONS_FIGURES_H

#include "figures.h"
#include <vector>
#include <iostream>

using std::vector;

namespace figures {

class Point;

class Figure;

class Line;

class Circle;

class PolyLine;

class Point {
 public:
  void setX(double x);

  void setY(double y);

  double getX();

  double getY();

  Point(double x, double y);

 private:
  double x_;
  double y_;
};

class Figure {
 public:
  virtual int length() = 0;

  virtual vector<Point> intersect(Line &) = 0;

  virtual vector<Point> intersect(Circle &) = 0;

  virtual vector<Point> intersect(PolyLine &) = 0;

  virtual vector<Point> intersect(Figure &) = 0;
};

class Line : Figure {
 public:
  void setPointA(Point p);

  void setPointB(Point p);

  Point getPointA();

  Point getPointB();

  vector<double> getParametrs();

  void setParametrs();

  virtual int length();

  virtual vector<Point> intersect(Line &);

  virtual vector<Point> intersect(Circle &);

  virtual vector<Point> intersect(PolyLine &);

  virtual vector<Point> intersect(Figure &);

  Line(Point a, Point b);

 private:
  Point pa_;
  Point pb_;
  double a_, b_, c_;
};

class Circle : Figure {
 public:
  void setCenter(Point c);

  void setRadius(Point c);

  Point getCenter();

  int getRadius();

  virtual int length();

  virtual vector<Point> intersect(Line &);

  virtual vector<Point> intersect(Circle &);

  virtual vector<Point> intersect(PolyLine &);

  virtual vector<Point> intersect(Figure &);

 private:
  Point c_;
  int r_;

};

class Polyline : Figure {

};
} // end namespace figures
#endif //INTERSECTIONS_FIGURES_H