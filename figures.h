#ifndef INTERSECTIONS_FIGURES_H
#define INTERSECTIONS_FIGURES_H

#include "figures.h"
#include <vector>
#include <iostream>

using std::vector;

namespace figures {

class Point;

class Figure;

class Segment;

class Circle;

class PolyLine;

class Point {
 public:
	void setX(const double& x);

	void setY(const double& y);

	double getX() const;

	double getY() const;

	Point(const double& x, const double& y);

 private:
	double x_;
	double y_;
};

class Figure {
 public:
	virtual double length() const = 0;

	virtual vector<Point> intersect(const Segment&) const = 0;

	virtual vector<Point> intersect(const Circle&) const = 0;

	virtual vector<Point> intersect(const PolyLine&) const = 0;

	virtual vector<Point> intersect(const Figure&) const = 0;

	virtual ~Figure() = default;
};

class Segment : Figure {
 public:
	void setPointA(const Point& p);

	void setPointB(const Point& p);

	Point getPointA() const;

	Point getPointB() const;

	vector<double> getParametrs() const;

	double length() const override;

	vector<Point> intersect(const Segment&) const override;

	vector<Point> intersect(const Circle&) const override;

	vector<Point> intersect(const PolyLine&) const override;

	vector<Point> intersect(const Figure&) const override;

	Segment(const Point& a, const Point& b);

 private:
	void recalculateParametrs();

	void recalculateLength();

	Point pa_;
	Point pb_;
	double a_, b_, c_, len_;
};

class Circle : Figure {
 public:
	void setCenter(const Point& c);

	Point getCenter() const;

	void setRadius(const double& r);

	double getRadius() const;

	double length() const override;

	vector<Point> intersect(const Segment&) const override;

	vector<Point> intersect(const Circle&) const override;

	vector<Point> intersect(const PolyLine&) const override;

	vector<Point> intersect(const Figure&) const override;

	Circle(const Point& c, const double& r);

 private:
	void recalculateLength();

	Point c_;
	double r_;
	double len_;
};

class PolyLine : Figure {
 public:
	double length() const override;

	vector<Point> *getPoints() const;

	vector<Point> intersect(const Segment&) const override;

	vector<Point> intersect(const Circle&) const override;

	vector<Point> intersect(const PolyLine&) const override;

	vector<Point> intersect(const Figure&) const override;

	explicit PolyLine(vector<Point> *points);
 private:
	void setLength();
	vector<Point> *points_;
	double len_;
};
} // end namespace figures
#endif //INTERSECTIONS_FIGURES_H