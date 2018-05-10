#ifndef INTERSECTIONS_FIGURES_H
#define INTERSECTIONS_FIGURES_H

#include <cstdlib>
#include <vector>

using std::vector;

namespace figures {

struct Point {
	Point();
	Point(double x, double y);
	double x, y;
	const bool operator<(const Point& rhs) const;
	const bool operator==(const Point& rhs) const;
};

class Figure;
class Segment;
class Circle;
class Polyline;

class Figure {
 public:
	virtual double length() const = 0;
	virtual bool belong(const Point&) const = 0;
	virtual vector<Point> intersect(const Figure&) const = 0;
	virtual vector<Point> intersect(const Segment&) const = 0;
	virtual vector<Point> intersect(const Circle&) const = 0;
	virtual vector<Point> intersect(const Polyline&) const = 0;
	virtual ~Figure() = default;
};

class Segment : Figure {
 public:
	Segment() = default;
	Segment(const Point& a, const Point& b);
	Segment(const double& x1, const double& x2, const double& y1, const double& y2);
	vector<double> getParameters() const;
	Point getA() const;
	Point getB() const;
	void setA(Point& a);
	void setB(Point& b);
	double length() const override;
	bool belong(const Point&) const override;
  vector<Point> intersect(const Figure&) const override;
	vector<Point> intersect(const Segment&)  const override;
	vector<Point> intersect(const Circle&) const override;
	vector<Point> intersect(const Polyline&) const override;
 private:
	void recalculateParameters();
	void recalculateLength();
	Point pa_, pb_;
	double a_ = 0, b_ = 0, c_ = 0, length_ = 0;
};

class Circle: Figure {
 public:
	Circle() = default;
	Circle(const Point& c, const double& r);
	double getRadius() const;
	Point getCenter() const;
	void setRadius(const double& r);
	void setCenter(const Point& c);
	double length() const override;
	bool belong(const Point&) const override;
	vector<Point> intersect(const Figure&) const override;
	vector<Point> intersect(const Segment&) const override;
	vector<Point> intersect(const Circle&) const override;
	vector<Point> intersect(const Polyline&) const override;
 private:
	void recalculateLength();
	Point c_;
	double r_ = 0, length_ = 0;
};

class Polyline: Figure {
 public:
  Polyline() = default;
  explicit Polyline(const vector<Point>& points);
  vector<Point> getPoints() const;
  double length() const override;
  bool belong(const Point&) const override;
  vector<Point> intersect(const Figure&) const override;
  vector<Point> intersect(const Segment&) const override;
  vector<Point> intersect(const Circle&) const override;
  vector<Point> intersect(const Polyline&) const override;
 private:
  void recalculateLength();
  vector<Point> points_;
  double length_ = 0;
};
} // end namespace figures

#endif //INTERSECTIONS_FIGURES_H
