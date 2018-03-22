#include "figures.h"


using namespace figures;

void Point::setX(double x) {
	this->x_ = x;
}

void Point::setY(double y) {
	this->y_ = y;
}

double Point::getX() {
	return this->x_;
}

double Point::getY() {
	return this->y_;
}

Point::Point(double x, double y) : x_{x}, y_{y} {}

void Line::setPointA(Point p) {
	this->pa_ = p;
}

void Line::setPointB(Point p) {
	this->pb_ = p;
}

Point Line::getPointA() {
	return this->pa_;
}

Point Line::getPointB() {
	return this->pb_;
}

std::vector<double> Line::getParametrs() {
	std::vector<double> v;
	v.push_back(this->a_);
	v.push_back(this->b_);
	v.push_back(this->c_);
	return v;
}

void Line::setParametrs() {
	this->a_ = this->getPointA().getY() - this->getPointB().getY();
	this->b_ = this->getPointB().getX() - this->getPointA().getX();
	this->c_ = this->getPointA().getX() * this->getPointB().getY()
			   - this->getPointB().getX() * this->getPointA().getY();
}

std::vector<Point> Line::intersect(Figure &_cf) {
	return _cf.intersect(*this);
};

std::vector<Point> Line::intersect(Line &_cl) {
	std::cout << "Line intersects Line" << std::endl;

	std::vector<Point> answer;

	std::vector<double> v = this->getParametrs();
	double a1 = v[0], b1 = v[1], c1 = v[2];
	v.clear();

	v = _cl.getParametrs();
	double a2 = v[0], b2 = v[1], c2 = v[2];
	v.clear();

	if ((a1 == a2) && (b1 == b2)) {
		if (c1 == c2) {
			// Lines belong to one beeline, so there is no intersection points

			return answer;
		}
		// Lines are parallel, so there is no intersection points
		return answer;
	}

	double d = a2 * b1 - a1 * b2, x = (b2 * c1 - b1 * c2) / d, y = (a1 * c2 - a2 * c1) / d;

	if ((a1 * x + b1 * y + c1) == (a2 * x + b2 * y + c2)) {
		answer.emplace_back(Point(x, y));
	}

	return answer;
}

std::vector<Point> Line::intersect(Circle &_cc) {
	std::cout << "Line intersects Circle" << std::endl;
	std::vector<Point> v;
	return v;
}

std::vector<Point> Line::intersect(PolyLine &_cp) {
	std::cout << "Line intersects PolyLine" << std::endl;
	std::vector<Point> v;
	return v;
}

Line::Line(Point a, Point b) : pa_{a}, pb_{b} {
	this->setParametrs();
}