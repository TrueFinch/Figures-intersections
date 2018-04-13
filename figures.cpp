#include "figures.h"
#include <cmath>
#include <utility>
#include <set>
#include <numeric>

using std::cout;
using std::set;
using namespace figures;

struct comparePoints {
	bool operator()(const Point& a, const Point& b) {
		return ((a.getX() < b.getX()) && (a.getY() < b.getY()));
	}
};

void Point::setX(const double& x) {
	this->x_ = x;
}

double Point::getX() const {
	return this->x_;
}

void Point::setY(const double& y) {
	this->y_ = y;
}

double Point::getY() const {
	return this->y_;
}

bool Point::belong(const Point& _cp) const {
	return *this == _cp;
}

bool Point::belong(const Segment& _cs) const {
	return this->intersect(_cs)[0] == *this;
}

bool Point::belong(const Circle& _cc) const {
	return this->intersect(_cc)[0] == *this;
}

bool Point::belong(const PolyLine& _cp) const {
	return this->intersect(_cp)[0] == *this;
}

bool Point::belong(const Figure& _cf) const {
	return this->intersect(_cf)[0] == *this;
}

Point::Point(const double& x, const double& y) : x_{x}, y_{y} {}

const bool Point::operator==(const Point& rhs) const {
	return ((this->getX() == rhs.getX()) && (this->getY() == rhs.getY()));
}

vector<Point> Point::intersect(const figures::Point& _cp) const {
	vector<Point> answer;
	if (*this == _cp) {
		answer.emplace_back(*this);
	} else {
		answer.emplace_back(Point(this->getX() * (-1), this->getY() * (-1)));
	}
	return answer;
}

vector<Point> Point::intersect(const Segment& _cs) const {
	vector<Point> answer;
	vector<double> p = _cs.getParameters();
	double a = p[0], b = p[1], c = p[2];
	if ((a * this->getX() + b * this->getY() + c) == 0){
		answer.emplace_back(*this);
	} else {
		answer.emplace_back(Point(this->getX() * (-1), this->getY() * (-1)));
	}
	return answer;
}

vector<Point> Point::intersect(const Circle& _cc) const {
	vector<Point> answer;
	Point c = _cc.getCenter();
	double r = _cc.getRadius();
	if ((pow(this->getX() - c.getX(), 2) + pow(this->getY() - c.getY(), 2) == pow(r, 2))) {
		answer.emplace_back(*this);
	} else {
		answer.emplace_back(Point(this->getX() * (-1), this->getY() * (-1)));
	}
	return answer;
}

vector<Point> Point::intersect(const PolyLine& _cp) const {
	vector<Point> answer, points = _cp.getPoints();
	for (int i = 0; i < points.size() - 1; ++i) {
		if (this->intersect(Segment(points[i], points[i + 1]))[0] == *this) {
			continue;
		} else {
			answer.emplace_back(Point(this->getX() * (-1), this->getY() * (-1)));
			return answer;
		}
	}
	answer.emplace_back(*this);
	return answer;
}

vector<Point> Point::intersect(const Figure& _cf) const {
	return _cf.intersect(*this);
}

double Point::length() const {
	return 42;
}

void Segment::setA(const Point& p) {
	this->pa_ = p;
	this->recalculateParameters();
	this->recalculateLength();
}

Point Segment::getA() const {
	return this->pa_;
}

void Segment::setB(const Point& p) {
	this->pb_ = p;
	this->recalculateParameters();
	this->recalculateLength();
}

Point Segment::getB() const {
	return this->pb_;
}

void Segment::recalculateParameters() {
	this->a_ = this->getA().getY() - this->getB().getY();
	this->b_ = this->getB().getX() - this->getA().getX();
	this->c_ = this->getA().getX() * this->getB().getY() - this->getB().getX() * this->getA().getY();
	if (this->a_ < 0) {
		this->a_ *= (-1);
		this->b_ *= (-1);
		this->c_ *= (-1);
	} else if ((this->b_ < 0) && (this->a_ == 0)) {
		this->b_ *= (-1);
		this->c_ *= (-1);
	}
}

vector<double> Segment::getParameters() const {
	vector<double> v;
	v.push_back(this->a_);
	v.push_back(this->b_);
	v.push_back(this->c_);
	return v;
}

void Segment::recalculateLength() {
	this->len_ = sqrt(pow(this->getB().getX() - this->getA().getX(), 2)
												+ pow(this->getB().getY() - this->getA().getY(), 2));
}

double Segment::length() const {
	return this->len_;
}

vector<Point> Segment::intersect(const Point& _cp) const {
	return _cp.intersect(*this);
}

vector<Point> Segment::intersect(const Segment& _cs) const {
	vector<Point> answer;

	vector<double> p = this->getParameters();
	double a1 = p[0], b1 = p[1], c1 = p[2];
	p.clear();

	p = _cs.getParameters();
	double a2 = p[0], b2 = p[1], c2 = p[2];
	p.clear();

	if ((a1 == a2) || (b1 == b2)) {
		if (c1 == c2) {
			// Segment are located on same line
			double
					x1 = this->getA().getX(), y1 = this->getA().getY(),
					x2 = this->getB().getX(), y2 = this->getB().getY(),
					x3 = _cs.getA().getX(), y3 = _cs.getA().getY(),
					x4 = _cs.getB().getX(), y4 = _cs.getB().getY();
			if (((x3 - x1) * (x4 - x1) + (y3 - y1) * (y4 - y1)) <= 0) {
				answer.emplace_back(Point(x1, y1));
			} else if (((x3 - x2) * (x4 - x2) + (y3 - y2) * (y4 - y2)) <= 0) {
				answer.emplace_back(Point(x2, y2));
			}
			if (((x1 - x3) * (x2 - x3) + (y1 - y3) * (y2 - y3)) <= 0) {
				answer.emplace_back(Point(x3, y3));
			} else if (((x1 - x4) * (x2 - x4) + (y1 - y4) * (y2 - y4)) <= 0) {
				answer.emplace_back(Point(x4, y4));
			}
			return answer;
		}
		// Segments are parallel, so there is no intersection points
		return answer;
	}

	double d = a2 * b1 - a1 * b2, x = (b2 * c1 - b1 * c2) / d, y = (a1 * c2 - a2 * c1) / d;

	if ((a1 * x + b1 * y + c1) == (a2 * x + b2 * y + c2)) {
		answer.emplace_back(Point(x, y));
	}

	return answer;
}

vector<Point> Segment::intersect(const Circle& _cc) const {
	vector<Point> answer;

	vector<double> v = this->getParameters();
	double a1 = v[0], b1 = v[1], c1 = v[2];

	double
			p = _cc.getCenter().getX(),
			q = _cc.getCenter().getY(),
			r = _cc.getRadius(),
			d = pow(fabs(a1 * p + b1 * q + c1), 2) / (a1 * a1 + b1 * b1);

	if (d <= r) {
		double
				k = (-1) * a1 / b1,
				b = (-1) * c1 / b1,
				t = q - b,
				x = (p + k * t) / (1 + k * k),
				y = k * x + b;
		if (d == r) {
			answer.emplace_back(Point(x, y));
		} else if (d < r) {
			double
					u = sqrt((2 * k + t) * t - (p * p - r * r) * k * k + r * r) / (1 + k * k),
					x1 = x + u,
					y1 = k * x1,
					x2 = x - u,
					y2 = k * x2;
			answer.emplace_back(Point(x1, y1));
			answer.emplace_back(Point(x2, y2));
		}
	}
	return answer;
}

vector<Point> Segment::intersect(const PolyLine& _cp) const {
	vector<Point> answer, tmp, points = _cp.getPoints();
	for (int i = 0; i < points.size() - 1; ++i) {
		tmp = Segment(points[i], points[i + 1]).intersect(*this);
		for (auto j : tmp) {
			answer.push_back(j);
		}
		tmp.clear();
	}
	return answer;
}

vector<Point> Segment::intersect(const Figure& _cf) const {
	return _cf.intersect(*this);
};

Segment::Segment(const Point& a, const Point& b) : pa_{a}, pb_{b} {
	this->recalculateParameters();
	this->recalculateLength();
}

void Circle::setCenter(const Point& c) {
	this->c_ = c;
}

Point Circle::getCenter() const {
	return this->c_;
}

void Circle::setRadius(const double& r) {
	this->r_ = r;
	this->recalculateLength();
}

double Circle::getRadius() const {
	return this->r_;
}

void Circle::recalculateLength() {
	this->len_ = 2 * M_PI * this->getRadius();
}

double Circle::length() const {
	return this->len_;
}

vector<Point> Circle::intersect(const Point& _cp) const {
	return _cp.intersect(*this);
}

vector<Point> Circle::intersect(const Segment& _cs) const {
	return _cs.intersect(*this);
};

vector<Point> Circle::intersect(const Circle& _cc) const {
	double
			cx1 = this->getCenter().getX(),
			cy1 = this->getCenter().getY(),
			cx2 = _cc.getCenter().getX(),
			cy2 = _cc.getCenter().getY(),
			r1 = this->getRadius(),
			r2 = _cc.getRadius(),
			d = Segment(Point(cx1, cy1), Point(cx2, cy2)).length();
	vector<Point> answer;

	if ((d > fabs(r1 + r2) or (d < fabs(r1 - r2)))) {
//		circles lie separately or one circle is inside the other
		return answer;
	}
	if ((cx1 == cx2) and (cy1 == cy2)) {
		if (r1 == r2) {
//			circles coincide
			answer.emplace_back(Point(this->getRadius() + cx1, cy1));
			answer.emplace_back(Point((-1) * this->getRadius() + cx1, cy1));
			return answer;
		} else {
//			one circle is inside the other
			answer.clear();
			return answer;
		}
	}

	cx2 -= cx1;
	cy2 -= cy1;

	double
			u = (cx2 * cx2 + cy2 * cy2),
			v = (r2 * r2 - r1 * r1),
			root = sqrt(8 * r1 * r1 * cx2 * cx2 - 4 * u * v - pow(u + v, 2)),
			x1 = (u - v + root) / (4 * cx2),
			x2 = (u - v - root) / (4 * cx2),
			y1 = (u - v) / (2 * cx2) + x1,
			y2 = (u - v) / (2 * cx2) + x2;

	answer.emplace_back(Point(x1, y1));
	if ((x1 != x2) and (y1 != y2)) {
		answer.emplace_back(Point(x2, y2));
	}
	return answer;
};

vector<Point> Circle::intersect(const PolyLine& _cp) const {
	vector<Point> answer, tmp, points = _cp.getPoints();
	for (int i = 0; i < points.size() - 1; ++i) {
		tmp = Segment(points[i], points[i + 1]).intersect(*this);
		for (auto j : tmp) {
			answer.push_back(j);
		}
		tmp.clear();
	}
	return answer;
};

vector<Point> Circle::intersect(const Figure& _cf) const {
	return _cf.intersect(*this);
};

Circle::Circle(const Point& c, const double& r) : c_{c}, r_{r} {
	this->recalculateLength();
}

double PolyLine::length() const {
	return this->len_;
}

void PolyLine::setLength() {
	this->len_ = 0;
	for (int i = 0; i < this->points_.size() - 1; ++i) {
		this->len_ += Segment(this->points_[i], this->points_[i + 1]).length();
	}
}

vector<Point> PolyLine::getPoints() const {
	return this->points_;
}

vector<Point> PolyLine::intersect(const Point& _cp) const {
	return _cp.intersect(*this);
}

vector<Point> PolyLine::intersect(const Segment& _cs) const {
	return _cs.intersect(*this);
}

vector<Point> PolyLine::intersect(const Circle& _cc) const {
	return _cc.intersect(*this);
}

vector<Point> PolyLine::intersect(const PolyLine& _cp) const {
	set<Point, comparePoints> tmp_set;
	vector<Point> tmp_vector, points1 = this->getPoints(), points2 = _cp.getPoints(), answer;

	for (auto i = 0; i < points1.size() - 1; ++i) {
		Segment segment = Segment(points1[i], points1[i + 1]);
		for (auto j = 0; j < points2.size() - 1; ++j) {
			tmp_vector = segment.intersect(Segment(points2[j], points2[j + 1]));
			for (const auto& point : tmp_vector) {
				tmp_set.insert(point);
			}
		}
	}

	if (tmp_set.empty()) {
		return answer;
	} else {
		for (auto it : tmp_set) {
			answer.push_back(it);
		}
		return answer;
	}
}

vector<Point> PolyLine::intersect(const Figure& _cf) const {
	return _cf.intersect(*this);
}

PolyLine::PolyLine(const vector<Point>& points) : points_{points} {
	this->setLength();
}