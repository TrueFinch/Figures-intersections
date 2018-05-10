#include "figures.h"
#include <cmath>

double det(double a, double b, double c, double d) {
  return a * d - b * c;
}

void swap(double* a, double* b) {
  double buff = *a;
  *a = *b;
  *b = buff;
}

using namespace figures;

Point::Point() {
  this->x = 42.0;
  this->y = 42.0;
}

Point::Point(double x, double y) : x{x}, y{y} {}

const bool Point::operator<(const Point& rhs) const {
  return ((this->x < rhs.x) && (this->y < rhs.y));
}

const bool Point::operator==(const Point& rhs) const {
  return ((this->x == rhs.x) && (this->y == rhs.y));
}

Segment::Segment(const Point& a, const Point& b) : pa_{a}, pb_{b} {
  this->recalculateParameters();
  this->recalculateLength();
}

Segment::Segment(const double& x1, const double& y1, const double& x2, const double& y2) : pa_{x1, y1}, pb_{x2, y2} {
  this->recalculateParameters();
  this->recalculateLength();
}

vector<double> Segment::getParameters() const {
  vector<double> ans;
  ans.push_back(this->a_);
  ans.push_back(this->b_);
  ans.push_back(this->c_);
  return ans;
}

Point Segment::getA() const { return this->pa_; }

Point Segment::getB() const { return this->pb_; }

void Segment::setA(Point& a) {
  this->pa_ = a;
  this->recalculateParameters();
  this->recalculateLength();
}

void Segment::setB(Point& b) {
  this->pb_ = b;
  this->recalculateParameters();
  this->recalculateLength();
}

double Segment::length() const { return this->length_; }

bool Segment::belong(const Point& _cp) const {
  vector<double> p = this->getParameters();
  double x1 = this->getA().x, y1 = this->getA().y, x2 = this->getB().x, y2 = this->getB().y;
  if (x2 > x1) { swap(&x1, &x2); }
  if (y2 > y1) { swap(&y1, &y2); }
  return (((_cp.x * p[0] + _cp.y * p[1] + p[2]) == 0)
      && (x2 <= _cp.x) && (_cp.x <= x1)
      && (y2 <= _cp.y) && (_cp.y <= y1));
}

vector<Point> Segment::intersect(const Figure& _cf) const { return _cf.intersect(*this); }

vector<Point> Segment::intersect(const Segment& _cs) const {
  vector<Point> answer;

  vector<double> p = this->getParameters();
  double a1 = p[0], b1 = p[1], c1 = p[2];
  p = _cs.getParameters();
  double a2 = p[0], b2 = p[1], c2 = p[2];

  double zn1 = det(a1, b1, a2, b2), zn2 = det(a1, c1, a2, c2), zn3 = det(b1, c1, b2, c2);

  if (zn1 == 0) {
    if ((zn2 == 0) and (zn3 == 0)) {
      if (_cs.belong(this->getA())) { answer.push_back(this->getA()); }
      if (_cs.belong(this->getB())) { answer.push_back(this->getB()); }
      if (answer.size() == 2) { return answer; }
      if ((this->belong(_cs.getA()))) { answer.push_back(_cs.getA()); }
      if ((this->belong(_cs.getB()))) { answer.push_back(_cs.getB()); }
    }
    return answer;
  }

  double d = a2 * b1 - a1 * b2, x = (b2 * c1 - b1 * c2) / d, y = (a1 * c2 - a2 * c1) / d;
  Point point = {x, y};
  if (this->belong(point) && _cs.belong(point)) { answer.push_back(point); }
  return answer;
}

vector<Point> Segment::intersect(const figures::Circle& _cc) const {
  vector<Point> answer;
  vector<double> v(this->getParameters());
  double a = v[0], b = v[1], c = v[2], p = _cc.getCenter().x, q = _cc.getCenter().y, r = _cc.getRadius(),
      d = fabs(a * p + b * q + c) / sqrt(a * a + b * b);

  if (d <= r) {
    double x1, y1, x2, y2;
    if (b == 0) {
      y1 = q + sqrt((r - c / a - p) * (r + c / a + p));
      x1 = (-1) * c / a;
      y2 = q - sqrt((r - c / a - p) * (r + c / a + p));
      x2 = (-1) * c / a;
    } else {
      double k = b * b * p - a * (b * q + c), s = a * a + b * b, t = b * sqrt(r * r * s - pow(a * p + b * q + c, 2));
      x1 = (k + t) / s;
      y1 = (-1) * (a * x1 + c) / b;
      x2 = (k - t) / s;
      y2 = (-1) * (a * x2 + c) / b;
    }
    if (d == r) {
      Point point(x1, y1);
      if (this->belong(point)) { answer.push_back(point); }
    } else if (d < r) {
      Point p1(x1, y1), p2(x2, y2);
      if (this->belong(p1)) { answer.push_back(p1); }
      if (this->belong(p2)) { answer.push_back(p2); }
    }
  }
  return answer;
}

vector<Point> Segment::intersect(const figures::Polyline& _cp) const {
  vector<Point> answer, tmp, points{_cp.getPoints()};
  for (int i = 0; i < points.size() - 1; ++i) {
    tmp = this->intersect(Segment(points[i], points[i + 1]));
    for (auto& item : tmp) {
      bool flag = true;
      for (int k = 0; flag && (k < answer.size()); ++k) {
        if (item == answer[k]) { flag = false; }
      }
      if (flag) { answer.push_back(item); }
    }
    tmp.clear();
  }
  return answer;
}

void Segment::recalculateParameters() {
  this->a_ = this->getA().y - this->getB().y;
  this->b_ = this->getB().x - this->getA().x;
  this->c_ = this->getA().x * this->getB().y - this->getB().x * this->getA().y;
  if (this->a_ < 0) {
    this->a_ *= (-1);
    this->b_ *= (-1);
    this->c_ *= (-1);
  } else if ((this->b_ < 0) && (this->a_ == 0)) {
    this->b_ *= (-1);
    this->c_ *= (-1);
  }
}

void Segment::recalculateLength() { this->length_ = sqrt(pow(pb_.x - pa_.x, 2) + pow(pb_.y - pa_.y, 2)); }

Circle::Circle(const Point& c, const double& r) : c_{c}, r_{r} { this->recalculateLength(); }

double Circle::getRadius() const { return r_; }

Point Circle::getCenter() const { return c_; }

void Circle::setRadius(const double& r) {
  r_ = r;
  recalculateLength();
}

void Circle::setCenter(const Point& c) {
  c_ = c;
  recalculateLength();
}

double Circle::length() const { return length_; }

bool Circle::belong(const figures::Point& p) const { return (pow(p.x - c_.x, 2) + pow(p.y - c_.y, 2)) == pow(r_, 2); }

vector<Point> Circle::intersect(const figures::Figure& _cf) const { return _cf.intersect(*this); }

vector<Point> Circle::intersect(const figures::Segment& _cs) const { return _cs.intersect(*this); }

vector<Point> Circle::intersect(const figures::Circle& _cc) const {
  double
      cx1 = getCenter().x, cy1 = getCenter().y, cx2 = _cc.getCenter().x, cy2 = _cc.getCenter().y,
      r1 = getRadius(), r2 = _cc.getRadius(), d = Segment(Point(cx1, cy1), Point(cx2, cy2)).length();
  vector<Point> answer;

  if ((d > fabs(r1 + r2)) or (d < fabs(r1 - r2))) {
//    circles lie separately or one circle is inside the other
    return answer;
  }
  if ((cx1 == cx2) and (cy1 == cy2)) {
    if (r1 == r2) {
//      circles are same to each other
      answer.emplace_back(Point(cx1 + getRadius(), cy1));
      answer.emplace_back(Point(cx1 - getRadius(), cy1));
    }
  } else {
    cx2 -= cx1;
    cy2 -= cy1;
    double
        u = (pow(cx2, 2) + pow(cy2, 2)),
        v = (pow(r2, 2) - pow(r1, 2)),
        root = sqrt(8 * pow(r1, 2) * pow(cx2, 2) - 4 * u * v - pow(u + v, 2)),
        x1 = (u - v + root) / (4 * cx2), y1 = (u - v) / (2 * cx2) + x1,
        x2 = (u - v - root) / (4 * cx2), y2 = (u - v) / (2 * cx2) + x2;
    answer.emplace_back(Point(x1, y1));
    if ((x1 != x2) and (y1 != y2)) { answer.emplace_back(Point(x2, y2)); }
  }
  return answer;
}

vector<Point> Circle::intersect(const figures::Polyline& _cp) const {
  vector<Point> answer, points{_cp.getPoints()};
  for (int i = 0; i < points.size() - 1; ++i) {
    vector<Point> tmp = this->intersect(Segment(points[i], points[i + 1]));
    for (auto& item : tmp) {
      bool flag = true;
      for (int k = 0; flag && (k < answer.size()); ++k) {
        if (item == answer[k]) {
          flag = false;
        }
      }
      if (flag) {
        answer.push_back(item);
      }
    }
  }
  return answer;
}

void Circle::recalculateLength() { length_ = 2 * M_PI * r_; }

Polyline::Polyline(const vector<figures::Point>& points) : points_{points} {
  recalculateLength();
}

vector<Point> Polyline::getPoints() const { return points_; }

double Polyline::length() const { return length_; }

bool Polyline::belong(const figures::Point& point) const {
  bool res;
  for (int i = 0; i < points_.size() - 1; ++i) {
    res = Segment(points_[i], points_[i + 1]).belong(point);
    if (res) {
      return true;
    }
  }
  return false;
}

vector<Point> Polyline::intersect(const Figure& _cf) const { return _cf.intersect(*this); }

vector<Point> Polyline::intersect(const Segment& _cs) const { return _cs.intersect(*this); }

vector<Point> Polyline::intersect(const Circle& _cc) const { return _cc.intersect(*this); }

vector<Point> Polyline::intersect(const Polyline& _cp) const {
  vector<Point> answer, points1(_cp.getPoints()), points2(this->getPoints());
  for (int i = 0; i < points2.size() - 1; ++i) {
    Segment segment(points2[i], points2[i + 1]);
    for (int j = 0; j < points1.size() - 1; ++j) {
      vector<Point> tmp = segment.intersect(Segment(points1[j], points1[j + 1]));
      for (auto& item : tmp) {
        bool flag = true;
        for (int k = 0; flag && (k < answer.size()); ++k) {
          if (item == answer[k]) {
            flag = false;
          }
        }
        if (flag) {
          answer.push_back(item);
        }
      }
    }
//    tmp.clear();
  }
  return answer;
}

void Polyline::recalculateLength() {
  double len = 0.0;
  for (int i = 0; i < points_.size() - 1; ++i) {
    len += Segment(points_[i], points_[i + 1]).length();
  }
  length_ = len;
}