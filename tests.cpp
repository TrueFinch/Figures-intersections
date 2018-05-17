#include <iostream>
#include <vector>
#include "catch.hpp"
#include "figures.h"

using namespace std;
using namespace figures;
TEST_CASE("Point's tests", "[]") {
  SECTION("Operator <") {
    Point p1 = {1.0, 1.0}, p2 = {42.0, 42.0}, p3 = {23.08, 1999.0};
    REQUIRE(p1 < p2);
    REQUIRE(p1 < p3);
    REQUIRE(!(p2 < p3));
    REQUIRE(!(p2 < p1));
    REQUIRE(!(p3 < p1));
    REQUIRE(!(p3 < p2));
  }
  SECTION("Operator ==") {
    Point p1(1.0, 1.0), p2(1.0, 1.0), p3(2.0, 1.0);
    REQUIRE(p1 == p2);
    REQUIRE(!(p1 == p3));
    REQUIRE(!(p2 == p3));
  }
}

TEST_CASE("Segment's tests", "[]") {
  SECTION("Getters&Setters") {
    Segment s1;
    Point A(42.0, 42.0), B(9.11, 5);
    vector<double> exp_params{37, -32.89, -172.62};
    s1.setA(A);
    s1.setB(B);
    vector<double> params = s1.getParameters();
    for (int i = 0; i < params.size(); ++i) {
      REQUIRE(params[i] == exp_params[i]);
    }
  }
  SECTION("Length") {
    Segment s1(Point(-12.0, 0.0), Point(12.0, 0.0)), s2(Point(0.0, 0.0), Point(0.0, 0.0));
    double exp_len1 = 24, exp_len2 = 0;
    REQUIRE(s1.length() == exp_len1);
    REQUIRE(s2.length() == exp_len2);
  }
  SECTION("Belong") {
    Segment segment(7.0, 3.0, 13.0, 3.0);
    vector<Point> belong_points{Point(7.0, 3.0), Point(13.0, 3.0), Point(10.0, 3.0)};
    vector<Point> not_belong_points{Point(0.0, 0.0), Point(200.0, 200.0), Point(6.0, 3.0), Point(20.0, 3.0)};
    for (auto point : belong_points) {
      REQUIRE(segment.belong(point));
    }
    for (auto point : not_belong_points) {
      REQUIRE(!segment.belong(point));
    }
  }
  SECTION("Segment segment intersection 1") {
    Segment s1(Point(0.0, 0.0), Point(2.0, 2.0)), s2(Point(0.0, 2.0), Point(2.0, 0.0));
    Point exp_p(1.0, 1.0);
    int exp_vec_size = 1;
    vector<Point> res = s1.intersect(s2);
    REQUIRE(exp_vec_size == res.size());
    REQUIRE(res[0] == exp_p);
  }
  SECTION("Segment segment intersection 2") {
    Segment s1(Point(6.0, 1.0), Point(6.0, 5.0)), s2(Point(6.0, 3.0), Point(6.0, 7.0));
    vector<Point> exp_points{Point(6.0, 5.0), Point(6.0, 3.0)};
    auto exp_vec_size = (int) exp_points.size();
    vector<Point> res = s1.intersect(s2);
    REQUIRE(res.size() == exp_vec_size);
    for (int i = 0; i < exp_vec_size; ++i) {
      REQUIRE(res[i] == exp_points[i]);
    }
  }
  SECTION("Segment segment intersection 3") {
    Segment s1(Point(2.0, -2.0), Point(4.0, 0.0)), s2(Point(2.0, -3.0), Point(4.0, -3.0));
    int exp_vec_size = 0;
    vector<Point> res = s1.intersect(s2);
    REQUIRE(res.size() == exp_vec_size);
  }
  SECTION("Segment segment intersection 4") {
    Segment s1(Point(-3.0, 1.0), Point(-3.0, 3.0)), s2(Point(-4.0, 3.0), Point(-4.0, 5.0));
    int exp_vec_size = 0;
    vector<Point> res = s1.intersect(s2);
    REQUIRE(res.size() == exp_vec_size);
  }
  SECTION("Segment segment intersection 4") {
    Segment s1(Point(-3.0, 1.0), Point(-3.0, 3.0)), s2(Point(-4.0, 3.0), Point(-4.0, 5.0));
    int exp_vec_size = 0;
    vector<Point> res = s1.intersect(s2);
    REQUIRE(res.size() == exp_vec_size);
  }
  Circle circle(Point(3.0, 3.0), 3.0);
  SECTION("Segment circle intersection 1") {
    Segment segment(1.0, -2.0, 8.0, 5.0);
    vector<Point> exp_points{Point(3.0, 0.0), Point(6.0, 3.0)};
    auto exp_vec_size = (int) exp_points.size();
    vector<Point> res(segment.intersect(circle));
    REQUIRE(res.size() == exp_vec_size);
    for (int i = 0; i < exp_vec_size; ++i) {
      REQUIRE(res[i] == exp_points[i]);
    }
  }
  SECTION("Segment circle intersection 2") {
    Segment segment(4.0, -2.0, 7.0, -2.0);
    int exp_vec_size = 0;
    vector<Point> res(segment.intersect(circle));
    REQUIRE(res.size() == exp_vec_size);
  }
  SECTION("Segment circle intersection 3") {
    Segment segment(3.0, -2.0, 3.0, -4.0);
    int exp_vec_size = 0;
    vector<Point> res(segment.intersect(circle));
    REQUIRE(res.size() == exp_vec_size);
  }
  SECTION("Segment circle intersection 4") {
    Segment segment(-4.0, 3.0, -2.0, 3.0);
    int exp_vec_size = 0;
    vector<Point> res(segment.intersect(circle));
    REQUIRE(res.size() == exp_vec_size);
  }
  SECTION("Segment circle intersection 5") {
    Segment segment(-2.0, 6.0, 3.0, 6.0);
    vector<Point> exp_points{Point(3.0, 6.0)};
    auto exp_vec_size = (int) exp_points.size();
    vector<Point> res(segment.intersect(circle));
    REQUIRE(res.size() == exp_vec_size);
    for (int i = 0; i < exp_vec_size; ++i) {
      REQUIRE(res[i] == exp_points[i]);
    }
  }
  SECTION("Segment circle intersection 6") {
    Segment segment(2.0, 2.0, 4.0, 3.0);
    int exp_vec_size = 0;
    vector<Point> res(segment.intersect(circle));
    REQUIRE(res.size() == exp_vec_size);
  }
  Polyline polyline(vector<Point>{Point(-6, 1), Point(-4, 3), Point(-3, 2), Point(-1, 4), Point(0, 2)});
  SECTION("Segment polyline intersection 1") {
    Segment segment(1.0, 3.0, 3.0, 3.0);
    double exp_vec_size = 0;
    vector<Point> res(segment.intersect(polyline));
    REQUIRE(res.size() == exp_vec_size);
  }
  SECTION("Segment polyline intersection 2") {
    Segment segment(-7.0, 2.0, 3.0, 2.0);
    vector<Point> exp_points{Point(-5.0, 2.0), Point(-3.0, 2.0), Point(0.0, 2.0)};
    auto exp_vec_size = (int) exp_points.size();
    vector<Point> res(segment.intersect(polyline));
    REQUIRE(res.size() == exp_vec_size);
    for (int i = 0; i < exp_vec_size; ++i) {
      REQUIRE(res[i] == exp_points[i]);
    }
  }
  SECTION("Segment polyline intersection 3") {
    Segment segment(-3.0, -2.0, -1.0, 0.0);
    double exp_vec_size = 0;
    vector<Point> res(segment.intersect(polyline));
    REQUIRE(res.size() == exp_vec_size);
  }
}

TEST_CASE("Circle's tests", "[]") {
  SECTION("Getters&Setters") {
    Circle circle;
    circle.setCenter(Point(1.0, 1.0));
    circle.setRadius(42.0);
    Point exp_point(1.0, 1.0);
    REQUIRE(circle.getCenter() == exp_point);
    double exp_r = 42.0;
    REQUIRE(circle.getRadius() == exp_r);
  }
  SECTION("Length") {
    Circle circle(Point(123.0, 123.0), 42.0);
    double exp_len = 2 * M_PI * 42.0;
    REQUIRE(circle.length() == exp_len);
  }
  SECTION("Belong") {
    Circle circle(Point(-1.0, 4.0), 3.0);
    vector<Point> belong_points{Point(-1.0, 7.0), Point(-4.0, 4.0), Point(2.0, 4.0), Point(-1.0, 1.0)};
    vector<Point> not_belong_points{Point(-1.0, 4.0), Point(0.0, 0.0), Point(9.0, 5.0), Point(18.0, 8.0)};
    for (auto point : belong_points) {
      REQUIRE(circle.belong(point));
    }
    for (auto point : not_belong_points) {
      REQUIRE(!circle.belong(point));
    }
  }
  SECTION("Circle circle intersection 1") {
    Circle c1(Point(-3.0, 0.0), 3.0), c2(Point(0.0, 3.0), 3.0);
    vector<Point> exp_points{Point(-3.0, 3.0), Point(0.0, 0.0)};
    auto exp_vec_size = (int) exp_points.size();
    vector<Point> res(c1.intersect(c2));
    REQUIRE(res.size() == exp_vec_size);
    for (int i = 0; i < exp_vec_size; ++i) {
      REQUIRE(res[i] == exp_points[i]);
    }
  }
  SECTION("Circle circle intersection 2") {
    Circle c1(Point(0.0, 3.0), 3.0), c2(Point(6.0, 3.0), 3.0);
    vector<Point> exp_points{Point(3.0, 3.0)};
    auto exp_vec_size = (int) exp_points.size();
    vector<Point> res(c1.intersect(c2));
    REQUIRE(res.size() == exp_vec_size);
    for (int i = 0; i < exp_vec_size; ++i) {
      REQUIRE(res[i] == exp_points[i]);
    }
  }
  SECTION("Circle circle intersection 3") {
    Circle c1(Point(-3.0, 0.0), 3.0), c2(Point(6.0, 3.0), 3.0);
    int exp_vec_size = 0;
    vector<Point> res(c1.intersect(c2));
    REQUIRE(res.size() == exp_vec_size);
  }
  Circle circle(Point(-1.0, 4.0), 3.0);
  SECTION("Circle polyline intersection 1") {
    Polyline polyline(vector<Point>{Point(-5.0, 4.0), Point(-3.0, 4.0), Point(-1.0, 6.0),
                                    Point(-1.0, 8.0), Point(2.0, 8.0), Point(2.0, 2.0),
                                    Point(-1.0, -1.0), Point(-1.0, 2.0)});
    vector<Point> exp_points{Point(-4.0, 4.0), Point(-1, 7.0), Point(2.0, 4.0), Point(-1.0, 1.0)};
    auto exp_vec_size = (int) exp_points.size();
    vector<Point> res(circle.intersect(polyline));
    REQUIRE(res.size() == exp_vec_size);
    for (int i = 0; i < exp_vec_size; ++i) {
      REQUIRE(res[i] == exp_points[i]);
    }
  }
}
TEST_CASE("Polyline's tests", "[]") {
  SECTION("Getter") {
    Polyline polyline(vector<Point>{Point(-4.0, 0.0), Point(-2.0, 2.0), Point(0.0, 0.0), Point(2.0, 2.0),
                                    Point(4.0, 0.0), Point(2.0, -2.0), Point(-2.0, -2.0), Point(-4.0, 0.0)});
    vector<Point> exp_points{Point(-4.0, 0.0), Point(-2.0, 2.0), Point(0.0, 0.0), Point(2.0, 2.0),
                             Point(4.0, 0.0), Point(2.0, -2.0), Point(-2.0, -2.0), Point(-4.0, 0.0)};
    vector<Point> points(polyline.getPoints());
    for (auto i = 0; i < exp_points.size(); ++i) {
      REQUIRE(points[i] == exp_points[i]);
    }
  }
  SECTION("Length") {
    Polyline polyline(vector<Point>{Point(0, 0.0), Point(0.0, 1.0), Point(1.0, 1.0), Point(1.0, 0.0), Point(0.0, 0.0)});
    double exp_length = 4.0;
    REQUIRE(polyline.length() == exp_length);
  }
  SECTION("Belong") {
    Polyline polyline(vector<Point>{Point(0, 0.0), Point(0.0, 1.0), Point(1.0, 1.0), Point(1.0, 0.0), Point(0.0, 0.0)});
    vector<Point> belong_points{Point(0.0, 0.5), Point(0.5, 1.0), Point(1.0, 0.5), Point(0.5, 0.0)};
    vector<Point> not_belong_points{Point(0.5, 0.5), Point(2.0, 1.0), Point(0.9, 0.5), Point(0.5, -0.1)};
    for (auto point : belong_points) {
      REQUIRE(polyline.belong(point));
    }
    for (auto point : not_belong_points) {
      REQUIRE(!polyline.belong(point));
    }
  }
  SECTION("Polyline segment intersection 1") {
    Segment segment(0.0, 0.0, 2.0, 2.0);
    Polyline polyline(vector<Point>{Point(0.0, 4.0), Point(4.0, 0.0), Point(4.0, 2.0), Point(0.0, 2.0)});
    vector<Point> exp_points{Point(2.0, 2.0)};
    auto exp_vec_size = (int) exp_points.size();
    vector<Point> res(polyline.intersect(segment));
    REQUIRE(res.size() == exp_vec_size);
    for (int i = 0; i < exp_vec_size; ++i) {
      REQUIRE(res[i] == exp_points[i]);
    }
  }
  SECTION("Polyline circle intersection 1") {
    Circle circle(Point(0.0, 0.0), 4.0);
    Polyline polyline(vector<Point>{Point(-4.0, 4.0), Point(4.0, 4.0), Point(0.0, 8.0), Point(0.0, 0.0)});
    vector<Point> exp_points{Point(0.0, 4.0)};
    auto exp_vec_size = (int) exp_points.size();
    vector<Point> res(polyline.intersect(circle));
    REQUIRE(res.size() == exp_vec_size);
    for (int i = 0; i < exp_vec_size; ++i) {
      REQUIRE(res[i] == exp_points[i]);
    }
  }
  SECTION("Polyline polyline intersection 1") {
    Polyline p1(vector<Point>{Point(-4.0, 0.0), Point(-2.0, 2.0), Point(0.0, 0.0), Point(2.0, 2.0),
                              Point(4.0, 0.0), Point(2.0, -2.0), Point(-2.0, -2.0), Point(-4.0, 0.0)}),
        p2(vector<Point>{Point(-2.0, 1.0), Point(-2.0, 40), Point(2.0, 4.0), Point(2.0, 1.0), Point(0.0, -1.0)});
    vector<Point> exp_points{Point(-2.0, 2.0), Point(2.0, 2.0)};
    auto exp_vec_size = (int) exp_points.size();
    vector<Point> res(p1.intersect(p2));
    REQUIRE(res.size() == exp_vec_size);
    for (int i = 0; i < exp_vec_size; ++i) {
      REQUIRE(res[i] == exp_points[i]);
    }
  }

}
TEST_CASE("Other's tests", "[]") {
  vector<shared_ptr<Figure>> figures;
  figures.emplace_back((Figure*) new Segment(-1.0, 6.0, -1.0, 8.0)); //make_shared<Figure>(Segment(-1.0, 6.0, -1.0, 8.0))
  figures.emplace_back((Figure*) new Circle(Point(-1.0, 4.0), 3.0));
  figures.emplace_back((Figure*) new Polyline({Point(-1.0, 1.0), Point(-1.0, 7.0), Point(-4.0, 4.0), Point(2.0, 4.0)}));
  vector<Point> exp_points{Point(-1.0, 7.0)};
  auto exp_vec_size = (int) exp_points.size();
  vector<Point> res(figures[0]->intersect(*figures[1]));
  REQUIRE(res.size() == exp_vec_size);
  for (int i = 0; i < exp_vec_size; ++i) {
    REQUIRE(res[i] == exp_points[i]);
  }
  res = figures[1]->intersect(*figures[0]);
  REQUIRE(res.size() == exp_vec_size);
  for (int i = 0; i < exp_vec_size; ++i) {
    REQUIRE(res[i] == exp_points[i]);
  }
  exp_points.clear();
  exp_points.emplace_back(-1.0, 7.0);
  exp_points.emplace_back(Point(-1.0, 1.0));
  exp_points.emplace_back(Point(-4.0, 4.0));
  exp_points.emplace_back(Point(2.0, 4.0));
  exp_vec_size = (int) exp_points.size();
  res = figures[2]->intersect(*figures[1]);
  REQUIRE(res.size() == exp_vec_size);
  for (int i = 0; i < exp_vec_size; ++i) {
    REQUIRE(res[i] == exp_points[i]);
  }
}