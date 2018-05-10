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
}

TEST_CASE("Segment's tests", "[]") {
  SECTION("Getters&Setters") {
    Segment s1;
    Point A(42.0, 42.0), B(9.11, 5);
    vector<double> exp_params {37, -32.89, -172.62};
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
  SECTION("Belong") {}
  SECTION("Segment segment intersection 1") {
    Segment s1(Point(0.0, 0.0), Point(2.0, 2.0)), s2(Point(0.0, 2.0), Point(2.0, 0.0));
    Point exp_p(1.0, 1.0);
    int exp_vec_size = 1;
    vector<Point> res = s1.intersect(s2);
    REQUIRE(exp_vec_size == res.size());
    REQUIRE(res[0] == exp_p);
//    REQUIRE(((exp_p.x == res[res.size() - 1].x) && (exp_p.y == res[res.size() - 1].y)));
  }
  Circle circle(Point(3.0, 3.0), 3.0);
  SECTION("Segment circle intersection 1") {
    Segment segment(1.0, -2.0, 8.0, 5.0);
    vector<Point> exp_points{Point(3.0, 0.0), Point(6.0, 3.0)};
    int exp_vec_size = (int) exp_points.size();
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
    int exp_vec_size = (int) exp_points.size();
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

}

TEST_CASE("Circle's tests", "[]") {
  SECTION("Getters&Setters") {}
  SECTION("Length") {}
  SECTION("Belong") {}
}