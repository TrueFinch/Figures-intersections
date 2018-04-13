//
// Created by Glushkov Vladimir on 30.03.2018.
//
//#define CATCH_CONFIG_FAST_COMPILE
#include "catch.hpp"
#include "figures.h"

using namespace figures;

TEST_CASE("Segment tests", "[]") {
	SECTION("Check not intersection methods of class 'Segment'") {
		SECTION("Check getters, setters and recalculation of parameters") {
			Segment s(Point(0.0, 0.0), Point(1.0, 1.0));
			s.setA(Point(12.0, 42.0));
			REQUIRE(s.getA().getX() == 12.0);
			REQUIRE(s.getA().getY() == 42.0);

		}
		SECTION("Check segment length") {
			Segment s1(Point(0.0, 0.0), Point(4.0, 3.0));
			REQUIRE((s1.length() - 5.0) < EPS);
//			REQUIRE(1 == 2);
		}
		SECTION("Check parameters") {
			Segment s1(Point(0.0, 0.0), Point(1.0, 1.0)), s2(Point(1.0, 1.0), Point(0.0, 0.0));
			vector<double> v1 = s1.getParameters();
			vector<double> v2 = s2.getParameters();
			REQUIRE(v1.size() == v2.size());
			for (auto i = 0; i < v1.size(); ++i) {
				REQUIRE(v1[i] == ((-1) * v2[i]));
			}
		}
	}
	SECTION("Intersections tests") {
		Segment s1(Point(2.0, 2.0), Point(4.0, 10.0)), s2(Point(2.0, 7.0), Point(4.0, 5.0));
		Segment s3(Point(7.0, 3.0), Point(13.0, 3.0)), s4(Point(8.0, 8.0), Point(8.0, 1.0));
		Segment s5(Point(-7.0, 3.0), Point(-7.0, -4.0)), s6(Point(-7.0, -2.0), Point(-7.0, -7.0));
		Segment s7(Point(-3.0, -7.0), Point(3.0, -7.0)), s8(Point(-3.0, -9.0), Point(3.0, -9.0));
	}
}

