//
// Created by Glushkov Vladimir on 30.03.2018.
//
#include "figures.h"
#include "catch.hpp"

using namespace figures;

TEST_CASE("Segment tests", "[]") {
//	SECTION("Check not intersection methods of class 'Segment'") {
//		SECTION("Check getters, setters and recalculation of parameters") {
//			Segment s(Point(0.0, 0.0), Point(1.0, 1.0));
//			s.setA(Point(12.0, 42.0));
//			REQUIRE(s.getA().getX() == 12.0);
//			REQUIRE(s.getA().getY() == 42.0);
//
//		}
		SECTION("Check segment length") {
			Segment s1(Point(0.0, 0.0), Point(4.0, 3.0));
			REQUIRE((s1.length() - 5.0) < EPS);
			REQUIRE(1 == 2);
		}
//	}
}

