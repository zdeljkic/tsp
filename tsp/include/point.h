#ifndef POINT_H
#define POINT_H

#include <iostream>
#include <cmath>

class Point
{
	public:
		Point(float x = 0, float y = 0);

		static unsigned int dist(const Point &a, const Point &b);

		friend std::istream &operator>>(std::istream &stream, Point &point);
		friend std::ostream &operator<<(std::ostream &stream, const Point &point);

	private:
		float x;
		float y;
};

#endif // POINT_H
