#include "point.h"

Point::Point(float x, float y) :
x(x),
y(y)
{
}

uint32_t Point::dist(const Point& a, const Point& b)
{
	float dx2 = (a.x - b.x) * (a.x - b.x);
	float dy2 = (a.y - b.y) * (a.y - b.y);

	return (uint32_t) (std::sqrt(dx2 + dy2) + 0.5f);
}

std::istream &operator>>(std::istream &stream, Point &point)
{
	stream >> point.x >> point.y;

	return stream;
}

std::ostream &operator<<(std::ostream &stream, const Point &point)
{
	stream << point.x << " " << point.y;

	return stream;
}
