#include <iostream>
#include <cmath>
#include <vector>
#include <queue>
#include <cstring>
#include <chrono>
#include <limits>
#include <algorithm>
#include <cstdint>
#include "point.h"

#define MAXN 1000
#define NEAREST 100
#define RANGE 10

uint16_t N;
Point P[MAXN];
uint32_t D[MAXN][MAXN];

std::chrono::time_point<std::chrono::system_clock> START_TIME;

void outputRoute(const std::vector<uint16_t> &route)
{
	for (uint16_t i : route) {
		std::cout << i << "\n";
    }
}

double elapsedTime()
{
	std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed = now - START_TIME;

	return elapsed.count();
}

uint32_t calculateRoute(const std::vector<uint16_t> &route)
{
	uint32_t sum = 0;

	for (uint16_t i = 0; i < N; ++i) {
		sum += D[route[i]][route[(i+1)%N]];
	}

	return sum;
}

std::vector<uint16_t> getNearestSimple(uint16_t p, uint16_t n)
{
	bool used[MAXN];
	std::vector<uint16_t> near;

	memset(used, 0, sizeof(bool) * N);
	used[p] = true;
	for (uint16_t k = 0; k < n; ++k) {
		uint16_t minInd;
		uint32_t minDist = std::numeric_limits<uint32_t>::max();

		for (uint16_t i = 0; i < N; ++i) {
			if (!used[i] && D[p][i] < minDist) {
				minInd = i;
				minDist = D[p][i];
			}
		}

		used[minInd] = true;
		near.push_back(minInd);
	}

	return near;
}

std::vector<uint16_t> getNearestQuick(uint16_t p, uint16_t n)
{
	std::vector<uint16_t> near;
	for (uint16_t i = 0; i < N; ++i)
		if (i != p)
			near.push_back(i);

	if (N <= n)
		return near;

	uint16_t left = 0;
	uint16_t right = N;

	while (right - left > RANGE) {
		uint16_t ind = left + std::rand() % (right - left);

		auto closerThan = [p, ind, &near](uint16_t a) {
			return D[p][a] < D[p][near[ind]];
		};
		auto it = std::partition(near.begin() + left, near.begin() + right, closerThan);

		ind = it - near.begin(); // ind == number of smallest elements on the left

		if (ind > n) {
			right = ind;
		} else if (ind < n) {
			left = ind;
		} else {
			left = n;
			break;
		}
	}

	// now, in the vector 'near' we have:
	// - all the smallest elements from the start up until 'left'-1 (exactly 'left' elements)
	// - all the largest elements from 'right' until the end
	// - the middle elements from left to right-1
	// also, left <= n < right
	// and right - left <= RANGE

	// we already have 'left' smallest elements at the start
	// the rest are between 'left' and 'right-1'
	// we get them using a simple algorithm
	for (uint16_t i = left; i < n; ++i) {
		uint16_t minJ = -1;
		uint32_t minDist = std::numeric_limits<uint32_t>::max();

		for (uint16_t j = i; j < right; ++j) {
			if (D[p][near[j]] < minDist) {
				minJ = j;
				minDist = D[p][near[j]];
			}
		}

		std::swap(near[i], near[minJ]);
	}

	// resize 'near' to the first n
	near.resize(n);

	return near;
}

std::vector<uint16_t> constructRouteNN()
{
	std::vector<uint16_t> route;
	bool used[MAXN];

	std::memset(used, 0, sizeof(used));
	uint16_t p = 0;
	for (uint16_t i = 0; i < N; ++i) {
        used[p] = true;
        route.push_back(p);

        uint32_t minD = std::numeric_limits<uint32_t>::max();
        uint16_t minJ = -1;
        for (uint16_t j = 0; j < N; ++j) {
			if (!used[j] && D[p][j] < minD) {
				minJ = j;
				minD = D[p][j];
			}
        }

        p = minJ;
	}

	return route;
}

std::vector<uint16_t> constructRouteGreedy()
{
	struct _Edge
	{
		uint16_t a, b;
		_Edge(uint16_t a, uint16_t b) : a(a), b(b) {}

		bool operator>(const _Edge &other) const
		{
			return D[a][b] > D[other.a][other.b];
		}
	};

	struct _Point
	{
		unsigned char deg;
		uint16_t e[2];
	};

	// fill the priority queue with edges
	std::priority_queue<_Edge, std::vector<_Edge>, std::greater<_Edge> > pq;
	for (uint16_t i = 0; i < N-1; ++i) {
		for (uint16_t j = i+1; j < N; ++j) {
            pq.emplace(i, j);
		}
	}

	_Point points[MAXN];
	std::memset(points, 0, sizeof(_Point) * N);
	for (uint16_t nLeft = N; nLeft > 0; pq.pop()) {
		_Edge e = pq.top();
		_Point &pa = points[e.a];
		_Point &pb = points[e.b];

		if (pa.deg > 1 || pb.deg > 1) {
			continue;
		}

		// check if it will create a cycle of length < N
		// by going from a and trying to get to b
		if (pa.deg == 1 && pb.deg == 1 && nLeft > 1) {
			uint16_t last;
			uint16_t now = e.a;
			uint16_t next = pa.e[0];
			while (points[next].deg == 2) {
				last = now;
				now = next;
				next = points[next].e[0] == last ? points[next].e[1] : points[next].e[0];
			}

			if (next == e.b) {
				continue;
			}
		}

		pa.e[pa.deg++] = e.b;
		pb.e[pb.deg++] = e.a;
		--nLeft;
	}

	// reconstruct path from edges, starting with node 0
	std::vector<uint16_t> route;
	uint16_t last;
	uint16_t now = 0;
	uint16_t next = points[0].e[0];
	for (uint16_t nLeft = N; nLeft > 0; --nLeft) {
		route.push_back(now);

		last = now;
		now = next;
		next = points[next].e[0] == last ? points[next].e[1] : points[next].e[0];
	}

	return route;
}

void improveRoute2opt(std::vector<uint16_t> &route)
{
	struct _Point {
		uint16_t ind;
		uint16_t prev, next;
		uint32_t dprev, dnext;
	} points[MAXN];

//	double DT;
//	double DT2;
//	double DTpre = 0;
//	double DTloop = 0;
//	double DTnear = 0;
//	double DTrev = 0;
//	double DTtot = elapsedTime();

//	DT = elapsedTime();
	std::vector<std::vector<uint16_t> > near;
	for (uint16_t i = 0; i < N; ++i) {
		near.push_back(getNearestQuick(i, NEAREST));
	}
//	DTnear += elapsedTime() - DT;

	while (true) {
		int32_t bestImprovement = 0;
		int32_t improvement;
		uint16_t impI, impJ;

//		DT = elapsedTime();
		for (uint16_t i = 0; i < N; ++i) {
			uint16_t prev = route[ i == 0 ? N-1 : i-1 ];
			uint16_t now = route[i];
			uint16_t next = route[ i == N-1 ? 0 : i+1 ];

			points[now].ind = i;
			points[now].prev = prev;
			points[now].next = next;
			points[now].dprev = D[now][prev];
			points[now].dnext = D[now][next];
		}
//		DTpre += elapsedTime() - DT;

//		DT = elapsedTime();
		for (uint16_t i : route) {
			for (uint16_t n : near[i]) {
				const _Point &pi = points[i];
				const _Point &pn = points[n];

				if (n == pi.prev || n == pi.next)
					continue;

				if (D[i][n] < std::max(pi.dprev, pi.dnext)) {
					improvement = pi.dprev + pn.dprev - D[i][n] - D[pi.prev][pn.prev];
					if (improvement > bestImprovement) {
						bestImprovement = improvement;
						impI = points[pi.prev].ind;
						impJ = points[pn.prev].ind;
					}

					improvement = pi.dnext + pn.dnext - D[i][n] - D[pi.next][pn.next];
					if (improvement > bestImprovement) {
						bestImprovement = improvement;
						impI = pi.ind;
						impJ = pn.ind;
					}
				}
			}
		}
//		DTloop += elapsedTime() - DT;

//		DT = elapsedTime();
		if (bestImprovement > 0){
			if (impJ < impI)
				std::swap(impI, impJ);

			auto rb = route.begin() + impI + 1;
			auto re = route.begin() + impJ + 1;
			std::reverse(rb, re);
		} else {
			break;
		}
//		DTrev += elapsedTime() - DT;

		if (elapsedTime() > 1.9)
			break;
	};

//	std::cerr << "near " << DTnear << std::endl;
//	std::cerr << "pre " << DTpre << std::endl;
//	std::cerr << "loop " << DTloop << std::endl;
//	std::cerr << "rev " << DTrev << std::endl;
//	std::cerr << "total 2opt time " << elapsedTime() - DTtot << std::endl;
}

int main(int argc, char **argv)
{
	// start the timer
	START_TIME = std::chrono::system_clock::now();

	// parse input arg choosing the algorithm
	char alg = '\0';
	if (argc == 1) {
		alg = 'b';
	} else if (argc == 2) {
		switch (argv[1][0]) {
		case 'g':
		case 'n':
		case '2':
		case '3':
		case 'N':
			alg = argv[1][0];
			break;
		}
	}

	if (!alg) {
		std::cout << "Usage: " << argv[0] << " [s|alg]\n";
		std::cout << "alg = n|g|2|3|N, no alg specified = use fastest, don't output distance (for kattis)\n";
		return 1;
	}

	// input the data
	std::cin >> N;
	for (uint16_t i = 0; i < N; ++i) {
		std::cin >> P[i];
	}

	// edge case where N == 1
	if (N == 1) {
		std::cout << "0\n";
		return 0;
	}

	// pre-calculate distances
    for (uint16_t i = 0; i < N; ++i) {
		for (uint16_t j = 0; j < N; ++j) {
			uint32_t d = Point::dist(P[i], P[j]);
			D[i][j] = d;
			D[j][i] = d;
		}
    }

    // calculate the route
    std::vector<uint16_t> route;
    switch (alg) {
    case 'n':
		route = constructRouteNN();
		break;
	case 'g':
		route = constructRouteGreedy();
		break;
	case '2':
	case 'b':
		route = constructRouteGreedy();
		improveRoute2opt(route);
		break;
	case 'N':
		route = constructRouteNN();
		improveRoute2opt(route);
		break;
	default:
		std::cerr << "not implemented yet\n";
		return 2;
		break;
	}

    // output the route
    outputRoute(route);

    // output the distance
    if (alg != 'b')
		std::cout << "Distance: " << calculateRoute(route) << std::endl;
}
