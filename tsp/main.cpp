#include <iostream>
#include <cmath>
#include <vector>
#include <queue>
#include <cstring>
#include <chrono>
#include <limits>
#include "point.h"

#define MAXN 1000

unsigned short N;
Point P[MAXN];
unsigned int D[MAXN][MAXN];

std::chrono::time_point<std::chrono::system_clock> START_TIME;

double elapsedTime()
{
	std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed = now - START_TIME;

	return elapsed.count();
}

unsigned int calculateRoute(const std::vector<unsigned short> &route)
{
	unsigned int sum = 0;

	for (unsigned short i = 0; i < N; ++i) {
		sum += D[route[i]][route[(i+1)%N]];
	}

	return sum;
}

std::vector<unsigned short> constructRouteNN()
{
	std::vector<unsigned short> route;
	bool used[MAXN];

	std::memset(used, 0, sizeof(used));
	unsigned short p = 0;
	for (unsigned short i = 0; i < N; ++i) {
        used[p] = true;
        route.push_back(p);

        unsigned int minD = std::numeric_limits<unsigned int>::max();
        unsigned short minJ = -1;
        for (unsigned short j = 0; j < N; ++j) {
			if (!used[j] && D[p][j] < minD) {
				minJ = j;
				minD = D[p][j];
			}
        }

        p = minJ;
	}

	return route;
}

std::vector<unsigned short> constructRouteGreedy()
{
	struct _Edge
	{
		unsigned short a, b;
		_Edge(unsigned short a, unsigned short b) : a(a), b(b) {}

		bool operator>(const _Edge &other) const
		{
			return D[a][b] > D[other.a][other.b];
		}
	};

	struct _Point
	{
		unsigned char deg;
		unsigned short e[2];
	};

	//std::cerr << "ulazim u greedy: " << elapsedTime() << std::endl; //DEBUG

	// fill the priority queue with edges
	std::priority_queue<_Edge, std::vector<_Edge>, std::greater<_Edge> > pq;
	for (unsigned short i = 0; i < N-1; ++i) {
		for (unsigned short j = i+1; j < N; ++j) {
            pq.emplace(i, j);
		}
	}

//	std::cerr << "napunio priority queue: " << elapsedTime() << std::endl; //DEBUG
//	std::cerr << "unutra ih je: " << pq.size() << std::endl; //DEBUG
//	double DEB_CTOT = 0.0;
//	double DEB_TRAZ = elapsedTime();

	_Point points[MAXN];
	std::memset(points, 0, sizeof(points));
	for (unsigned short nLeft = N; nLeft > 0; pq.pop()) {
		_Edge e = pq.top();
		_Point &pa = points[e.a];
		_Point &pb = points[e.b];
		//std::cerr << e.a << " " << e.b << " " << D[e.a][e.b];

		if (pa.deg > 1 || pb.deg > 1) {
			//std::cerr << " NE - DEG" << std::endl;
			continue;
		}

		// check if it will create a cycle of length < N
		// by going from a and trying to get to b
		if (pa.deg == 1 && pb.deg == 1 && nLeft > 1) {
			//std::cerr << "--trazim ciklus: " << elapsedTime() << std::endl; //DEBUG
//			double DEB_C = elapsedTime();//DEBUG
			unsigned short last;
			unsigned short now = e.a;
			unsigned short next = pa.e[0];
			while (points[next].deg == 2) {
				last = now;
				now = next;
				next = points[next].e[0] == last ? points[next].e[1] : points[next].e[0];
			}
//			DEB_CTOT += elapsedTime() - DEB_C;

			if (next == e.b) {
				//std::cerr << " NE - CIKLUS" << std::endl;
				continue;
			}
			//std::cerr << "--naso ciklus: " << elapsedTime() << std::endl; //DEBUG
		}

		pa.e[pa.deg++] = e.b;
		pb.e[pb.deg++] = e.a;
		--nLeft;

		//std::cerr << " DA " << nLeft << std::endl;
	}

//	std::cout << "naso rutu: " << elapsedTime() << std::endl; //DEBUG
//	std::cout << "vrijeme u ciklusima: " << DEB_CTOT << std::endl; //DEBUG
//	std::cout << "vrijeme trazenja: " << elapsedTime() - DEB_TRAZ << std::endl; //DEBUG

	// reconstruct path from edges, starting with node 0
	std::vector<unsigned short> route;
	unsigned short last;
	unsigned short now = 0;
	unsigned short next = points[0].e[0];
	for (unsigned short nLeft = N; nLeft > 0; --nLeft) {
		route.push_back(now);

		last = now;
		now = next;
		next = points[next].e[0] == last ? points[next].e[1] : points[next].e[0];
	}

//	std::cout << "rekonstruiro rutu: " << elapsedTime() << std::endl; //DEBUG

	return route;
}

void improveRoute2opt(std::vector<unsigned short> &route)
{
	for (unsigned short i : route) {
		for (unsigned short i : route) {
			std::cout << i;
		}
	}
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
			alg = argv[1][0];
			break;
		}
	}

	if (!alg) {
		std::cout << "Usage: " << argv[0] << " [s|alg]\n";
		std::cout << "alg = n|g|2|3, no alg specified = use fastest\n";
		return 1;
	}

	// input the data
	std::cin >> N;
	for (unsigned short i = 0; i < N; ++i) {
		std::cin >> P[i];
	}

	// edge case where N == 1
	if (N == 1) {
		std::cout << "0\n";
		return 0;
	}

	// pre-calculate distances
    for (unsigned short i = 0; i < N; ++i) {
		for (unsigned short j = 0; j < N; ++j) {
			unsigned int d = Point::dist(P[i], P[j]);
			D[i][j] = d;
			D[j][i] = d;
		}
    }

    // calculate the route
    std::vector<unsigned short> route;
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

		break;
	default:
		std::cout << "not implemented yet\n";
		return 2;
		break;
	}
    // output the route
    for (unsigned short i : route) {
		std::cout << i << "\n";
    }

    // output the distance
    if (alg != 'b')
		std::cout << "Distance: " << calculateRoute(route) << std::endl;

    // output the time
    //std::cout << elapsedTime() << std::endl;
}
