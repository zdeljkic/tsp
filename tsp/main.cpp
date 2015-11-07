#include <iostream>
#include <cmath>
#include <vector>
#include <queue>
#include <cstring>
#include <chrono>
#include <limits>
#include <algorithm>
#include <cstdint>
#include <ctime>
#include "point.h"

/*
 * Algorithm parameters:
 *  > MAXN
 *    - maximum number of cities
 *  > NEAREST
 *    - 2opt should should search NEAREST nearest neighbors of a city for a good move
 *    - recommended range: between 5 and 200
 *    - time vs tour quality tradeoff, less neighbors -> faster but worse quality
 *  > PERCENTAGE
 *    - in greedy, choose the best edge with PERCENTAGE chance, and the second best with 1 - PERCENTAGE
 *    - needed when different starting tours need to be generated
 *    - recommended range: around ??? seems to work well
 *    - lower value -> more variations, but of lesser quality, also slower!
 *  > RANGE
 *    - at what range to jump out of quickselect and use a simple linear algorithm
 *    - recommended range: 2 - 20
 */
#define MAXN 1000
#define NEAREST (N / 5)
#define PERCENTAGE 1.0 - 1.0/20.0
#define RANGE 10

/*
 * Algorithm time parameters:
 *  > TIME_MAX
 *    - maximum time the program is allowed to run
 *  > TIME_2OPT_SAFETY
 *    - if there is still TIME_2OPT_SAFETY time left, do another iteration of 2opt
 *    - otherwise leave and output the result
 *  > TIME_ITER_SAFETY_FAC & _CONST
 *    - iteration safety factor, when iterating the average iteration time is calculated
 *    - if there is still ..._FAC * average_iteration_time + ..._CONST left, do another iteration
 */
#define TIME_MAX 2.0
#define TIME_2OPT_SAFETY 0.1
#define TIME_ITER_SAFETY_FAC 1.2
#define TIME_ITER_SAFETY_CONST 0.25

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
	uint16_t right = N - 1;

	while (right - left > RANGE) {
		uint16_t ind = left + std::rand() % (right - left);

		uint32_t indDist = D[p][near[ind]]; // omfg ko bi se sjetio da se mijenja za vrijeme partitiona
		auto closerThan = [p, indDist, ind, &near](uint16_t a) {
			// If all elements between left and right are equal, we will get stuck in an infinite loop
			// that's why we need a way of ordering elements with equal distance as the index,
			// one that does not give the same result for all elements.
			// It is enough to differentiate between the index and others.
			// Using this specifically will make it so that elements equally distant as the index
			// are "smaller" than the index, while the index itself isn't smaller than itself
			if (D[p][a] == indDist)
				return a != near[ind];
			else
				return D[p][a] < indDist;
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

std::vector<uint16_t> constructRouteGreedy(double p)
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

	_Edge other = pq.top();
	pq.pop();

	for (uint16_t nLeft = N; nLeft > 0; ) {
		_Edge e = pq.top();
		pq.pop();

		// randomization:
		// - with chance p we pick the best edge
		// - and with chance 1-p the second best
		if (((double) std::rand() / RAND_MAX) > p) {
			_Edge tmp = e;
			e = other;
			other = tmp;
		}

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

	static std::vector<std::vector<uint16_t> > near;
	if (near.empty()) {
		for (uint16_t i = 0; i < N; ++i) {
			near.push_back(getNearestQuick(i, NEAREST));
		}
	}

	while (true) {
		int32_t bestImprovement = 0;
		int32_t improvement;
		uint16_t impI, impJ;

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

		if (bestImprovement > 0){
			if (impJ < impI)
				std::swap(impI, impJ);

			auto rb = route.begin() + impI + 1;
			auto re = route.begin() + impJ + 1;
			std::reverse(rb, re);
		} else {
			break;
		}

		if (elapsedTime() > TIME_MAX - TIME_2OPT_SAFETY)
			break;
	};
}

std::vector<uint16_t> iterate2opt()
{
	double avgIterTime = 0;
	uint32_t iterN = 0;

	std::vector<uint16_t> bestRoute;
	uint32_t bestDist = std::numeric_limits<uint32_t>::max();

	while (elapsedTime() + TIME_ITER_SAFETY_FAC * avgIterTime  + TIME_ITER_SAFETY_CONST < TIME_MAX ) {
		double time = elapsedTime();

		std::vector<uint16_t> routeNew(constructRouteGreedy(PERCENTAGE));

		improveRoute2opt(routeNew);

		uint32_t dist = calculateRoute(routeNew);
		++iterN;

		if (dist < bestDist) {
			bestRoute = routeNew;
			bestDist = dist;
//			std::cerr << "naso bolju rutu u iteraciji " << iterN << std::endl;
		}

		time = elapsedTime() - time;
		avgIterTime = ((iterN - 1)*avgIterTime + time) / iterN;
	}

//	std::cerr << "Broj iteracija: " << iterN << std::endl;

	return bestRoute;
}

int main(int argc, char **argv)
{
	// start the timer
	START_TIME = std::chrono::system_clock::now();

	// srand
	//std::srand(std::time(0));

	// parse input arg choosing the algorithm
	char alg = '\0';
	if (argc == 1) {
		alg = 'b';
	} else if (argc == 2) {
		switch (argv[1][0]) {
		case 'n':
		case 'g':
		case 'G':
		case '1':
		case '2':
		case '3':
		case 'N':
			alg = argv[1][0];
			break;
		}
	}

	if (!alg) {
		std::cout << "Usage: " << argv[0] << " [s|alg]\n";
		std::cout << "alg = n|g|G|1|2|3|N, no alg specified = use fastest, don't output distance (for kattis)\n";
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
		route = constructRouteGreedy(0.0);
		break;
	case 'G':
		route = constructRouteGreedy(PERCENTAGE);
		break;
	case '1':
		route = constructRouteGreedy(0.0);
		improveRoute2opt(route);
		break;
	case '2':
		route = constructRouteGreedy(PERCENTAGE);
		improveRoute2opt(route);
		break;
	case '3':
	case 'b':
		route = iterate2opt();
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

    // output the route or distance
    if (alg == 'b')
		outputRoute(route);
    else
		std::cout << "Distance: " << calculateRoute(route) << std::endl;
}
