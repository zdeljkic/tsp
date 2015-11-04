#include <random>
#include <iostream>
#include <cstdlib>

#define RANGE 1000000.0f

int main(int argc, char **argv)
{
	unsigned int n;
	if (argc == 1) {
		std::cin >> n;
	} else {
		n = std::atoi(argv[1]);
	}

	std::cout << n << "\n";

	std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dis(0.0f, RANGE);

	for (unsigned int i = 0; i < n; ++i) {
		float rand1 = dis(gen);
		float rand2 = dis(gen);

		std::cout << rand1 << " " << rand2 << "\n";
	}
}
