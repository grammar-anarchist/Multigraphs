#include <iostream>
#include "../z2_homologies.hpp"

int main() {
	Filtration f;
	int simpleces_num;
	std::cin >> simpleces_num;
	int dim, dim_num;
	for (int i = 0; i != simpleces_num; i += dim_num) {
		std::cin >> dim >> dim_num;
		int curr;
		double time;
		for (int j = 0; j != dim_num; ++j) {
			std::vector<int> vertices;
			for (int k = 0; k != dim; ++k) {
				std::cin >> curr;
				vertices.push_back(curr);
			}
			std::cin >> time;
			f.push(Simplex{vertices, time, (int)vertices.size() - 1});
		}
	}
	auto ans = f.persistent_homologies();
	for (size_t i = 0; i != ans.size(); ++i) {
		std::cout << "dim" << i << "\n";
		for (size_t j = 0; j != ans[i].size(); ++j) {
			std::cout << f[ans[i][j].ind_birth] << " ";
			double inf = std::numeric_limits<double>::infinity();
			if (ans[i][j].death != inf) {
				std::cout << f[ans[i][j].ind_death] << "\n";
			} else {
				std::cout <<"immortal\n";
			}
		}
	}
	std::cout << "finished\n";
}
