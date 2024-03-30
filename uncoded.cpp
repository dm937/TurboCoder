#include "conv_enc.h"

int main() {
	int wl = 200;
	std::string u_in = "";
	std::vector <double> x_vec;
	for (int i = 0; i < wl; i++) x_vec.push_back(-1.0);
	int max_it = int(pow(10, 9));
	std::vector <double> Ebs = {-3, -2, -1, 0, 1, 2, 3};
	int error_cap = 1000;
	std::ofstream myfile;
	std::string filename = "uncoded_" + std::to_string(Ebs[0]) + "_to_" +std::to_string(Ebs[Ebs.size()-1]) + "_wl_" + std::to_string(wl) + ".txt";
	std::cout << filename << std::endl;
	myfile.open(filename);

	for (int e = 0; e < Ebs.size(); e++) {
		int error_count = 0;
		int i = 0;
		while (i < max_it) {
			if (i % int(max_it / 100000) == 0) std::cout << "iters " << max_it << " comp \% " << i * 100.0 / max_it << " error count  " << (1.0 * error_count) / (1.0 * i * wl) << "  Ebn0  " << Ebs[e] << std::endl;
			std::vector <double> y_vec = x_vec;
			AWGN(&y_vec, pow(10, -Ebs[e] / 20.0) / sqrt(2.0));
			//print_v(&y_vec);
			std::vector <double> result;
			for (int j = 0; j < y_vec.size(); j++) {
				if (y_vec[j] < 0) result.push_back(0);
				else result.push_back(1.0);
			}
			//print_v(&result);
			error_count += sumv(&result);
			if (error_count > error_cap) {
				myfile << (1.0 * error_count) / (1.0 * i * wl) << " " << Ebs[e] << std::endl;
				std::cout << (1.0 * error_count) / (1.0 * i * wl) << " " << Ebs[e] << std::endl;
				i = max_it;
			}
			if (i == max_it - 2) {
				//JIC don't get enough in time
				myfile << (1.0 * error_count) / (1.0 * i * wl) << " " << Ebs[e] << " premature " << error_count << std::endl;
			}
			i++;
		}
	}
	myfile.close();
	return 1;
}