#include "conv_enc.h"

void print_v(std::vector <double>* v) {
	for (int i = 0; i < v->size(); i++) {
		std::cout << " " << (*v)[i];
	}
	std::cout << std::endl;
}

// now want to look at implementing BCJR
double maxs(double a, double b) {
	double diff = abs(a - b);
	double corrterm = 0;
	if (diff < 7) corrterm = log(1 + exp(-diff));
	if (a >= b) {
		return a +corrterm;
	}
	return b + corrterm;
}

double maxsv(std::vector <double> v) {
	if (v.size() == 1) {
		return v[0];
	}
	else if (v.size() == 2) {
		return maxs(v[0], v[1]);
	}
	else {
		std::vector <double> v2;
		v2 = std::vector <double>(v.begin() + 1, v.end());
		return maxs(v[0], maxsv(v2));
	}
}

std::vector <double> edge_bin2double(tl* edge) {
	std::string op = edge->op;
	std::vector <double> output;
	for (int i = 0; i < op.length(); i++) {
		output.push_back(2.0 * c2i(op[i]) - 1.0);
	}
	return output;
}

void compute_y(trellis* Trellis, std::vector <double> y, double Ebn0_dB, std::vector<double> L_ip) {
	// Using the log MAP rule
	int n = Trellis->stages[1].nodes[0].edges[0].op.length(); //gives rate as 1/n
	double Rc = 1.0 / (1.0 * n);
	double Lc = 4 * Rc * pow(10, Ebn0_dB/10.0);
	// assuming P(uk = -1) = P(uk= 1) so L(uk) = 0
	// ==> Y(s,s') = Lc/2 sum_{l=1 to n}(x_kl * y_kl)
	// Lc/2 = 4 * Ebno * Rc 
	//stage indexing starts at 1
	if (y.size() != (Trellis->length - 1) * n) {
		printf("non matching trellis and y sizes");
		return;
	}

	for (int L = Trellis->length - 1; L > 0; L--) {
		// slice y for corresponding stage
		std::vector <double> y_slice = std::vector <double>(y.begin() + (L -1)* n, y.begin() + L *n);
		//print_v(&y_slice);
		// go to trellis stage and compute gamma for each edge
		for (int i = 0; i < Trellis->stages[L].nodes.size(); i++) {
			for (int j = 0; j < Trellis->stages[L].nodes[i].edges.size(); j++) {
				std::vector <double> x = edge_bin2double(&Trellis->stages[L].nodes[i].edges[j]);
				//print_v(&x);
				double Gamma = 0.0;
				for (int k = 0; k < n; k++) {
					Gamma += x[k] * y_slice[k];
				}
				Trellis->stages[L].nodes[i].edges[j].gamma = Gamma * 0.5 * Lc;
			}
			//print_trellis_node(&Trellis->stages[L].nodes[i]);
		}
	}
}

void compute_a(trellis* Trellis) {
	// initialise A's 0 if stage = "00..0" -inf otherwise
	int M = int(pow(2, Trellis->stages[0].nodes[0].state.length())); //this gives maximum number of nodes in stage: 2^(#delays)

	for (int L = 0; L < Trellis->stages.size(); L++) {
		//-1 as need to index last node
		for (int i = 0; i < Trellis->stages[L].nodes.size(); i++) {
			if (b2n(Trellis->stages[L].nodes[i].state, 2) == 0) {
				Trellis->stages[L].nodes[i].a = 0;
				Trellis->stages[L].nodes[i].b = 0;
			}
			else {
				Trellis->stages[L].nodes[i].a = -3.4E20; // -inf
				Trellis->stages[L].nodes[i].b = -3.4E20;
			}
		}
	}
	// now need to feed forward through trellis to compute a's
	std::vector <double> prev_as;
	for (int i = 0; i < M; i++) prev_as.push_back(0.0); // just declaring to use indexing
	for (int L = 1; L < Trellis->stages.size(); L++) {
		// first place previous a's into ordered list -- usefull for searching 
		for (int i = 0; i < Trellis->stages[L - 1].nodes.size(); i++) {
			// state of previous node gives index -- corresponding alpha
			prev_as[b2n(Trellis->stages[L - 1].nodes[i].state, 2)] = Trellis->stages[L - 1].nodes[i].a;
		}
		//printf("reached 1 \n");
		for (int i = 0; i < Trellis->stages[L].nodes.size(); i++) {
			std::vector <double> a_plus_gammas;
			for (int j = 0; j < Trellis->stages[L].nodes[i].edges.size(); j++) {
				int prev_state = b2n(Trellis->stages[L].nodes[i].edges[j].prev,2);
				//printf("reached 2, %d\n", prev_state);
				a_plus_gammas.push_back((prev_as[prev_state] + Trellis->stages[L].nodes[i].edges[j].gamma));
			}
			// now have vector of A+Y to use in max*vector()
			//printf("reached 3"); 
			Trellis->stages[L].nodes[i].a = maxsv(a_plus_gammas);
			//print_trellis_stage(&Trellis->stages[L]);
		}
	}
}
//agree with alpha calculations so far


void compute_b(trellis* Trellis) {
	// need to work backward which is bit trickier because it means going back a stage and searching forward -- tl only 1 way info
	int M = int(pow(2, Trellis->stages[0].nodes[0].state.length())); //this gives maximum number of nodes in stage: 2^(#delays)

	//B's already initialised from compute_a
	for (int L = Trellis->stages.size() - 2; L > -1; L--) {
		// make ordered vector of vectors -- containing B+Y for each state
		std::vector <std::vector <double>> BYarr;
		for (int i = 0; i < M; i++) {
			// loading with -inf vectors won't affect max*v as max*(-inf,a) = a + e^-inf = a
			BYarr.push_back({-3.4E20 });
		}
		// now loop through forward stage and collect Y+B indexed using tl.prev
		for (int i = 0; i < Trellis->stages[L + 1].nodes.size(); i++) {
			for (int j = 0; j < Trellis->stages[L + 1].nodes[i].edges.size(); j++) {
				int index = b2n(Trellis->stages[L + 1].nodes[i].edges[j].prev, 2);
				BYarr[index].push_back((Trellis->stages[L + 1].nodes[i].edges[j].gamma + Trellis->stages[L + 1].nodes[i].b));
			}
		}
		//now have collected B+Y going into next nodes
		// use max*v to calculate b's of current row
		for (int i = 0; i < Trellis->stages[L].nodes.size(); i++) {
			//need to use state for indexing
			int index = b2n(Trellis->stages[L].nodes[i].state,2);
			Trellis->stages[L].nodes[i].b = maxsv(BYarr[index]);
		}
		//print_trellis_stage(&Trellis->stages[L]);
	}
}

void compute_L(trellis* Trellis) {
	// L = max*R1[A_k-1(s') + Yk(s',s) + Bk(s)] - max*R0[ as before]
	// work backward -- condense A(s') into ordered list
	// have R1 be 1 list, R2 be another
	// L = maxsv(R1) - maxsv(R0)
	int M = int(pow(2, Trellis->stages[0].nodes[0].state.length())); //this gives maximum number of nodes in stage: 2^(#delays)
	std::vector <double> A_prev;
	for (int i = 0; i < M; i++) {
		A_prev.push_back(0.0);
	}
	std::vector <double> R1;
	std::vector <double> R0;
	for (int L = Trellis->stages.size() - 1; L >0; L--) {
		// go through links and generate AYB for each
		for (int i = 0; i < Trellis->stages[L-1].nodes.size(); i++) {
			// First go thru prev stage and set elements of A
			int index = b2n(Trellis->stages[L - 1].nodes[i].state, 2);
			A_prev[index] = Trellis->stages[L - 1].nodes[i].a;
		}
		//print_v(&A_prev);
		for (int i = 0; i < Trellis->stages[L].nodes.size(); i++) {
			//computing AYB for all edges
			for (int j = 0; j < Trellis->stages[L].nodes[i].edges.size(); j++) {
				double prev_A = A_prev[b2n(Trellis->stages[L].nodes[i].edges[j].prev, 2)];
				Trellis->stages[L].nodes[i].edges[j].AYB = prev_A + Trellis->stages[L].nodes[i].edges[j].gamma + Trellis->stages[L].nodes[i].b;
				// appending to R1 or R0 depending on R1l 
				if (Trellis->stages[L].nodes[i].edges[j].R1) {
					R1.push_back(Trellis->stages[L].nodes[i].edges[j].AYB);
				}
				else {
					R0.push_back(Trellis->stages[L].nodes[i].edges[j].AYB);
				}
			}
		}
		/*
		print_v(&R1);
		printf("next \n");
		print_v(&R0);
		printf("fin\n");
		*/
		if (R1.size() == 0) {
			//printf("R1 zero size\n");
			Trellis->stages[L].L = -3.4E20;
		}
		else if (R0.size() == 0) {
			// never case for systematic but could be for recursive
			Trellis->stages[L].L = 3.4E20;
		}
		else {
			//std::cout << "summation, R1max: " << maxsv(R1) << "  R0max: " << maxsv(R0) << std::endl;
			Trellis->stages[L].L = maxsv(R1) - maxsv(R0);
		}
		R1.clear();
		R0.clear();
	}
}

std::string gen_uest(trellis* Trellis) {
	std::string uest;
	for (int L = 1; L < Trellis->stages.size(); L++) {
		if (Trellis->stages[L].L <= 0) {
			uest.push_back('0');
		}
		else uest.push_back('1');
	}
	return uest;
}

void BCJR(trellis* Trellis, std::vector <double> y_vec, double Ebn0_dB) {
	compute_y(Trellis, y_vec, Ebn0_dB, {});
	compute_a(Trellis);
	compute_b(Trellis);
	compute_L(Trellis);
}


void AWGN(std::vector <double>* x, double sigma) {

	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::normal_distribution<double> distribution(0.0, sigma);
	for (int i = 0; i < x->size(); i++) {
		(*x)[i] += distribution(generator);
	}
}

int gen_sys_results() {
	conv_encoder* test_encoder = new conv_encoder();
	test_encoder->size = 2;
	test_encoder->gen_polys = { b2b("15", 8, 2), b2b("17", 8, 2) };
	sdiag test_sdiag = make_systematic_sdaig(test_encoder);
	std::string u;
	int word_length = 50;
	for (int i = 0; i < word_length; i++) {
		u += '0';
	}
	for (int i = 0; i < test_encoder->gen_polys[0].length() - 1; i++) {
		//u[u.length() - 1 - i] = '0';
	}
	std::string x = conv_enc_list(u, test_encoder, false);
	std::vector <double> x_vec;
	for (int i = 0; i < x.length(); i++) {
		if (x[i] == '1') {
			x_vec.push_back(1.0);
		}
		else {
			x_vec.push_back(-1.0);
		}
	}

	trellis test_trellis = make_trellis(u.size() + 1, &test_sdiag);
	std::cout << u.size() << std::endl;
	std::cout << test_trellis.length << std::endl;
	//terminate_trellis(&test_trellis);
	int error_count = 0;
	std::ofstream myfile;
	std::string filename;
	filename = "15_18_sys_BCJR_tailEbn0.txt";
	myfile.open(filename);
	int N_data_points = 2;
	double Ebn0_min = 3.6;
	double Ebn0_max = 4;
	std::vector <double> Ebs;
	double jump = (Ebn0_max - Ebn0_min) / N_data_points;
	for (int i = 0; i < N_data_points + 1; i++) {
		Ebs.push_back(Ebn0_min);
		Ebn0_min += jump;
	}
	std::string uest;
	for (int e= 0; e<Ebs.size() ; e++){
		int num_runs = pow(10, 8);
		error_count = 0;
		int i = 0;
		while (i < num_runs) {
			if (i % int(num_runs / 100) == 0) std::cout << i * 1.0 / num_runs << "error count  " << (1.0 * error_count) / (1.0 * i * word_length) << "  Ebn0  " << Ebs[e] << std::endl;
			std::vector <double> y_vec = x_vec;
			//print_v(&y_vec);
			AWGN(&y_vec, pow(10, -Ebs[e] / 20.0) / sqrt(2.0));
			//print_v(&y_vec);
			BCJR(&test_trellis, y_vec, Ebs[e]);
			uest = gen_uest(&test_trellis);
			for (int j = 0; j<word_length; j++) {
				if (uest[j] != '0') {
					error_count += 1;
					//printf("error!");
				}
			}
			if (error_count > 150) {
				std::cout << "final error rate :  " << (1.0 * error_count) / (1.0 * i * word_length) << " " << Ebs[e] << std::endl;
				myfile << (1.0 * error_count) / (1.0 * i * word_length) << " " << Ebs[e] << std::endl;
				i = num_runs;
			}
			i++;
		}
	}
	myfile.close();
	return 0;
}

int testing_conv_stuff() {
	conv_encoder* test_encoder = new conv_encoder();
	test_encoder->size = 2;
	test_encoder->gen_polys = { "111", "101"};
	sdiag test_sdiag = make_systematic_sdaig(test_encoder);
	trellis test_trellis = make_trellis(7, &test_sdiag);
	//std::cout << u.size() << std::endl;
	//std::cout << test_trellis.length << std::endl;
	//terminate_trellis(&test_trellis);
	std::vector <double> y = { 0.3, 0.1, -0.5, 0.2, 0.8, 0.5, -0.5, 0.3, 0.1, -0.7, 1.5, -0.4 };
	BCJR(&test_trellis, y, 3.9794);
	print_trellis(&test_trellis);
	std::cout << gen_uest(&test_trellis);
	printf("\n\n");
	y = { 0, 0, 0, 0, -0.8, 0.5, -0.5, -0.3, 0.1, 0.7, 2, -0.4 };
	BCJR(&test_trellis, y, 3.9794);
	print_trellis(&test_trellis);
	std::cout << gen_uest(&test_trellis);
	return 0;
}

