#include "conv_enc.h"


std::string rec_conv_enc_list(std::string inlist, conv_encoder* encoder) {
	std::string output;
	std::string state = "";
	for (int i = 0; i < int(encoder->gen_polys[0].length()) - 1; i++) {
		state += '0';
	}
	//std::cout << state << std::endl;
	for (int i = 0; i < inlist.length(); i++) {
		std::string sys_op = rec_conv_symbol_andVk(c2i(inlist[i]), encoder, state);
		state = sys_op[1] + state; // recursive so feedback output into state
		state.pop_back();
		sys_op.pop_back(); // remove Vk so can append to output
		output += sys_op;
	}
	return output;
}

std::string rec_conv_symbol_andVk(int i, conv_encoder* encoder, std::string state) {
	// first encoder->size bits will be the output, last bit
	std::string out;
	out.push_back(i2c(i));
	//  if we treat ip as 0 this will telll us the output of top + going into + with i/p to give Vk
	char top_plus = conv_encode_symbol(0, encoder, state)[0];
	char Vk;
	if (top_plus == i2c(i)) {
		Vk = '0';
	}
	else { 
		Vk = '1'; 
	}
	//assuming working with rate 1/2 coders
	out.push_back(conv_encode_symbol(c2i(Vk), encoder, state)[1]);
	out.push_back(Vk);
	return out;
}

sdiag make_rec_sdaig(conv_encoder* encoder) {
	// should work to make trellis
	int L = encoder->gen_polys[0].length() - 1; //state length
	sdiag state_diagram;
	std::vector <sdn> states;
	std::vector <sdl> Links;
	for (int i = 0; i< int(pow(2, L)); i++) {
		sdn currstate;
		sdl link;
		currstate.state = n2b(i, L, 2, true);
		//std::cout << "check 1 " << currstate.state << std::endl;
		//std::cout << currstate.state << std::endl;
		link.from = currstate.state;

		for (int j = 0; j < 2; j++) {
			// now to will be defined from ip and state as recursive
			std::string out = rec_conv_symbol_andVk(j, encoder, currstate.state); // from this we get recursive
			//std::cout <<"check 2 "<< out << std::endl;
			std::string to;
			to = out[out.size() - 1] + currstate.state; // using Vk
			to.pop_back();
			out.pop_back();
			//std::cout <<"check 3 " <<  out << std::endl;
			if (j) {
				link.R1 = true;
			}
			else link.R1 = false;
			link.op = out;
			link.to = to;
			//print_link(&link);
			Links.push_back(link);
		}
		currstate.R1l = Links[2 * i + 1];
		currstate.R0l = Links[2 * i];
		states.push_back(currstate);
	}
	state_diagram.Links = Links;
	state_diagram.states = states;
	return state_diagram;
}


void print_vi(std::vector <int>* v) {
	for (int i = 0; i < v->size(); i++) {
		std::cout << " " << (*v)[i];
	}
	std::cout << std::endl;
}

double sumv(std::vector <double>* v) {
	double sum = 0;
	for (int i = 0; i < v->size(); i++) {
		sum += (*v)[i];
	}
	return sum;
}

std::vector <int> make_I(int L) {
	std::vector <int> P;
	for (int i = 0; i < L; i++) P.push_back(i);
	std::random_shuffle(P.begin(), P.end());
	//print_vi(&P);
	return P;
}

void interleave(std::vector <double>* v, std::vector<int>* I) {
	if (v->size() != I->size()) {
		std::cout << "permutation error non matching sizes" << std::endl;
		return;
	}
	std::vector <double> placeholder;
	placeholder = *v;
	for (int i = 0; i < v->size(); i++) {
		(*v)[i] = placeholder[(*I)[i]];
	}
}

void interleave_s(std::string* s, std::vector<int>* I) {
	if (s->size() != I->size()) {
		std::cout << "permutation error non matching sizes" << std::endl;
		return;
	}
	std::string placeholder;
	placeholder = *s;
	for (int i = 0; i < s->size(); i++) {
		(*s)[i] = placeholder[(*I)[i]];
	}
}

std::string make_x_trans(std::string* u, conv_encoder* encoder, std::vector <int>*I) {
	// will make x to be transmitted
	// generates random interleave pattern
	if (u->size() % 2 != 0) {
		std::cout << "make_x_trans odd sized u will result more parity bits from non-interleaved decoder" << std::endl;
	}
	std::string x = rec_conv_enc_list(*u, encoder);
	std::string u_interleaved;
	for (int i = 0; i < u->size(); i++) u_interleaved.push_back((*u)[i]);
	//std::cout << u_interleaved << std::endl;
	interleave_s(&u_interleaved, I);
	//std::cout << u_interleaved << std::endl;
	std::string x_interleaved = rec_conv_enc_list(u_interleaved, encoder);
	bool isX1 = true;
	// now puncturing -- collecting info bit then alternating parity bit from non then interleaved
	std::string x_trans;
	for (int i = 0; i < x.size(); i++) {
		if (i % 2 == 0) {
			x_trans.push_back((*u)[i / 2]);
		}
		else if (isX1) {
			x_trans.push_back(x[i]);
		}
		else {
			x_trans.push_back(x_interleaved[i]);
		}
	}
	return x_trans;
}

void deinterleave(std::vector <double>* v, std::vector<int>* I) {
	if (v->size() != I->size()) {
		std::cout << "DEpermutation error non matching sizes" << std::endl;
		return;
	}
	std::vector <double> placeholder;
	placeholder = *v;
	for (int i = 0; i < v->size(); i++) {
		(*v)[(*I)[i]] = placeholder[i];
	}
}


void compute_y_Luk(LLR_Decoder* decoder, double damp_term = 1.0) {
	// Using the log MAP rule
	int n = decoder->Trellis->stages[1].nodes[0].edges[0].op.length(); // gives rate as 1/n
	double Rc = 1.0 / (1.0 * n);
	double Lc = decoder->Lc;
	// assuming P(uk = -1) = P(uk= 1) so L(uk) = 0
	// ==> Y(s,s') = Lc/2 sum_{l=1 to n}(x_kl * y_kl)
	// Lc/2 = 4 * Ebno * Rc 
	//stage indexing starts at 1
	if (decoder->ydec.size() != (decoder->Trellis->stages.size()-1) * n) {
		printf("non matching trellis and y sizes\n");
		return;
	}

	if (decoder->Luk.size() != decoder->Trellis->stages.size() - 1) {
		printf("non matching Luk and trellis sizes\n");
		return;
	}

	for (int L = 1; L < decoder->Trellis->stages.size(); L++) {
		// slice y for corresponding stage
		std::vector <double> y_slice = std::vector <double>(decoder->ydec.begin() + (L - 1) * n, decoder->ydec.begin() + L * n);
		//print_v(&y_slice);
		// go to trellis stage and compute gamma for each edge
		for (int i = 0; i < decoder->Trellis->stages[L].nodes.size(); i++) {
			for (int j = 0; j < decoder->Trellis->stages[L].nodes[i].edges.size(); j++) {
				//std::cout << "edge " << decoder->Trellis->stages[L].nodes[i].edges[j].op << std::endl;
				std::vector <double> x = edge_bin2double(&decoder->Trellis->stages[L].nodes[i].edges[j]);
				//print_v(&x);
				double Gamma = 0.0;
				//TAKE CARE 0 OR 1 !!!!!
				for (int k = 1; k < n; k++) {
					Gamma += x[k] * y_slice[k];
				}
				double uk;
				if (decoder->Trellis->stages[L].nodes[i].edges[j].R1) {
					uk = 1.0;
				}
				else {
					uk = -1.0;
				}
				decoder->Trellis->stages[L].nodes[i].edges[j].gamma = Gamma * 0.5 * Lc + uk * 0.5 * (damp_term * decoder->Luk[L - 1] + y_slice[0] * Lc); //now taking soft information into account
			}
			//print_trellis_node(&Trellis->stages[L].nodes[i]);
		}
	}
}

void make_ys_dec(std::vector <double>* y, std::vector <double>* y1, std::vector <double>* y2, std::vector <int> I) {
	// this will assume puncturing switches intermittently
	// ie P = [1 1 0 0 1 0 0 1] from 3F4 FTR handout
	if ((y->size() != y1->size()) or (y->size() != y2->size())) {
		printf("make_ys incorrect sizing provided");
		return;
	}
	if (I.size() != y->size() / 2) {
		printf("make_ys incorrect I size relative to y");
		return;
	}
	std::vector <double> yk1I;
	for (int i = 0 ; i < y->size(); i++) {
		if (i % 2 == 0) {
			yk1I.push_back((*y)[i]);
		}
	}
	//std::cout << "reached 1" << std::endl;
	interleave(&yk1I, &I);
	bool intoy1 = true;
	//std::cout << "reached 2 " << std::endl;
	for (int i = 0; i < y->size(); i++) {
		//if i is even then we want to copy yk1 into both y1,2
		// if i is odd then initialy want to put ykp into y1 and 0 (erasure) into y2
		// change intoy1 value
		if (i % 2 == 0) {
			(*y1)[i] = (*y)[i];
			// need to interleave y2
			(*y2)[i] = yk1I[i/2];
		}
		else {
			if (intoy1) {
				(*y1)[i] = (*y)[i];
				(*y2)[i] = 0;
				intoy1 = false;
			}
			else {
				(*y2)[i] = (*y)[i];
				(*y1)[i] = 0;
				intoy1 = true;
			}
		}
	}
}


void run_dec(LLR_Decoder* dec, double damp_term = 1.0) {
	// Assumes inputs allready applied	
	// used to update the Lue, Luk|y (trellis) 
	//printf("run_dec \n");
	//print_dec(dec);
	compute_y_Luk(dec, damp_term);
	//print_dec(dec);
	//printf("comp y\n");
	compute_a(dec->Trellis);
	//printf("comp a\n");
	compute_b(dec->Trellis);
	//printf("comp b\n");
	compute_L(dec->Trellis);
	//this will have updated the Luk|y
	// now updating Lue
	for (int i = 0; i < dec->Lue.size(); i++) {
		//std::cout << "Luky  " << dec->Trellis->stages[i + 1].L << "   Luk  " << dec->Luk[i] << "   Lc yk1  " << dec->Lc * dec->ydec[2 * i] << std::endl;
		dec->Lue[i] = dec->Trellis->stages[i+1].L - dec->Luk[i] - dec->Lc * dec->ydec[2*i];
	}
}

void iterate_turbo_decoder(LLR_Decoder* dec1, LLR_Decoder* dec2, std::vector <int>* I, double damp_term = 1.0) {
	std::vector <double> Luk1;
	//round_Lue(dec2);
	Luk1 = dec2->Lue;
	deinterleave(&Luk1, I);
	//print_v(Luk1);
	//std::cout << "reached 1 iter td" << std::endl;
	//print_v(Luk1);
	dec1->Luk = Luk1;
	run_dec(dec1, damp_term);
	//round_Lue(dec1);
	//print_dec(dec1);
	std::vector <double> Luk2;
	for (int i = 0; i < dec1->Luk.size(); i++) {
		Luk2.push_back(dec1->Lue[i]);
	}
	//printf(" 111 resolve issue :: "); print_v(&Luk2);
	interleave(&Luk2, I);
	//printf(" 222 resolve issue :: "); print_v(&Luk2);
	dec2->Luk.clear();
	for (int i = 0; i < Luk2.size(); i++) {
		dec2->Luk.push_back(Luk2[i]);
	}
	//printf("333 res issue ::"); print_v(&dec2->Lue);
	run_dec(dec2, damp_term);
	//soft info will be copied next iteration
	//std::cout << std::endl << "dec 2" << std::endl;
	//print_dec(dec2);
	//std::cout << std::endl << std::endl;
}

double round_up(double value, int decimal_places) {
	const double multiplier = std::pow(10.0, decimal_places);
	return std::ceil(value * multiplier) / multiplier;
}

void round_Lue(LLR_Decoder* decoder) {
	for (int i = 0; i < decoder->Lue.size(); i++) {
		decoder->Lue[i] = round_up(decoder->Lue[i], 2);
	}
}

void norm_y_Luk(LLR_Decoder* decoder) {
	// Using the log MAP rule
	int n = decoder->Trellis->stages[1].nodes[0].edges[0].op.length(); // gives rate as 1/n
	double Rc = 1.0 / (1.0 * n);
	double Lc = decoder->Lc;
	// assuming P(uk = -1) = P(uk= 1) so L(uk) = 0
	// ==> Y(s,s') = Lc/2 sum_{l=1 to n}(x_kl * y_kl)
	// Lc/2 = 4 * Ebno * Rc 
	//stage indexing starts at 1
	if (decoder->ydec.size() != (decoder->Trellis->stages.size() - 1) * n) {
		printf("non matching trellis and y sizes\n");
		return;
	}

	if (decoder->Luk.size() != decoder->Trellis->stages.size() - 1) {
		printf("non matching Luk and trellis sizes\n");
		return;
	}
	//norm_Lus(decoder);
	for (int L = 1; L < decoder->Trellis->stages.size(); L++) {
		// slice y for corresponding stage
		std::vector <double> y_slice = std::vector <double>(decoder->ydec.begin() + (L - 1) * n, decoder->ydec.begin() + L * n);
		//print_v(&y_slice);
		// go to trellis stage and compute gamma for each edge
		for (int i = 0; i < decoder->Trellis->stages[L].nodes.size(); i++) {
			for (int j = 0; j < decoder->Trellis->stages[L].nodes[i].edges.size(); j++) {
				//std::cout << "edge " << decoder->Trellis->stages[L].nodes[i].edges[j].R1 << "/" << decoder->Trellis->stages[L].nodes[i].edges[j].op << std::endl;
				std::vector <double> x = edge_bin2double(&decoder->Trellis->stages[L].nodes[i].edges[j]);
				//print_v(&x);
				//print_v(&y_slice);
				double xi = 0.0;
				//TAKE CARE 0 OR 1 !!!!!
				double xdoty = 0.0;
				for (int k = 0; k < n; k++) {
					if (k != 0) xi += x[k] * y_slice[k];
					xdoty += x[k] * y_slice[k];
				}

				double euclid_dist = 0;
				for (int k = 0; k < n; k++) {
					euclid_dist += pow(x[k] - y_slice[k], 2);
				}
				double uk;
				if (decoder->Trellis->stages[L].nodes[i].edges[j].R1) {
					uk = 1;
					//decoder->Trellis->stages[L].nodes[i].edges[j].gamma = exp(decoder->Luk[L - 1] - decoder->Ebn0 * euclid_dist);
				}
				else {
					uk = -1;
					//decoder->Trellis->stages[L].nodes[i].edges[j].gamma = exp( - decoder->Ebn0 * euclid_dist);
				}
				//std::cout << decoder->Luk[L - 1] << std::endl;
				//print_v(&y_slice);
				// basically dont want numbers growing bigger then 1E20
				//std::cout << decoder->Luk[L - 1] << std::endl;
				//std::cout << xi * 0.5 * Lc + uk * 0.5 * (decoder->Luk[L - 1] + y_slice[0] * Lc) << std::endl;
				double Puk = exp(uk * decoder->Luk[L - 1]) / (1 + exp(uk * decoder->Luk[L - 1]));
				if (((xi * 0.5 * Lc) + uk * 0.5 * (decoder->Luk[L - 1] + y_slice[0] * Lc)) > 100){
					//printf("combat num ins\n");
					decoder->Trellis->stages[L].nodes[i].edges[j].gamma = exp(100);
				}
				if (((xi * 0.5 * Lc) + uk * 0.5 * (decoder->Luk[L - 1] + y_slice[0] * Lc)) < - 100) {
					//printf("combat num ins\n");
					decoder->Trellis->stages[L].nodes[i].edges[j].gamma = exp(-100);
				}
				else {
					decoder->Trellis->stages[L].nodes[i].edges[j].gamma = exp(xi * 0.5 * Lc) * exp(uk * 0.5 * (0.7*decoder->Luk[L - 1] + y_slice[0] * Lc)); //now taking soft information into account
				}
				//std::cout << Puk << std::endl;
				// note Puk more unstable than ukLuk/2 due to very small values or 1 being generated
				//decoder->Trellis->stages[L].nodes[i].edges[j].gamma = Puk * exp(-decoder->Ebn0 * euclid_dist);
				
			}
			//print_trellis_node(&Trellis->stages[L].nodes[i]);
		}
	}
}


void comp_norm_a(trellis* Trellis) {
	// initialise A's 0 if stage = "00..0" -inf otherwise
	int M = int(pow(2, Trellis->stages[0].nodes[0].state.length())); //this gives maximum number of nodes in stage: 2^(#delays)

	for (int L = 0; L < Trellis->stages.size(); L++) {
		//-1 as need to index last node
		for (int i = 0; i < Trellis->stages[L].nodes.size(); i++) {
			if (b2n(Trellis->stages[L].nodes[i].state, 2) == 0) {
				Trellis->stages[L].nodes[i].a = 1;
				Trellis->stages[L].nodes[i].b = 1;
			}
			else {
				Trellis->stages[L].nodes[i].a = 0; // -inf
				Trellis->stages[L].nodes[i].b = 0;
			}
		}
	}
	// now need to feed forward through trellis to compute a's
	std::vector <double> prev_as;
	for (int i = 0; i < M; i++) prev_as.push_back(0.0); // just declaring to use indexing
	for (int L = 1; L < Trellis->stages.size(); L++) {
		double sumaks = 0;
		// first place previous a's into ordered list -- usefull for searching 
		for (int i = 0; i < Trellis->stages[L - 1].nodes.size(); i++) {
			// state of previous node gives index -- corresponding alpha
			prev_as[b2n(Trellis->stages[L - 1].nodes[i].state, 2)] = Trellis->stages[L - 1].nodes[i].a;
		}
		//printf("reached 1 \n");
		for (int i = 0; i < Trellis->stages[L].nodes.size(); i++) {
			std::vector <double> a_plus_gammas;
			for (int j = 0; j < Trellis->stages[L].nodes[i].edges.size(); j++) {
				int prev_state = b2n(Trellis->stages[L].nodes[i].edges[j].prev, 2);
				//printf("reached 2, %d\n", prev_state);
				a_plus_gammas.push_back((prev_as[prev_state] * Trellis->stages[L].nodes[i].edges[j].gamma));
			}
			// now have vector of A+Y to use in max*vector()
			//printf("reached 3"); 
			Trellis->stages[L].nodes[i].a = sumv(&a_plus_gammas);
			sumaks += sumv(&a_plus_gammas);
			//print_trellis_stage(&Trellis->stages[L]);
		}
		for (int i = 0; i < Trellis->stages[L].nodes.size(); i++) {
			// normalising
			Trellis->stages[L].nodes[i].a /= sumaks;
		}
	}
}


void comp_norm_b(trellis* Trellis) {
	// need to work backward which is bit trickier because it means going back a stage and searching forward -- tl only 1 way info
	int M = int(pow(2, Trellis->stages[0].nodes[0].state.length())); //this gives maximum number of nodes in stage: 2^(#delays)

	//B's already initialised from compute_a
	for (int L = Trellis->stages.size() - 2; L > -1; L--) {
		// make ordered vector of vectors -- containing B+Y for each state
		std::vector <double> BYarr;
		for (int i = 0; i < M; i++) {
			// loading with -inf vectors won't affect max*v as max*(-inf,a) = a + e^-inf = a
			BYarr.push_back(0.0);
		}
		// now loop through forward stage and collect Y+B indexed using tl.prev
		for (int i = 0; i < Trellis->stages[L + 1].nodes.size(); i++) {
			for (int j = 0; j < Trellis->stages[L + 1].nodes[i].edges.size(); j++) {
				int index = b2n(Trellis->stages[L + 1].nodes[i].edges[j].prev, 2);
				BYarr[index] += Trellis->stages[L + 1].nodes[i].edges[j].gamma * Trellis->stages[L + 1].nodes[i].b;
			}
		}
		//now have collected B+Y going into next nodes
		// use max*v to calculate b's of current row
		for (int i = 0; i < Trellis->stages[L].nodes.size(); i++) {
			//need to use state for indexing
			int index = b2n(Trellis->stages[L].nodes[i].state, 2);
			Trellis->stages[L].nodes[i].b = BYarr[index] / sumv(&BYarr);
		}
		//print_trellis_stage(&Trellis->stages[L]);
	}
}


void comp_norm_L(trellis* Trellis) {
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
	for (int L = Trellis->stages.size() - 1; L > 0; L--) {
		// go through links and generate AYB for each
		for (int i = 0; i < Trellis->stages[L - 1].nodes.size(); i++) {
			// First go thru prev stage and set elements of A
			int index = b2n(Trellis->stages[L - 1].nodes[i].state, 2);
			A_prev[index] = Trellis->stages[L - 1].nodes[i].a;
		}
		//print_v(&A_prev);
		for (int i = 0; i < Trellis->stages[L].nodes.size(); i++) {
			//computing AYB for all edges
			for (int j = 0; j < Trellis->stages[L].nodes[i].edges.size(); j++) {
				double prev_A = A_prev[b2n(Trellis->stages[L].nodes[i].edges[j].prev, 2)];
				Trellis->stages[L].nodes[i].edges[j].AYB = prev_A * Trellis->stages[L].nodes[i].edges[j].gamma * Trellis->stages[L].nodes[i].b;
				// appending to R1 or R0 depending on R1l 
				if (Trellis->stages[L].nodes[i].edges[j].R1) {
					R1.push_back(Trellis->stages[L].nodes[i].edges[j].AYB);
				}
				else {
					R0.push_back(Trellis->stages[L].nodes[i].edges[j].AYB);
				}
			}
		}
		double sumPk = sumv(&R1) + sumv(&R0);
		for (int i = 0; i < Trellis->stages[L].nodes.size(); i++) {
			//normalising P(s,s',y)
			for (int j = 0; j < Trellis->stages[L].nodes[i].edges.size(); j++) {
				Trellis->stages[L].nodes[i].edges[j].AYB /= sumPk;
				// appending to R1 or R0 depending on R1l 
			}
		}
		double sumR1 = sumv(&R1)/sumPk;
		double sumR0 = sumv(&R0)/sumPk;
		/*
		print_v(&R1);
		printf("next \n");
		print_v(&R0);
		printf("fin\n");
		*/
		if (sumR1 == 0) {
			//printf("R1 zero size\n");
			Trellis->stages[L].L = log(1E-100);
		}
		else if (sumR0 == 0) {
			// never case for systematic but could be for recursive
			Trellis->stages[L].L = log(1E100);
		}
		else if (sumR1/sumR0 < 1E-100) {
			Trellis->stages[L].L = log(1E-100);
		}
		else if (sumR1 / sumR0 > 1E100) {
			Trellis->stages[L].L = log(1E100);
		}
		else {
			//std::cout << "summation, R1max: " << maxsv(R1) << "  R0max: " << maxsv(R0) << std::endl;
			Trellis->stages[L].L = log(sumR1/sumR0);
		}
		R1.clear();
		R0.clear();
	}
}


void print_dec(LLR_Decoder* dec) {
	std::cout << "L(uk|y) hard decisions: " << std::endl;
	for (int i = 1; i < dec->Trellis->stages.size(); i++) {
		std::cout << dec->Trellis->stages[i].L << "  ";
	}
	std::cout << std::endl << "Lu info provided from previous" << std::endl;
	for (int i = 0; i < dec->Luk.size(); i++) {
		std::cout << dec->Luk[i] << "  ";
	}
	int n = dec->Trellis->stages[1].nodes[0].edges[0].op.size();
	std::cout << std::endl << "Lcyk1" << std::endl;
	for (int i = 0; i < dec->ydec.size(); i++) {
		if (i % n == 0) {
			std::cout << dec->ydec[i] << "  ";
		}
	}
	std::cout << std::endl << "Lcykp" << std::endl;
	for (int i = 0; i < dec->ydec.size(); i++) {
		if (i % n != 0) {
			std::cout << dec->ydec[i] << "  ";
		}
	}
	std::cout << std::endl << "Lue  soft info passed to next decoder (de/interleaved AFTER) " << std::endl;
	for (int i = 0; i < dec->Lue.size(); i++) {
		std::cout << dec->Lue[i] << "  ";
	}
	std::cout << std::endl;
}

std::string turbo_uest(LLR_Decoder* dec2, std::vector <int>* I) {
	std::vector <double> Luky_final;
	for (int i = 1; i < dec2->Trellis->stages.size(); i++) {
		Luky_final.push_back(dec2->Trellis->stages[i].L);
	}
	deinterleave(&Luky_final, I);
	//printf(" turbo uest Luk final after deinterleaved Luky of Dec2: ");
	//print_v(&Luky_final);
	std::string uest;
	for (int i = 0; i < Luky_final.size(); i++) {
		if (Luky_final[i] < 0) {
			uest.push_back('0');
		}
		else {
			uest.push_back('1');
		}
	}
	return uest;
}

void init_LLR_Decoder(LLR_Decoder* dec, trellis* Trellis, double Ebn0) {
	//not y still needs to be made and assigned
	dec->Trellis = Trellis;
	dec->Ebn0 = Ebn0;
	std::vector <double> Luk, Lue;
	for (int i = 0; i < Trellis->stages.size() - 1; i++) {
		Luk.push_back(0.0);
		Lue.push_back(0.0);
	}
	for (int i = 0; i < Trellis->stages.size(); i++) {
		Trellis->stages[i].L = 0;
	}
	dec->Lue = Lue;
	dec->Luk = Luk;
	dec->Lc = 4 * pow(10, Ebn0/10) / Trellis->stages[1].nodes[0].edges[0].op.size(); //4R*Ebn0 = 4 Ebn0 / output size -- Ebn0 is linear
}

void clear_trellis(trellis* Trellis) {
	for (int i = 0; i < Trellis->stages.size(); i++) {
		Trellis->stages[i].L = 0;
	}
}

void equate_encoder_polys(conv_encoder* encoder) {
	if (encoder->gen_polys[0].size() != encoder->gen_polys[1].size()) {
		int L0 = encoder->gen_polys[0].size();
		int L1 = encoder->gen_polys[1].size();
		if (L0 < L1) {
			encoder->gen_polys[0] = n2b(b2n(encoder->gen_polys[0], 2), L1, 2, true);
		}
		else encoder->gen_polys[1] = n2b(b2n(encoder->gen_polys[1], 2), L0, 2, true);
	}
}

//call compare theory
int ct(){//double Ebn0, int iterations, int wordlength) {
	//std::ofstream myfile;
	//std::string filename; 
	//filename = "Showing_turbo_strengthening_EBn0_" + std::to_string(Ebn0) + ".txt";
	//myfile.open(filename);
	std::cout << log(1E100);
	double Ebn0 = -3.979;
	std::vector <double> y;
	//for (int i = 0; i <2* wordlength; i++) y.push_back(-1.0);
	//AWGN(&y, pow(10, - Ebn0 / 20.0) / sqrt(2.0));
	y = { 0.3, -4.0, -1.9, -2.0, -2.4, -1.3, 1.2 ,-1.1, 0.7, -2.0 ,-1.0, -2.1, -0.2 ,-1.4, -0.3, -0.1, -1.1 ,0.3 };
	//y.clear();
	//for (int i = 0; i <2* 9; i++) y.push_back(-1.0);
	//AWGN(&y, pow(10, - Ebn0 / 20.0) / sqrt(2.0));
	//std::vector <int> I = make_I(wordlength); // if revert need to -1
	std::vector <int> I = { 1, 4, 7, 2, 5, 9, 3, 6, 8 };
	for (int i = 0; i < I.size(); i++) I[i] -= 1;
	std::vector <double> y1, y2;
	for (int i = 0; i < y.size(); i++) {
		y1.push_back(0.0);
		y2.push_back(0.0);
	}
	make_ys_dec(&y, &y1, &y2, I);
	conv_encoder encoder;
	encoder.size = 2;
	encoder.gen_polys = { "111", "101" };
	equate_encoder_polys(&encoder);
	sdiag state_diagram = make_rec_sdaig(&encoder);
	trellis trell_d1 = make_trellis(y.size() / 2 + 1, &state_diagram);
	trellis trell_d2 = make_trellis(y.size() / 2 + 1, &state_diagram);
	LLR_Decoder dec1;
	LLR_Decoder dec2;
	dec1.ydec = y1;
	dec2.ydec = y2;
	//std::cout << "reach 1 main" << std::endl;
	init_LLR_Decoder(&dec1, &trell_d1, Ebn0);
	init_LLR_Decoder(&dec2, &trell_d2, Ebn0);
	dec1.Lc = 1;
	dec2.Lc = 1;
	//std::cout << "reach 2 main" << std::endl;
	std::vector <double> Luky;
	for (int i = 0; i < 5; i++) {
		iterate_turbo_decoder(&dec1, &dec2, &I);
		std::cout << "dec1" << std::endl;
		print_dec(&dec1);
		/*
		for (int i = 0; i < y.size()/2; i++) {
			// writing the L1(uk|y) first
			std::cout << dec1.Trellis->stages[i + 1].L << " ";
		}
		std::cout << std::endl;
		for (int i = 0; i < y.size() / 2; i++) {
			// next L(uk|y) which is final output 
			Luky.push_back(dec2.Trellis->stages[i + 1].L);
		}
		deinterleave(&Luky, &I);
		for (int i = 0; i < y.size() / 2; i++) {
			// next L(uk|y) which is final output 
			std::cout << Luky[i] << " ";
		}
		Luky.clear();
		//myfile << std::endl;
		*/
		std::cout << std::endl << std::endl << "dec2" << std::endl;
		print_dec(&dec2);
		std::cout << std::endl << "deinterleaved: " << std::endl;
		std::vector <double> deinter = {};
		for (int i = 1; i < dec2.Trellis->stages.size(); i++) {
			deinter.push_back(dec2.Trellis->stages[i].L);
		}
		deinterleave(&deinter, &I);
		print_v(&deinter);
	}
	//myfile.close();
	int error_count=  0;
	for (int i = 1; i < dec2.Trellis->stages.size(); i++) {
		//std::cout << dec2.Trellis->stages[i].L;
		if (dec2.Trellis->stages[i].L > 0.0) {
			error_count += 1;
		}
	}
	/*
	if (error_count > 0) {
		std::cout << "hmmm" << std::endl;
		for (int i = 1; i < dec2.Trellis->stages.size(); i++) {
			std::cout << dec2.Trellis->stages[i].L << " ";
		}
		std::cout << std::endl;
	}
	*/
	return error_count;
}



//call test damping ratio
int test_dr() {
	int N = 500;
	double dr_min_log = -8;
	double dr_max_log = 0;
	std::vector <double> dr;
	for (int i = 0; i < N; i++) dr.push_back(dr_min_log + i * (dr_max_log - dr_min_log) / N);
	for (int i = 0; i < N; i++) dr[i] = pow(10, dr[i] );
	int wordlength = 200;
	int max_it = int(pow(10, 9));
	std::vector <double> Ebs = {4, 8, 12, 20};
	int i = 0;
	int error_count = 0;
	int error_cap = 300;

	std::string u;
	for (int j = 0; j < wordlength; j++) u.push_back('0');

	conv_encoder encoder;
	encoder.size = 2;
	encoder.gen_polys = { b2b("15", 8, 2), b2b("17", 8, 2) };
	equate_encoder_polys(&encoder);
	sdiag state_diagram = make_rec_sdaig(&encoder);

	LLR_Decoder dec1;
	LLR_Decoder dec2;
	trellis trell_d1 = make_trellis(u.size() + 1, &state_diagram);
	trellis trell_d2 = make_trellis(u.size() + 1, &state_diagram);
	
	double min_err_rate = 1.0;
	double damp_term;
	double good_dr = 1.0;
	for (int e = 0; e < Ebs.size(); e++) {
		//std::ofstream myfile;
		std::string filename = "DR_ebn0_" + std::to_string(Ebs[e]) + "_dp_" +std::to_string(N)+ "_dr_" + std::to_string(dr_min_log) + "to" + std::to_string(dr_max_log) + "_wl_" + std::to_string(wordlength) + "_ec_" + std::to_string(error_cap) + ".txt";
		std::cout << filename << std::endl;
		//myfile.open(filename);
		for (int damp_ratio = 0; damp_ratio < dr.size(); damp_ratio++) {
			std::cout << min_err_rate << "   " << good_dr << std::endl;
			i = 0;
			error_count = 0;
			if (damp_ratio % (dr.size() / 100) == 0) std::cout << "\% thru " << damp_ratio * 1.0 / (dr.size() / 100) << std::endl;
			while (i < max_it) {
				//if (i%(max_it/10000) == 0) std::cout << "Ebn0: " << Ebs[e] << "\ % complete " << (i * 100.0) / max_it << " error rate " << error_count * 1.0 / (i * wordlength * 1.0) << " damp ratio : " << dr[damp_ratio] << std::endl;
				std::vector <double> make_u;
				for (int j = 0; j < wordlength; j++) make_u.push_back(0.0);
				AWGN(&make_u, 0.4);
				u = "";
				for (int j = 0; j < wordlength; j++) {
					if (make_u[j] < 0) {
						u.push_back('0');
					}
					else u.push_back('1');
				}
				//std::cout << u << std::endl;
				std::string x = rec_conv_enc_list(u, &encoder);

				//std::cout << "X: " << x << std::endl;
				std::vector <double> x_vec;
				for (int j = 0; j < x.size(); j++) {
					if (x[i] == '0') {
						x_vec.push_back(-1.0);
					}
					else {
						x_vec.push_back(1.0);
					}
				}

				std::vector <int> I = make_I(u.size());
				//print_vi(&I);
				std::vector <double> y1;
				std::vector <double> y2;
				//std::cout << x_trans << std::endl;
				std::vector <double> y_vec;
				for (int it = 0; it < x_vec.size(); it++) {
					y_vec.push_back(x_vec[it]);
				}
				//std::cout << y_vec.size() << std::endl;
				//print_v(&y_vec);
				//AWGN(&y_vec, pow(10, -Ebs[e] / 20.0) / sqrt(2.0));
				//print_v(&y_vec);
				// now need to make y's to go into each channel
				for (int k = 0; k < y_vec.size(); k++) {
					y1.push_back(0.0);
					y2.push_back(0.0);
				}
				make_ys_dec(&y_vec, &y1, &y2, I);
				dec1.ydec = y1;
				dec2.ydec = y2;
				init_LLR_Decoder(&dec1, &trell_d1, Ebs[e]);
				init_LLR_Decoder(&dec2, &trell_d2, Ebs[e]);
				//std::cout <<std::endl << dec1.Lc << "   " << dec2.Lc << std::endl;
				double damp_term = dr[damp_ratio]; //absolute
				for (int iter = 0; iter < 10; iter++) {
					//std::cout << "iterated \n";
					iterate_turbo_decoder(&dec1, &dec2, &I, damp_term);//, damp_term);
					damp_term *= 1.04;
				}
				std::string u_est = turbo_uest(&dec2, &I);

				for (int j = 0; j < u_est.size(); j++) {
					if (u_est[j] != u[j]) {
						error_count += 1;
					}
				}
				if (error_count > error_cap) {
					if ((1.0 * error_count) / (1.0 * i * wordlength) < min_err_rate) {
						min_err_rate = (1.0 * error_count) / (1.0 * i * wordlength);
						good_dr = damp_term;
					}
					//myfile << (1.0 * error_count) / (1.0 * i * wordlength) << " " << dr[damp_ratio] << std::endl;
					if((1.0 * error_count)/ (1.0*i*wordlength) < 0.5) std::cout <<"we got something! "<< (1.0 * error_count) / (1.0 * i * wordlength) << " " << dr[damp_ratio] << std::endl;
					i = max_it;
				}

				i++;
			}
		}
		//myfile.close();
	}
	char a;
	std::cin >> a;
	return 0;
}

//call brown fox
int brownfox() {
	int N = 20;
	double dr_min = 0.001;
	double dr_max = 0.003;
	std::vector <double> dr;
	for (int i = 0; i < N; i++) dr.push_back(dr_min + i * (dr_max - dr_min) / N);
	/*
	int wordlength = 200;
	std::vector <double> make_u;
	for (int i = 0; i < wordlength; i++) make_u.push_back(0.0);
	AWGN(&make_u, 0.4);
	std::string u = "";
	for (int i = 0; i < wordlength; i++) {
		if (make_u[i] < 0) {
			u.push_back('0');
		}
		else u.push_back('1');
	}
	*/
	std::string test = "The quick brown fox jumped over the lazy dog";
	int l = 7;
	std::string u = "";
	for (int i = 0; i < test.size(); i++) {
		u = u + n2b(int(test[i]), l, 2, true);
	}
	int nullpad = 0;//e*40;
	for (int i = 0; i < nullpad; i++) {
		u = '0' + u;
	}
	conv_encoder encoder;
	encoder.size = 2;
	encoder.gen_polys = { b2b("7", 8, 2), b2b("5", 8, 2) };
	equate_encoder_polys(&encoder);
	sdiag state_diagram = make_rec_sdaig(&encoder);

	LLR_Decoder dec1;
	LLR_Decoder dec2;
	trellis trell_d1 = make_trellis(u.size() + 1, &state_diagram);
	trellis trell_d2 = make_trellis(u.size() + 1, &state_diagram);
	std::string x = rec_conv_enc_list(u, &encoder);

	//std::cout << "X: " << x << std::endl;
	std::vector <double> x_vec;
	for (int i = 0; i < x.size(); i++) {
		if (x[i] == '0') {
			x_vec.push_back(-1.0);
		}
		else {
			x_vec.push_back(1.0);
		}
	}
	for (int i = 0; i < dr.size(); i++) {
		std::cout << "Using damping ratio: " << dr[i] << std::endl;
		std::vector <double> Ebs = { 10.0 };
		for (int e = 0; e < 10; e++) {
			//u = "00000000000000000111";
			std::vector <int> I = make_I(u.size());
			//print_vi(&I);
			std::vector <double> y1;
			std::vector <double> y2;
			//std::cout << x_trans << std::endl;
			std::vector <double> y_vec;
			for (int it = 0; it < x_vec.size(); it++) {
				y_vec.push_back(x_vec[it]);
			}
			//print_v(&y_vec);
			//print_v(&x_vec);
			//print_v(&y_vec);
			//std::cout << y_vec.size() << std::endl;
			//print_v(&y_vec);
			AWGN(&y_vec, pow(10, -Ebs[0] / 20.0) / sqrt(2.0));
			//print_v(&y_vec);
			// now need to make y's to go into each channel
			for (int k = 0; k < y_vec.size(); k++) {
				y1.push_back(0.0);
				y2.push_back(0.0);
			}
			make_ys_dec(&y_vec, &y1, &y2, I);
			dec1.ydec = y1;
			dec2.ydec = y2;
			init_LLR_Decoder(&dec1, &trell_d1, Ebs[0]);
			init_LLR_Decoder(&dec2, &trell_d2, Ebs[0]);
			//std::cout <<std::endl << dec1.Lc << "   " << dec2.Lc << std::endl;
			double damp_term = dec1.Lc * dr[i];
			for (int iter = 0; iter < 10; iter++) {
				//std::cout << "iterated \n";
				iterate_turbo_decoder(&dec1, &dec2, &I, damp_term);
				damp_term *= 1.0;
			}
			std::string u_est = turbo_uest(&dec2, &I);
			//below for when null padded
			//u_est.erase(u_est.begin(), u_est.begin() + nullpad);
			std::string toprint = "";

			for (int j = 0; j < u_est.size() / l; j++) {
				std::string block;
				for (int k = 0; k < l; k++) {
					block.push_back(u_est[j * l + k]);
				}
				toprint.push_back(char(b2n(block, 2)));
			}
			//std::cout << u_est << std::endl;
			std::cout << "next red: \n" << toprint << std::endl;
		}
	}
	return 0;
}

//call turbo_res
int turbo_res() {
	// generating some error rate for turbo decoder
	//comment on EbN0 for longer runtimes
	std::vector <int> wls = { 10,20, 50, 100, 200, 500 };
	for (int w = 0; w < wls.size(); w++) {
		conv_encoder encoder;
		encoder.size = 2;
		encoder.gen_polys = { b2b("15", 8, 2), b2b("17", 8, 2) };//{ b2b("15", 8, 2), b2b("17", 8, 2)};
		sdiag state_diagram = make_rec_sdaig(&encoder);

		int word_length = wls[w];
		//generating all 0 word
		std::string u;
		for (int i = 0; i < word_length; i++) {
			u.push_back('0'); // using all 0 codeword for simulations
		}

		LLR_Decoder dec1;
		LLR_Decoder dec2;
		trellis trell_d1 = make_trellis(u.size() + 1, &state_diagram);
		trellis trell_d2 = make_trellis(u.size() + 1, &state_diagram);

		std::string u_est;
		int num_runs = int(pow(10, 6));
		int error_count;
		int N_data_points = 2;
		double Ebn0_min = 0.1;
		double Ebn0_max = 0.2;
		std::string Ebn0_min_s = std::to_string(Ebn0_min);
		std::string Ebn0_max_s = std::to_string(Ebn0_max);
		for (int i = 0; i < 4; i++) {
			Ebn0_min_s.pop_back();
			Ebn0_max_s.pop_back();
		}
		double Ebn0_current = Ebn0_min;
		std::vector <double> Ebs;
		double jump = (Ebn0_max - Ebn0_min) / N_data_points;
		Ebs = { -3, -2, -1, -0.5, 0, 0.15 };
		std::vector <double> x_vec;
		for (int it = 0; it < u.size() * encoder.size; it++) {
			x_vec.push_back(-1.0);
		}
		int error_cap = 300;
		int num_iterations;
		// test 5 first as will have the longest runtimes for smallest Ebs
		std::vector <int> num_iters = {10 };

		for (int iter_index = 0; iter_index < num_iters.size(); iter_index++) {
			num_iterations = num_iters[iter_index];
			std::ofstream myfile;
			std::string filename = "iters_" + std::to_string(num_iterations) + "_ebn0_" + Ebn0_min_s + "to" + Ebn0_max_s + "_wl_" + std::to_string(word_length) + "_ec_" + std::to_string(error_cap) + ".txt";
			std::cout << filename << std::endl;
			myfile.open(filename);
			for (int e = 0; e < Ebs.size(); e++) {
				error_count = 0;
				int i = 0;
				num_runs = int(pow(10, 9)); // cant be bigger than this! -- 200 word length so BER i can observe is 10^-9 (need 100 observations)
				while (i < num_runs) {
					if (i % int(num_runs / 100000) == 0) std::cout << "iters " << num_iterations << " comp \% " << i * 100.0 / num_runs << " error count  " << (1.0 * error_count) / (1.0 * i * word_length) << "  Ebn0  " << Ebs[e] << std::endl;
					// generating random interleaver every time
					std::vector <int> I = make_I(u.size());
					//print_vi(&I);
					std::vector <double> y1;
					std::vector <double> y2;
					//std::cout << x_trans << std::endl;
					std::vector <double> y_vec;
					for (int it = 0; it < x_vec.size(); it++) {
						y_vec.push_back(x_vec[it]);
					}
					//std::cout << y_vec.size() << std::endl;
					//print_v(&y_vec);
					AWGN(&y_vec, pow(10, -Ebs[e] / 20.0) / sqrt(2.0));
					//print_v(&y_vec);
					//print_v(&y_vec);
					// now need to make y's to go into each channel
					for (int k = 0; k < y_vec.size(); k++) {
						y1.push_back(0.0);
						y2.push_back(0.0);
					}
					make_ys_dec(&y_vec, &y1, &y2, I);
					dec1.ydec = y1;
					dec2.ydec = y2;
					init_LLR_Decoder(&dec1, &trell_d1, Ebs[e]);
					init_LLR_Decoder(&dec2, &trell_d2, Ebs[e]);

					for (int iter = 0; iter < num_iterations; iter++) {
						//std::cout << "iterated \n";
						iterate_turbo_decoder(&dec1, &dec2, &I);
					}
					u_est = turbo_uest(&dec2, &I);
					//std::cout << u_est << std::endl;
					//std::cout << u_est;
					for (int j = 0; j < word_length; j++) {
						if (u_est[j] != '0') {
							error_count += 1;
							//printf("error!");
						}
					}
					if (error_count > error_cap) {
						myfile << (1.0 * error_count) / (1.0 * i * word_length) << " " << Ebs[e] << std::endl;
						std::cout << (1.0 * error_count) / (1.0 * i * word_length) << " " << Ebs[e] << std::endl;
						i = num_runs;
					}
					if (i == num_runs - 2) {
						//JIC don't get enough in time
						myfile << (1.0 * error_count) / (1.0 * i * word_length) << " " << Ebs[e] << " premature " << error_count << std::endl;
					}
					i++;
				}
			}
			myfile.close();
		}
	}
	return 1;
	// us generated from dec2 Luky

}

void nnmain() {
	int num_runs = int(pow(10, 4));
	int error_count;
	int N_data_points = 20;
	double Ebn0_min = -1;
	double Ebn0_max = 1;
	std::vector <double> Ebs;
	double jump = (Ebn0_max - Ebn0_min) / N_data_points;
	for (int i = 0; i < N_data_points + 1; i++) {
		Ebs.push_back(Ebn0_min);
		Ebn0_min += jump;
	}
	error_count = 0;
	int num_iter = 5;
	int wordlength = 20;
	for (int e = 0; e < Ebs.size(); e++) {
		for (int i = 0; i < num_runs; i++) {
			if (i % (num_runs / 100) == 0) {
				std::cout << (1.0 * i) / (num_runs) <<"  error rate " << (1.0 * error_count) / (i * 1.0)<< " Ebn0" << Ebs[e] << std::endl;
			}
			//error_count += run_exp(Ebs[e], num_iter, wordlength);
			if (error_count > 100) {
				std::cout << (1.0 * error_count) / (i * 1.0) << std::endl; 
				error_count = 0;
				i = num_runs;
			}
		}
	}
}



int checking_normal() {
	std::vector <double> x;
	int dp = 1000000;
	for (int i = 0; i < dp; i++) x.push_back(0.0);
	AWGN(&x, 0.4);
	int N = 100;
	double x_min = -0.5;
	double x_max = 0.5;
	double jump = (x_max - x_min / N);
	std::vector <double> range;
	for (int i = 0; i < N + 1; i++) {
		range.push_back(x_min + i * jump);
	}
	std::vector <int> bins;
	for (int j = 0; j < dp + 1; j++) bins.push_back(0);

	for (int i = 0; i < dp; i++) {
		for (int j = 0; j < N + 1; j++) {
			if (j == 0) {
				if (x[i] < range[0]) {
					bins[j] += 1;
				}
			}
			else if (j == N) {
				if (x[i] > range[N]) {
					bins[j] += 1;
				}
			}
			else if ((x[i] > range[j]) && (x[i] < range[j + 1])) {
				bins[j] += 1;
			}
		}
	}
	for (int j = 0; j < N; j++) {
		bins[j] /= 700;
		for (int i = 0; i < bins[j]; i++) {
			std::cout << "#";
		}
		if (bins[j] != 0) {
			std::cout << std::endl;
		}
	}
	return 1;
}