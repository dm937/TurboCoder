#pragma once
# define I_M_PI           0.3183098862  /* pi */
#include <stdio.h>
#include <malloc.h>
#include <random>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <float.h>
#include <errno.h>
#include <fenv.h>
#include <math.h>  
#include <chrono>
#include <iomanip>
#include <ctime>


struct conv_encoder {
	int size; // output size, R = 1/size
	std::vector <std::string> gen_polys; //there are |size| gen_polys
};


struct sdl {
	bool R1;
	std::string op;
	std::string to;
	std::string from;
};

struct sdn {
	std::string state;
	sdl R1l;
	sdl R0l;
};

struct sdiag {
	std::vector <sdn> states;
	std::vector <sdl> Links;
};

struct tl { 
	//Trellis link
	bool R1; // R1 true then 1/...
	std::string prev; // contains previous state -- as current node has tl in vector
	std::string op; // /op
	double gamma = 0; // gamma
	double AYB = 0.0; // A_k-1 + Y_k + B_k
};

struct tn { 
	// Trellis node and vector of edges connected to node
	std::string state; //string containing state 
	double a = 0; //alpha
	double b = 0; //beta
	std::vector <tl> edges; // vector containing edges into this node
};


struct trellis_stage {
	int stage;
	std::vector <tn> nodes;
	double L = 0; // likelihood -- used by all but 0th stage to generate uk
};

struct trellis {
	int length;
	bool terminated = true;
	std::vector <trellis_stage> stages;
};

struct a_state { // this is just helpful for states
	std::string state;
	int a;
};

struct LLR_Decoder {
	trellis* Trellis; // This will be used for L(uk|y) calculation -- Stores aswell in trellis stages giving hard estimate
	std::vector <double> Luk; //Also used in  L(uk|y) calculation (counts as input)
	std::vector <double> ydec; // This will contain the yk1 and y(1/2)kp used for input -- assign by make_ys(y, &ydec1, &ydec2)
	std::vector <double> Lue; // output used for soft estimate -- Lue = L(uk|y) - Lcyk1 - Luk
	double Ebn0; //Channel value used to calculate Lc
	double Lc; //Assign from Ebn0 when initiating decoder
};

int c2i(char c);

char i2c(int i);

std::string conv_encode_symbol(int insymb, conv_encoder* encoder, std::string state);

std::string conv_enc_list(std::string inlist, conv_encoder* encoder, bool terminated = true);

int b2n(std::string num, int base);

std::string n2b(int num, int length, int base = 2, bool nl = false);

std::string b2b(std::string num, int ib, int nb);

void print_link(sdl* link);

void print_sys_state_diag(sdiag* sys_state_diag);

sdiag make_systematic_sdaig(conv_encoder* encoder);

int test_encoding();

tl sdl2tl(sdl* sdl_link);

void print_trellis_node(tn* node);

void print_trellis_stage(trellis_stage* stage);

void print_trellis(trellis* Trellis);

trellis_stage make_next_stage(trellis_stage* prev_stage, sdiag* state_diag, int s_dash);

trellis make_trellis(int length, sdiag* sys_state_diag);

int test_trellis();

void terminate_trellis(trellis* Trellis);

void print_v(std::vector <double>* v);

double maxs(double a, double b);

double maxsv(std::vector <double> v);

std::vector <double> edge_bin2double(tl* edge);

void compute_y(trellis* Trellis, std::vector <double> y, double Ebn0_dB, std::vector<double> L_ip);

void compute_a(trellis* Trellis);

void compute_b(trellis* Trellis);

void compute_L(trellis* Trellis);

std::string gen_uest(trellis* Trellis);

void BCJR(trellis* Trellis, std::vector <double> y_vec, double Ebn0_dB);

void AWGN(std::vector <double>* x, double sigma);

int gen_sys_results();

int testing_conv_stuff();

std::string rec_conv_enc_list(std::string inlist, conv_encoder* encoder);

std::string rec_conv_symbol_andVk(int i, conv_encoder* encoder, std::string state);

sdiag make_rec_sdaig(conv_encoder* encoder);

void print_vi(std::vector <int>* v);

std::vector <int> make_I(int L);

void interleave(std::vector <double>* v, std::vector<int>* P);

void deinterleave(std::vector <double>* v, std::vector<int>* P);

void print_dec(LLR_Decoder* dec);

double sumv(std::vector <double>* v);

void comp_norm_L(trellis* Trellis);

void comp_norm_b(trellis* Trellis);

void comp_norm_a(trellis* Trellis);

void norm_y_Luk(LLR_Decoder* decoder);

void round_Lue(LLR_Decoder* decoder);

void init_LLR_Decoder(LLR_Decoder* dec, trellis* Trellis, double Ebn0);

void equate_encoder_polys(conv_encoder* encoder);
