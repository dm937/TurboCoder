#include "conv_enc.h"

//I think I'll keep binary numbers as strings


int c2i(char c) {
	return (int(c) - int('0'));
}

char i2c(int i) {
	if (i > 9) {
		std::cout << "too large int i2c" << std::endl;
	}
	return char(i + int('0'));
}

//function that given binary generators will generate outputs
std::string conv_encode_symbol(int insymb, conv_encoder* encoder, std::string state) {
	std::string out = "";
	int tot = 0;
	int state_len = encoder->gen_polys[0].length();
	for (int i = 0; i < encoder->size; i++) {
		for (int j = 0; j < state_len; j++) {
			int curr = c2i(encoder->gen_polys[i][j]);
			if (j == 0) {
				tot += curr * insymb;
			}
			else {
				tot += curr * c2i(state[j-1]);
			}
		}
		//std::cout <<std::endl << "total is " << tot << std::endl;
		out += char(tot % 2 + int('0'));
		tot = 0;
	}
	//std::cout << out;
	return out;
}

std::string conv_enc_list(std::string inlist, conv_encoder* encoder, bool terminated) {
	std::string output;
	std::string state = "";
	for (int i = 0; i < int(encoder->gen_polys[0].length())-1; i++) {
		state += '0';
	}
	//std::cout << state << std::endl;
	for (int i = 0; i < inlist.length(); i++) {
		output += conv_encode_symbol(c2i(inlist[i]), encoder, state);
		state = inlist[i] + state;
		state.pop_back();
	}
	if (terminated) {
		for (int i = 0; i < state.length(); i++) {
			output += conv_encode_symbol(0, encoder, state);
			state = '0' + state;
			state.pop_back();
		}
	}
	return output;
}

int b2n(std::string num, int base) {
	int number = 0;
	int q = 1;
	for (int i = 0; i < num.length(); i++) {
		number += c2i(num[num.length() - i - 1]) * q;
		q *= base;
	}
	return number;
}

std::string n2b(int num, int length, int base, bool nl) {
	std::string reverse;
	std::string Bnum;

	if (num == 0 && nl) {
		for (int i = 0; i < length; i++) {
			Bnum += '0';
		}
		return Bnum;
	}

	while (num) {
		reverse += std::to_string(num % base);
		num /= base;
	}

	if (nl == true) {
		int diff = length - reverse.length();
		if (diff < 0) {
			std::cout << "incorrect length";
		}
		else {
			for (int i = 0; i < diff; i++) {
				reverse += '0';
			}
		}
	}

	for (int i = 0; i < reverse.length(); i++) {
		Bnum += reverse[reverse.length() - i - 1];
	}
	
	return Bnum;
}

std::string b2b(std::string num, int ib, int nb) {
	int number = b2n(num, ib);
	return n2b(number, nb);
}


