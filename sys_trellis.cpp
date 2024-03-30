#include "conv_enc.h"

void print_link(sdl* link) {
	std::cout << "From  " << link->from << "  to  " << link->to << "  op:  ";
	if (link->R1) {
		std::cout << "1/";
	}
	else { std::cout << "0/"; }
	std::cout << link->op << std::endl;
}

void print_sys_state_diag(sdiag* sys_state_diag) {
	for (int i = 0; i < sys_state_diag->states.size(); i++) {
		std::cout << "state: " << sys_state_diag->states[i].state << "   outputs:" << std::endl;
		print_link(&sys_state_diag->states[i].R1l);
		print_link(&sys_state_diag->states[i].R0l);
	}
}

sdiag make_systematic_sdaig(conv_encoder* encoder) {
	int L = encoder->gen_polys[0].length() - 1; //state length
	sdiag state_diagram;
	std::vector <sdn> states;
	std::vector <sdl> Links;
	for (int i = 0; i< int(pow(2, L)); i++) {
		sdn currstate;
		sdl link;
		currstate.state = n2b(i, L, 2, true);
		//std::cout << currstate.state << std::endl;
		link.from = currstate.state;

		for (int j = 0; j < 2; j++) {
			std::string to;
			to = currstate.state;
			to = i2c(j) + to;
			to.pop_back();
			std::string out = conv_encode_symbol(j, encoder, currstate.state);
			if (j) {
				link.R1 = true;
			}
			else link.R1 = false;
			link.op = out;
			link.to = to;
			//print_link(&link);
			Links.push_back(link);
		}
		currstate.R1l = Links[2*i + 1];
		currstate.R0l = Links[2*i];
		states.push_back(currstate);
	}
	state_diagram.Links = Links;
	state_diagram.states = states;
	return state_diagram;
}

int test_encoding() {
	conv_encoder* test_encoder = new conv_encoder();
	test_encoder->size = 2;
	test_encoder->gen_polys = { b2b("25", 8,2), b2b("37", 8, 2) };
	std::string states;
	

	std::string instr;
	int len = 10;
	for (int i = 0; i < pow(2, len); i++) {
		instr = n2b(i, len, 2, true);
		if (i % int(pow(2, len - 3)) == 0) {
				std::cout << 0.143 * i / pow(2, len - 3) << std::endl;
		}
		std::cout << instr << std::endl;
		std::cout << conv_enc_list(instr, test_encoder) << std::endl;
	}
	std::cout << conv_enc_list(instr, test_encoder) << std::endl;


	delete test_encoder;
	return 1;
}

tl sdl2tl(sdl* sdl_link) {
	tl trellis_link;
	trellis_link.prev = sdl_link->from;
	trellis_link.op = sdl_link->op;
	trellis_link.R1 = sdl_link->R1;
	return trellis_link;
}

void print_trellis_node(tn* node) {
	std::cout << "node state:  " << node->state << "  alpha:  " << node->a << "  beta:  " << node->b << std::endl;
	if (node->edges.size()) {
		for (int i = 0; i < node->edges.size(); i++) {
			std::cout << "from state:   " << node->edges[i].prev << "   ";
			if (node->edges[i].R1) {
				std::cout << "1/";
			}
			else { std::cout << "0/"; }
			std::cout << node->edges[i].op << "   gamma:  " << node->edges[i].gamma << std::endl;
		}
	}
	std::cout << std::endl;
}

void print_trellis_stage(trellis_stage* stage) {
	std::cout << "Trellis stage: " << stage->stage <<"   L(uk|y): " << stage->L<< std::endl;
	for (int i = 0; i < stage->nodes.size(); i++) {
		print_trellis_node(&stage->nodes[i]);
	}
}

void print_trellis(trellis* Trellis) {
	//std::cout << "called ";
	for (int i = 0; i < Trellis->stages.size(); i++) {
		print_trellis_stage(&(Trellis->stages[i]));
	}
}

trellis_stage make_next_stage(trellis_stage* prev_stage, sdiag* state_diag, int s_dash) {
	trellis_stage next_stage;
	next_stage.stage = s_dash;
	next_stage.nodes = {};
	std::vector <std::string> next_states;
	//printf("called : %d\n", prev_stage->nodes.size());
	for (int i = 0; i < prev_stage->nodes.size(); i++) {
		// loop through all nodes of prev stage
		int state = b2n(prev_stage->nodes[i].state, 2);
		//printf("state is : %d\n", state);
		// collect values of all the next states to make next nodes
		sdl R0l = state_diag->states[state].R0l;
		std::string next_R0_state = R0l.to;
		sdl R1l = state_diag->states[state].R1l;
		std::string next_R1_state = R1l.to;
		//std::cout << next_R1_state << " xxxx " << next_R0_state << std::endl;
		//getting the to's correct -- seem to be going to 0 however
		//need to make sure added if not allready acounted for
		if (!(std::find(next_states.begin(), next_states.end(), next_R0_state) != next_states.end())) {
			/* v doesn't contain x */
			next_states.push_back(next_R0_state);
			// trellis_stage nodes need added aswell
			// for some reason had to write out tn
			next_stage.nodes.push_back({ next_R0_state, 0, 0, {} });
		}
		if (!(std::find(next_states.begin(), next_states.end(), next_R1_state) != next_states.end())) {
			/* v doesn't contain x */
			next_states.push_back(next_R1_state);
			next_stage.nodes.push_back({ next_R1_state, 0, 0, {} });
		}
		
		//print_trellis_stage(&next_stage);
		// now need to update the states -- add the related R1 and R0 links
		for (int i = 0; i < next_stage.nodes.size(); i++) {
			if ((next_R0_state == next_stage.nodes[i].state) && 1) {
				// R0 link should be added to matching next_stage node
				tl tl_added = sdl2tl(&R0l);
				next_stage.nodes[i].edges.push_back(tl_added);
			}
			else if(next_R1_state == next_states[i]) {
				tl tl_added = sdl2tl(&R1l);
				next_stage.nodes[i].edges.push_back(tl_added);
			}
			
		}
	}
	//print_trellis_stage(&next_stage);
	return next_stage;
}

trellis make_trellis(int length, sdiag* sys_state_diag) {
	trellis Trellis;
	Trellis.length = length;
	trellis_stage init_stage;
	init_stage.stage = 0;
	tn init_node;
	init_node.state = sys_state_diag->states[0].state;
	// state having 0 size edges vector will identify as root
	init_node.edges = {};
	init_stage.nodes.push_back(init_node);
	Trellis.stages.push_back(init_stage);
	//print_trellis(&Trellis);
	//printf("started");
	for (int i = 1; i < length; i++) {
		trellis_stage next_stage = make_next_stage(&(Trellis.stages[i - 1]), sys_state_diag, i);
		//printf("adding next stage\n");
		//print_trellis_stage(&next_stage);
		Trellis.stages.push_back(next_stage);
		//print_trellis(&Trellis);
	}
	return Trellis;
}

void terminate_trellis(trellis* Trellis) {
	int n = Trellis->stages[1].nodes[0].state.length(); //gives rate as 1/n
	// start at 0 state at end -- then delete all other nodes in that stage
	// add states connected to previous -- delete all other nodes in that stage
	// keep going till length of list is 2^n
	int L = Trellis->length - 1; // used for indexing
	std::vector <std::string> collected = {};
	std::string first_collected_state;
	for (int i = 0; i < n; i++) {
		first_collected_state = '0' + first_collected_state;
	}
	collected.push_back(first_collected_state);
	//std::cout << collected[0] << std::endl;
	//print_trellis_stage(&Trellis->stages[L]);
	// now we want to delete all the nodes not in collected then add the linked nodes to collected
	while (collected.size() < Trellis->stages[L].nodes.size()) {
		for (int i = Trellis->stages[L].nodes.size()-1; i > -1; i--) {
			//std::cout << Trellis->stages[L].nodes.size() << std::endl;
			//going through all nodes
			std::string test_state = Trellis->stages[L].nodes[i].state;
			if (!(std::find(collected.begin(), collected.end(), test_state) != collected.end())) {
				// test_state not in collected => delete node from stage
				//also need to take into account nodes deleted when indexing
				//printf("reached 1 \n");
				Trellis->stages[L].nodes.erase(Trellis->stages[L].nodes.begin() + i);
				//printf("deleted\n");
			}
		}
		//printf("reached 2 \n");
		// now have deleted relevant nodes update collected with linked states
		collected.clear();
		for (int i = 0; i < Trellis->stages[L].nodes.size(); i++) {
			for (int j = 0; j < Trellis->stages[L].nodes[i].edges.size(); j++) {
				//For each edge add the 'prev' element to collected
				collected.push_back(Trellis->stages[L].nodes[i].edges[j].prev);
			}
		}
		L -= 1;
	}
}

int test_trellis() {
	conv_encoder* test_encoder = new conv_encoder();
	test_encoder->size = 2;
	test_encoder->gen_polys = { b2b("15", 8,2), b2b("17", 8, 2) };
	sdiag test_sdiag = make_systematic_sdaig(test_encoder);
	print_sys_state_diag(&test_sdiag);
	trellis test_trellis = make_trellis(6, &test_sdiag);
	print_trellis(&test_trellis);
	return 1;
}