#pragma once

#include "global.h"
#include "util.h"

#include <map>
#include <set>
#include <fstream>
#include <queue>


class BiGraph {
public:
	VI get_12hop(int u) const;
	int is_left_node(int u) const;
	VI get_kcore(int k) const;

	VVI adj_;

	int left_node_num_;
	int right_node_num_;
	int edge_num_;
};

BiGraph* from_text(std::string input_path);
BiGraph* from_edges(std::vector<PII> &edges);


class BiSubgraph {
public:
	BiSubgraph(const BiGraph& g, const VI& node_set);
	BiSubgraph(const BiGraph& g) : g_(g) {}

	void from_node_set(const VI& node_set);
	void clear();
	void deg_rule(int size_bound);
	void remove_nodes(const VI &nodes);

	bool is_node_exist(int u);
	
	const SI& get_remain_adj(int u);
	const SI& get_remain();
	VI intersect_P_with_2hop(const VI& P, int u);

private:
	const BiGraph& g_;
	SI remain_;
	std::vector<SI> remain_adj_;
};
