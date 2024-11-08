#pragma once

#include "bigraph.h"

#include <map>


struct SimilarBiclique {
public:
	VI L_;
	VI R_;
	double sim_;

	SimilarBiclique(VI& L, VI &R, double sim) : L_(L), R_(R), sim_(sim) {}
	SimilarBiclique() : L_(VI()), R_(VI()), sim_(0) {}

	void sort_vertices();
};

class SimilarityStore {
public:
	SimilarityStore(const VVI &adj) : adj_(adj) {}
	double get_sim(int u, int v);
	void set_sim(int u, int v, double sim);

private:
	const VVI &adj_;
	std::map<PII, double> sim_;
};

class LocalExact {
public:
	LocalExact(const BiGraph& g) : g_(g), sim_(g.adj_), size_(-1), subgraph_(g), cnt_(g.left_node_num_, 0) {}
	virtual SimilarBiclique local_exact_query(int q, int size, double init_sim = -1.0);
	SimilarBiclique local_exact_query_left(int q, int size, const VI& remain_left, double init_sim = -1.0);

protected:
	virtual void _sim_rule(int q, double sim);
	virtual void _deg_rule();
	virtual void _enum(const SimilarBiclique &cur, const VI& P, std::map<int, double> &last_sim);
	virtual void _simiarity_first_search_on_subgraph(int q);

	const BiGraph& g_;
	SimilarityStore sim_;
	VI cnt_;

	// temporary result, initialize upon every query
	SimilarBiclique result_;
	BiSubgraph subgraph_;
	int size_;

	VI _get_2hop_and_sim(int q, double sim);
	VI _get_2hop_and_sim(int q, const VI&remain_left, double sim);
};

class GlobalExact{
public:
	GlobalExact(const BiGraph& g) : g_(g) {}
	SimilarBiclique global_exact_query(int size);

private:
	const BiGraph& g_;
};

class GlobalApp {
public:
	GlobalApp(const BiGraph& g, int hash_num) : g_(g) { __init_min_hash(hash_num); }
	SimilarBiclique global_app_query(int size);

private:
	void __init_min_hash(int hash_num);

	const BiGraph& g_;
	int hash_num_;
	VVI min_hash_;
};

/// <summary>
/// Local Exact Abalation Test Classes
/// </summary>
class LocalExactNoDeg : public LocalExact {
public:
	LocalExactNoDeg(const BiGraph& g) : LocalExact(g) {}
protected:
	void _deg_rule() {
		return;
	}
};

class LocalExactNoHop : public LocalExact {
public:
	LocalExactNoHop(const BiGraph& g) : LocalExact(g) {}
	SimilarBiclique local_exact_query(int q, int size, double init_sim = -1.0);
	VI get_sim(int q, double sim);
};


class LocalExactNoSim : public LocalExact {
public:
	LocalExactNoSim(const BiGraph& g) : LocalExact(g) {}
	SimilarBiclique local_exact_query(int q, int size, double init_sim = -1.0);
protected:
	void _sim_rule(int q, double sim) {
		return;
	}
};

class LocalExactNoSFS : public LocalExact {
public:
	LocalExactNoSFS(const BiGraph& g) : LocalExact(g) {}
protected:
	void _enum(const SimilarBiclique& cur, const VI& P);
	void _simiarity_first_search_on_subgraph(int q);
};

class LocalExactNoSFS2 : public LocalExact {
public:
	LocalExactNoSFS2(const BiGraph& g) : LocalExact(g) {}
protected:
	void _simiarity_first_search_on_subgraph(int q);
};
