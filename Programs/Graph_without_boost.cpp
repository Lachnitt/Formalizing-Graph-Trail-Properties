#include <iostream>
#include <vector>
#include <tuple>
#include <map>
#include <algorithm>
#include <iomanip>
#include <numeric>

typedef unsigned int vertex;
typedef unsigned int weight;

typedef std::vector<vertex> vertexList;
typedef std::map<std::pair<vertex,vertex>,weight> weightedGraph;
typedef std::vector<std::pair<vertex, vertex>> edgeList;
typedef std::vector<std::tuple<vertex, vertex,int>> markedEdgeList;

/// <summary>
/// Splits up an odd number.
/// </summary>
/// <param name="n"> Number of vertices in the graph, has to be odd   </param>
/// <returns> Numbers m and m' such that n = 2*m+2*m'+1 </returns>
std::pair<int, int> findDivisionSize(int n) 
{
	n = --n;
	n = n / 2;
	if (n % 2) {
		return std::make_pair(n / 2 + 1, n / 2);
	}
	else {
		return std::make_pair(n / 2, n / 2);
	}
}

void addToEdgeList(std::vector<markedEdgeList> &G, int k, int i, int j, int mark) 
{
	G[k].push_back(std::make_tuple(i,j,mark));
}

/// <summary>
/// Recursively splits up a graph until it has 7,9, or 11 vertices and then merges the resulting subgraphs until 
/// n/2 decompositions with n-1 edges are generated
/// </summary>
/// <param name="V"> Set of vertices </param>
/// <returns> (G,A,B,Alpha,Beta), G is a list of edges which are marked with 0,1, or 2.
/// A, B, Alpha and Beta are lists of vertices </returns>
std::tuple<std::vector<markedEdgeList>, vertexList, vertexList, vertexList, vertexList> findDecompOdd(vertexList &V) 
{
	int n = V.size();

	std::vector<markedEdgeList> G;
	vertexList A;
	vertexList B;
	vertexList Alpha;
	vertexList Beta;

	if (n == 7) {
		G = std::vector<markedEdgeList>(3, markedEdgeList());

		addToEdgeList(G, 0, V[0], V[6],1);
		addToEdgeList(G, 0, V[1], V[6],2);
		addToEdgeList(G, 0, V[1], V[2],0);
		addToEdgeList(G, 0, V[1], V[3],0);
		addToEdgeList(G, 0, V[1], V[4],1);
		addToEdgeList(G, 0, V[0], V[4],2);
		addToEdgeList(G, 0, V[4], V[5],0);

		addToEdgeList(G, 1, V[0], V[1],0);
		addToEdgeList(G, 1, V[0], V[2],1);
		addToEdgeList(G, 1, V[0], V[3],2);
		addToEdgeList(G, 1, V[2], V[6],2);
		addToEdgeList(G, 1, V[3], V[6],1);
		addToEdgeList(G, 1, V[3], V[4],0);
		addToEdgeList(G, 1, V[3], V[5],0);

		addToEdgeList(G, 2, V[0], V[5],0);
		addToEdgeList(G, 2, V[1], V[5],0);
		addToEdgeList(G, 2, V[2], V[5],1);
		addToEdgeList(G, 2, V[2], V[3],0);
		addToEdgeList(G, 2, V[2], V[4],2);
		addToEdgeList(G, 2, V[4], V[6],1);
		addToEdgeList(G, 2, V[5], V[6],2);

		A.push_back(V[0]);
		A.push_back(V[1]);
		A.push_back(V[2]);

		B.push_back(V[3]);
		B.push_back(V[4]);
		B.push_back(V[5]);

		Alpha.push_back(V[1]);
		Alpha.push_back(V[0]);
		Alpha.push_back(V[2]);

		Beta.push_back(V[4]);
		Beta.push_back(V[3]);
		Beta.push_back(V[5]);

		return std::make_tuple(G,A,B,Alpha,Beta);
	}
	else if (n == 9) {
		G = std::vector<markedEdgeList>(4, markedEdgeList());

		addToEdgeList(G, 0, V[0], V[8], 0);
		addToEdgeList(G, 0, V[1], V[8], 0);
		addToEdgeList(G, 0, V[2], V[8], 0);
		addToEdgeList(G, 0, V[3], V[8], 1);
		addToEdgeList(G, 0, V[3], V[4], 2);
		addToEdgeList(G, 0, V[4], V[5], 1);
		addToEdgeList(G, 0, V[5], V[6], 2);
		addToEdgeList(G, 0, V[6], V[7], 1);
		addToEdgeList(G, 0, V[7], V[8], 2);

		addToEdgeList(G, 1, V[0], V[1], 1);
		addToEdgeList(G, 1, V[1], V[2], 2);
		addToEdgeList(G, 1, V[2], V[3], 1);
		addToEdgeList(G, 1, V[3], V[5], 2);
		addToEdgeList(G, 1, V[5], V[7], 1);
		addToEdgeList(G, 1, V[7], V[0], 2);
		addToEdgeList(G, 1, V[0], V[6], 0);
		addToEdgeList(G, 1, V[5], V[8], 0);
		addToEdgeList(G, 1, V[2], V[4], 0);

		addToEdgeList(G, 2, V[2], V[5], 1);
		addToEdgeList(G, 2, V[5], V[1], 2);
		addToEdgeList(G, 2, V[1], V[6], 1);
		addToEdgeList(G, 2, V[6], V[2], 2);
		addToEdgeList(G, 2, V[1], V[7], 0);
		addToEdgeList(G, 2, V[0], V[2], 0);
		addToEdgeList(G, 2, V[6], V[8], 0);
		addToEdgeList(G, 2, V[6], V[4], 0);
		addToEdgeList(G, 2, V[1], V[3], 0);

		addToEdgeList(G, 3, V[1], V[4], 0);
		addToEdgeList(G, 3, V[2], V[7], 0);
		addToEdgeList(G, 3, V[8], V[4], 0);
		addToEdgeList(G, 3, V[6], V[3], 0);
		addToEdgeList(G, 3, V[0], V[5], 0);
		addToEdgeList(G, 3, V[0], V[3], 1);
		addToEdgeList(G, 3, V[3], V[7], 2);
		addToEdgeList(G, 3, V[7], V[4], 1);
		addToEdgeList(G, 3, V[4], V[0], 2);

		A.push_back(V[0]);
		A.push_back(V[1]);
		A.push_back(V[2]);
		A.push_back(V[3]);

		B.push_back(V[4]);
		B.push_back(V[5]);
		B.push_back(V[6]);
		B.push_back(V[7]);

		Alpha.push_back(V[3]);
		Alpha.push_back(V[2]);
		Alpha.push_back(V[1]);
		Alpha.push_back(V[0]);

		Beta.push_back(V[7]);
		Beta.push_back(V[5]);
		Beta.push_back(V[6]);
		Beta.push_back(V[4]);

		return std::make_tuple(G, A, B, Alpha, Beta);
	}
	else if (n == 11) {
		G = std::vector<markedEdgeList>(5, markedEdgeList());

		addToEdgeList(G, 0, V[1], V[10], 1);
		addToEdgeList(G, 0, V[10], V[2], 2);
		addToEdgeList(G, 0, V[2], V[7], 1);
		addToEdgeList(G, 0, V[7], V[1], 2);
		addToEdgeList(G, 0, V[0], V[7], 0);
		addToEdgeList(G, 0, V[9], V[7], 0);
		addToEdgeList(G, 0, V[8], V[7], 0);
		addToEdgeList(G, 0, V[2], V[3], 0);
		addToEdgeList(G, 0, V[2], V[4], 0);
		addToEdgeList(G, 0, V[2], V[5], 0);
		addToEdgeList(G, 0, V[2], V[6], 0);

		addToEdgeList(G, 1, V[0], V[9], 0);
		addToEdgeList(G, 1, V[1], V[9], 0);
		addToEdgeList(G, 1, V[2], V[9], 0);
		addToEdgeList(G, 1, V[4], V[5], 0);
		addToEdgeList(G, 1, V[4], V[6], 0);
		addToEdgeList(G, 1, V[4], V[7], 0);
		addToEdgeList(G, 1, V[4], V[8], 0);
		addToEdgeList(G, 1, V[4], V[9], 1);
		addToEdgeList(G, 1, V[9], V[3], 2);
		addToEdgeList(G, 1, V[3], V[10], 1);
		addToEdgeList(G, 1, V[10], V[4], 2);

		addToEdgeList(G, 2, V[1], V[2], 0);
		addToEdgeList(G, 2, V[1], V[3], 0);
		addToEdgeList(G, 2, V[1], V[4], 0);
		addToEdgeList(G, 2, V[6], V[7], 0);
		addToEdgeList(G, 2, V[6], V[8], 0);
		addToEdgeList(G, 2, V[6], V[9], 0);
		addToEdgeList(G, 2, V[6], V[0], 0);
		addToEdgeList(G, 2, V[6], V[1], 1);
		addToEdgeList(G, 2, V[1], V[5], 2);
		addToEdgeList(G, 2, V[5], V[10], 1);
		addToEdgeList(G, 2, V[10], V[6], 2);

		addToEdgeList(G, 3, V[8], V[9], 0);
		addToEdgeList(G, 3, V[8], V[0], 0);
		addToEdgeList(G, 3, V[8], V[1], 0);
		addToEdgeList(G, 3, V[8], V[2], 0);
		addToEdgeList(G, 3, V[3], V[4], 0);
		addToEdgeList(G, 3, V[3], V[5], 0);
		addToEdgeList(G, 3, V[3], V[6], 0);
		addToEdgeList(G, 3, V[3], V[7], 1);
		addToEdgeList(G, 3, V[7], V[10], 2);
		addToEdgeList(G, 3, V[10], V[8], 1);
		addToEdgeList(G, 3, V[8], V[3], 2);

		addToEdgeList(G, 4, V[0], V[1], 0);
		addToEdgeList(G, 4, V[0], V[2], 0);
		addToEdgeList(G, 4, V[0], V[3], 0);
		addToEdgeList(G, 4, V[0], V[4], 0);
		addToEdgeList(G, 4, V[5], V[6], 0);
		addToEdgeList(G, 4, V[5], V[7], 0);
		addToEdgeList(G, 4, V[5], V[8], 0);
		addToEdgeList(G, 4, V[5], V[9], 1);
		addToEdgeList(G, 4, V[9], V[10], 2);
		addToEdgeList(G, 4, V[10], V[0], 1);
		addToEdgeList(G, 4, V[0], V[5], 2);


		A.push_back(V[0]);
		A.push_back(V[1]);
		A.push_back(V[2]);
		A.push_back(V[3]);
		A.push_back(V[4]);

		B.push_back(V[5]);
		B.push_back(V[6]);
		B.push_back(V[7]);
		B.push_back(V[8]);
		B.push_back(V[9]);

		Alpha.push_back(V[2]);
		Alpha.push_back(V[4]);
		Alpha.push_back(V[1]);
		Alpha.push_back(V[3]);
		Alpha.push_back(V[0]);

		Beta.push_back(V[7]);
		Beta.push_back(V[9]);
		Beta.push_back(V[6]);
		Beta.push_back(V[8]);
		Beta.push_back(V[5]);

		return std::make_tuple(G, A, B, Alpha, Beta);
	}

	std::pair<int, int> p = findDivisionSize(n);
	int m = p.first;
	int k = p.second;

	vertexList V1(V.begin(), V.begin() + 2 * m);
	V1.push_back(V.back());
	vertexList V2(V.begin() + 2 * m, V.end() - 1);
	V2.push_back(V.back());
	
	std::tuple<std::vector<markedEdgeList>, vertexList, vertexList, vertexList, vertexList> res1 = findDecompOdd(V1);
	std::tuple<std::vector<markedEdgeList>, vertexList, vertexList, vertexList, vertexList> res2 = findDecompOdd(V2);

	std::vector<markedEdgeList> G1 = std::get<0>(res1);
	std::vector<markedEdgeList> G2 = std::get<0>(res2);

	vertexList A1 = std::get<1>(res1);
	vertexList A2 = std::get<1>(res2);

	vertexList B1 = std::get<2>(res1);
	vertexList B2 = std::get<2>(res2);

	vertexList Alpha1 = std::get<3>(res1);
	vertexList Alpha2 = std::get<3>(res2);

	vertexList Beta1 = std::get<4>(res1);
	vertexList Beta2 = std::get<4>(res2);


	//Generate new tuple:
	A.insert(A.end(), A1.begin(), A1.end());
	A.insert(A.end(), A2.begin(), A2.end());

	B.insert(B.end(), B1.begin(), B1.end()); //Ich brauche eigentlich keinen neuen
	B.insert(B.end(), B2.begin(), B2.end());

	Alpha.insert(Alpha.end(), Alpha1.begin(), Alpha1.end());
	Alpha.insert(Alpha.end(), Alpha2.begin(), Alpha2.end());

	Beta.insert(Beta.end(),Beta1.begin(), Beta1.end());
	Beta.insert(Beta.end(), Beta2.begin(), Beta2.end());

	G.insert(G.end(),G1.begin(), G1.end());
	G.insert(G.end(), G2.begin(), G2.end());

	for (int i = 0; i < m;++i) {
		for (auto v : A2) {
			G[i].push_back(std::make_tuple(Alpha1[i],v,0));
		}
		for (auto v : B2) {
			G[i].push_back(std::make_tuple(Beta1[i], v,0));
		}
	}

	for (int i = m; i < m+k; ++i) {
		for (auto v : B1) {
			G[i].push_back(std::make_tuple(Alpha2[i-m], v,0));
			//G[i].push_back(std::make_tuple(v, Beta2[i],0));
		}
		for (auto v : A1) {
			G[i].push_back(std::make_tuple(Beta2[i-m], v,0));
			//G[i].push_back(std::make_tuple(v, Alpha2[i],0));
		}
	}

	return std::make_tuple(G, A, B, Alpha, Beta);

}

/// <summary>
/// Searches for edges in a graph by their weight
/// </summary>
/// <param name="G"> Edge-weighted graph </param>
/// <param name="w"> Weight </param>
/// <param name="n"> Number of vertices in the graph </param>
/// <returns> Edge that is labelled with weight w </returns>
std::pair<vertex, vertex> getWeight(weightedGraph &G, const weight &w, const int &n) {
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			if (G.find(std::make_pair(i,j)) ->second == w) {
				return std::make_pair(i, j);
			}
		}
	}

	return std::make_pair(0, 0);
}

/// <summary>
/// Adds an edge to a edge-weighted graph
/// </summary>
/// <param name="G"> Edge-weighted graph </param>
/// <param name="i"> First endpoint of edge </param>
/// <param name="j"> Second endpoint of edge </param>
/// <param name="w"> Weight the edge should be labelled with </param>
void addToGraph(weightedGraph &G,  vertex i, vertex j, weight w) {
	G.insert({ std::make_pair(i,j),w });
	G.insert({ std::make_pair(j,i),w });
}

/// <summary>
/// Checks if a vertex is an endpoint of a edge in an edge list
/// </summary>
/// <param name="decomp"> Current decomposition as an edge list </param>
/// <param name="v"> Vertex that should be found in edge list </param>
/// <returns> True, if v is an endpoint of an edge in decomp. False, otw. </returns>
bool isEndpointInVector(edgeList decomp, vertex v) {
	for (int i = 0; i < decomp.size(); ++i) {
		if (decomp[i].first == v || decomp[i].second == v) {
			return true;
		}
	}
	return false;
}

/// <summary>
/// Adds n/2 edges a graph that are were not in the graph before.
/// </summary>
/// <param name="n"> Number of vertices in the graph </param>
/// <param name="decomp"> List of edges, empty in the beginning </param>
/// <param name="G"> Edge-weighted graph </param>
/// <param name="L"> Labelling function, empty in the beginning </param>
/// <param name="r"> Index of current round </param>
/// <returns></returns>
edgeList findDecompositionEven(int &n, std::vector<std::pair<vertex,vertex>> &decomp, weightedGraph &G, std::vector<weight> &L, int &r) {
	if (decomp.size() == n/2) {
		return decomp;
	}

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			if (i != j && G.find(std::make_pair(i, j)) == G.end() && L[i] == r && L[j] == r
				&& !isEndpointInVector(decomp,i) && !isEndpointInVector(decomp, j)) { 
				//isEndpointInVector(decomp,i) can be moved outside of the loop
				decomp.push_back(std::make_pair(i, j));
				edgeList res = findDecompositionEven(n,decomp,G,L,r);

				if (res.size() == n/2 || res.size() > decomp.size()) {
					return res;
				}
				decomp.pop_back();
			}
		}
	}
}

/// <summary>
/// Creates a graph with n vertices w/o a trail of length n.
/// </summary>
/// <param name="n"> Number of vertices in the graph, has to be odd </param>
/// <returns> Edge-weighted graph w/o trail of length n </returns>
weightedGraph createOddGraph(int &n) {
	vertexList V(n, 0);
	std::iota(V.begin(), V.end(), 0);
	weightedGraph G;
	auto res = findDecompOdd(V);
	int weight = 1;

	for (auto decomp : std::get<0>(res)) {

		for (auto e : decomp) {
			if (std::get<2>(e) == 1) {
				addToGraph(G, std::get<0>(e), std::get<1>(e), weight);
				++weight;
			}
		}
		for (auto e : decomp) {
			if (std::get<2>(e) == 0) {
				addToGraph(G, std::get<0>(e), std::get<1>(e), weight);
				++weight;
			}
		}
		for (auto e : decomp) {
			if (std::get<2>(e) == 2) {
				addToGraph(G, std::get<0>(e), std::get<1>(e), weight);
				++weight;
			}
		}
	}

	return G;
}

/// <summary>
/// Creates a graph with n vertices w/o a trail of length n.
/// </summary>
/// <param name="n"> Number of vertices in the graph, has to be even </param>
/// <returns> Edge-weighted graph w/o trail of length n </returns>
weightedGraph createEvenGraph(int &n) {
	std::vector<weight> L(n, 0);
	weightedGraph G;
	edgeList curr_decomp;
	int weight = 1;

	for (int r = 0; r < n-1; ++r) {
		curr_decomp = findDecompositionEven(n, curr_decomp, G, L, r);

		for (auto e : curr_decomp) {
			++L[e.first];
			++L[e.second];
			addToGraph(G,e.first,e.second,weight);
			++weight;
		}
		curr_decomp.clear();

	}
	return G;
}

/// <summary>
/// Creates a witness of a graph with n vertices and no trail of length n.
/// </summary>
/// <param name="n"> Number of vertices in the graph </param>
/// <returns> Edge-weighted graph with n vertices w/o trail of length n </returns>
weightedGraph createGraph(int &n) {
	if (n % 2) {
		return createOddGraph(n);
	}
	return createEvenGraph(n);
}

/// <summary>
/// Finds the longest ordered trail in an edge-weighted graph
/// </summary>
/// <param name="G"> Edge-weighted graph </param>
/// <param name="n"> Number of vertices in the graph </param>
/// <returns> Length of the longest ordered trail in the graph </returns>
int findLongestOrderedTrail(weightedGraph &G, const int &n) {

	std::vector<weight> L(n, 0);

	for (int i = 0; i < G.size() / 2; ++i) {
		std::pair<vertex, vertex> e = getWeight(G, i + 1, n);

		vertex temp = L[e.first];
		L[e.first] = std::max(L[e.first], L[e.second] + 1);
		L[e.second] = std::max(temp + 1, L[e.second]);
	}

	return *max_element(std::begin(L), std::end(L));
}

/// <summary>
/// Pretty prints a graph.
/// </summary>
/// <param name="G"> Edge-weighted graph that is to be printed </param>
/// <param name="n"> Number of vertices in the graph </param>
void printGraph(weightedGraph &G, const int &n) {
	int digits = 0;
	int m = n * (n - 1) / 2;

	while (m != 0) {
		m = m / 10;
		++digits;
	}

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			if (G.find(std::make_pair(i, j)) != G.end()) {
				std::cout << std::setw(digits) << (G.find(std::make_pair(i, j)))->second << " ";
			}
			else {
				std::cout << std::setw(digits) << "-" << " ";
			}
		}
		std::cout << std::endl;
	}
}

int main()
{ 
	int n1 = 4;
	weightedGraph G1;
	addToGraph(G1, 0, 1, 1);
	addToGraph(G1, 0, 2, 2);
	addToGraph(G1, 0, 3, 6);
	addToGraph(G1, 1, 2, 3);
	addToGraph(G1, 1, 3, 4);
	addToGraph(G1, 2, 3, 5);
	printGraph(G1,n1);
	std::cout << "The longest ordered Trail has length: " << findLongestOrderedTrail(G1, n1) << std::endl << std::endl;

	int n2 = 8;
	weightedGraph G2 = createGraph(n2);
	printGraph(G2, n2);
	std::cout << "The longest ordered Trail has length: " << findLongestOrderedTrail(G2, n2) << std::endl << std::endl;

	int n3 = 27;
	weightedGraph G3 = createGraph(n3);
	printGraph(G3, n3);
	std::cout << "The longest ordered Trail has length: " << findLongestOrderedTrail(G3, n3) << std::endl;

	return 0;
}

