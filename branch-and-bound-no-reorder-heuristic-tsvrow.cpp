#include <limits>
#include <iostream>
#include <forward_list>
#include <vector>
#include <queue>
#include <stack>
#include <unordered_set>
#include <utility>
#include <algorithm>
#include <chrono>
#include <random>
#include "mbv-utils.h"

using namespace std;

int Backtrack(Graph &G, MBVSolutionUndo &sol, int &minimumSol, auto &start, int &depth, float &resortRatio)
{
	depth++;
	bool resort = (depth <= resortRatio*G.mEdges);
	int e = 0;
	int pruneBin = (depth*10.0)/(G.mEdges-1);
	//if(pruneBin > 9) cout<<"EITA "<<depth<<" "<<G.nVertices-1<<endl;
	// Prune by inviability
	if (sol.getBranchVertexCount() >= minimumSol || minimumSol == 0) {
		
	}
	// A valid solution was found
	else if (sol.getActiveEdgeCount() == G.nVertices-1) {
		if(sol.getBranchVertexCount() < minimumSol) {
			//cout<<"Solution "<<sol.getBranchVertexCount()<<"(";
			//auto end = chrono::steady_clock::now();
			//cout << chrono::duration <double, milli> (end-start).count() << " ms)" << endl;
			minimumSol = sol.getBranchVertexCount();
		}
	}
	// No solution yet
	else{
		for (int i = 0; i<G.mEdges; i++) {
			e = sol.edgeIndex[i];

			if(sol.getEdgeState(e) == 0) {
				if(sol.activateEdge(e)) {
					minimumSol = Backtrack(G, sol, minimumSol, start, depth, resortRatio);
					sol.undoActivateEdge();
				}

				sol.prohibitEdge(e);
				if(resort) sol.sortEdges();

				//Try kruskal ---
				int ke = 0;
				int opcount = 0;
				for (int k = 0; k < G.mEdges; ++k) {
					ke = sol.edgeIndex[k];
					if(sol.activateEdge(ke)) 
					{
						opcount++;
					}
					if(sol.getActiveEdgeCount() == G.nVertices - 1 ){
						break;
					}
				}

				if(sol.getBranchVertexCount() < minimumSol && sol.getActiveEdgeCount() == G.nVertices-1) {
					//cout<<"Solution on Kruskal "<<sol.getBranchVertexCount()<<"(";
					//auto end = chrono::steady_clock::now();
					//cout << chrono::duration <double, milli> (end-start).count() << " ms)" << endl;
					minimumSol = sol.getBranchVertexCount();
				}

				//Check if we should branch. Only branch if relaxed bound is < minimumSol and graph is connex
				bool doBranch = (sol.getRelaxedSolution() < minimumSol) && (sol.getActiveEdgeCount() == G.nVertices-1);

				for (int x = 0; x < opcount; ++x) {
					sol.undoActivateEdge();
				}

				//
				if(doBranch)
						minimumSol = Backtrack(G, sol, minimumSol, start, depth, resortRatio);

				sol.undoProhibitEdge(e);
				if(resort) sol.sortEdges();
				break;
			}

		}
	}
	depth--;
	return minimumSol;
}

void DFS (Graph &G, unordered_set<int>* neig, unordered_set<int> &result, int source, int destination) {

	//Auxiliary queue used during the algorithm
	std::stack<int> lifo;
	//For each vertex: store its parent
	int* parentEdge;
	parentEdge = new int[G.nVertices];
	//For each vertex: store its level
	int* levels;
	levels = new int[G.nVertices];
	//Fill arrays
	for (int i  = 0; i<G.nVertices; i++){
		parentEdge[i]=-2;
		levels[i]=-1;
	}
	//Set initial values
	parentEdge[source]=-1;
	levels[source]=0;
	lifo.push(source);

	int current = 0;
	//int iterations = 0;
	bool deleteFlag;
	int n = 0;
	while(!lifo.empty()){
		current = lifo.top();
		lifo.pop();
		for (int edge : neig[current]) {
			n = (G.edges[edge].first == current) ? G.edges[edge].second : G.edges[edge].first; 
			if (n == destination) {
				parentEdge[destination]=edge;
				levels[destination]=levels[current]+1;
				break;
			}
			else if (parentEdge[n] == -2) {
				lifo.push(n);
				parentEdge[n]=edge;
				levels[n]=levels[current]+1;
			}
		}
	}

	current = destination;
	int iterations = levels[destination];
	for (int i = 0; i < iterations; i++) {
		result.insert(parentEdge[current]);

		current = (G.edges[parentEdge[current]].first == current) ? G.edges[parentEdge[current]].second : G.edges[parentEdge[current]].first ;
	}

	delete [] parentEdge; delete [] levels;

	return;
}

int MBVGrasp (Graph &G, int maxIter, int randomMargin) {
	// Solution
	int minBranchVertexCount = numeric_limits<int>::max();

	// Random setup
	std::random_device r;
	std::default_random_engine engine(r());
	std::uniform_int_distribution<int> uDist(0, 0);

	// Create solution
	MBVSolutionUndo gS(G);

	// Make edges connecting vertices with degree 1 obligatory
	// Number of obligatory edges
	int obl = 0;
	//for (int i = 0; i < G.nVertices; ++i)
	//{
	//	if (G.vertexDegrees[i] == 1)
	//	{
	//		gS.activateEdge(G.vertices[i].front());
	//		obl++;
	//	}
	//}

	// Run maxIter iterations
	for (int iter = 0; iter < maxIter; ++iter)
	{
		// Find randomized greedy solution gS
		// Fill vector with edges
		vector<int> edges;
		edges.reserve(G.mEdges);
		for (int e = 0; e < G.mEdges; ++e) {
			edges.push_back(gS.edgeIndex[e]);
		}

		// Create vector for cycle-making edges
		queue<int> cycleEdges;

		// Auxiliar structure for local search (represents the tree)
		int treeDegrees[G.nVertices];
		unordered_set<int> neig[G.nVertices]; /// Fill this
		for (int i = 0; i < G.nVertices; ++i)
		{
			treeDegrees[i] = gS.vertexDegrees[i];
		}

		// Choose semi-randomly (randomMargin) with Kruskal
		int e = 0, eIndex = 0, u = 0, v = 0;
		for (int i = 0; i < G.mEdges; ++i) {
			// get random index (from end of vector to -randomMargin positions)
			eIndex = uDist(engine);
			// if there are less than randomMargin elements, fix index
			if (eIndex >= edges.size()) {
				eIndex = 0;
			}
			// get value of index on edges vector
			e = *(edges.begin()+edges.size()-1-eIndex);
			// try to add edge to solution
			if(!gS.activateEdge(e)) {
				cycleEdges.push(e);
			}
			else {
				// Get vertices connected by edge and their degrees onthe tree so far
				u = G.edges[e].first;
				v = G.edges[e].second;
				treeDegrees[u]++;
				treeDegrees[v]++;
				// Update neig array
				neig[u].insert(e);
				neig[v].insert(e);
			}

			// Remove edge from vector of edges left
			edges.erase(edges.begin()+edges.size()-1-eIndex); //something's wrong here

			// If kruskal is over
			if (gS.getActiveEdgeCount() == G.nVertices -1) {
				break;
			}
		}

		// Now, we have a solution on gS. 
		// On edges vector we have the edges that weren't used on Kruskal
		// On cycleEdges we have the ones that would form a cycle

		// Current branch vertex count
		int branchVertexCount = gS.getBranchVertexCount(); /// ALWAYS UPDATE THIS WHEN UPDATING TREE

		/**/
		// Add edges left from edges vector to cycleEdges vector
		while(!edges.empty()) {
				cycleEdges.push(edges.back());
				edges.pop_back();
		}

		// Number of iterations that will be executed to look for improvement
		int remaining = cycleEdges.size();
		// Number of edges reinserted on 
		int inserted = 0;

		// Local search starting on gS
		// Iterate over edges to evaluate local search neighbourhood
		e = 0, u = 0, v = 0;
		int uLine = 0, vLine = 0;
		while(remaining!=0) {
			// Remove last element from vector
			e = cycleEdges.front();
			cycleEdges.pop();

			// Get edge's vertices
			u = G.edges[e].first;
			v = G.edges[e].second;

			bool keepChange = false;
			// If edge is valid candidate
			if(treeDegrees[u] != 2 && treeDegrees[v] != 2){
				/// Virtually add e to tree, forming a cycle C (just in degrees for now)
				treeDegrees[u]++;
				treeDegrees[v]++;
				
				/// Get sequence SEQ of edges in cycle C (some kind of 'vector<int> DFS(vector<int>[] neig, u, v, size)' function)
				unordered_set<int> seq;
				DFS(G, neig, seq, u, v);

				//cout<<":: "<<u<<" --- "<<v<<endl;
				//for (int dummy : seq) cout<<G.edges[dummy].first<<","<<G.edges[dummy].second<<" -> "; cout<<endl;

				/// For each edge e'(i,j) on SEQ:
				for (int eLine : seq){
					uLine = G.edges[eLine].first;
					vLine = G.edges[eLine].second;
					/// If removing e' from tree reduces the degree of (i || j) from 3 to 2
					if (treeDegrees[uLine] == 3 || treeDegrees[vLine] == 3) {
						/// Remove e' from tree
						neig[uLine].erase(eLine);
						neig[vLine].erase(eLine);
						if(treeDegrees[uLine]-- == 3) branchVertexCount--;
						if(treeDegrees[vLine]-- == 3) branchVertexCount--;
						/// Add e' to cycleEdges
						cycleEdges.push(eLine);
						/// (remaining = inserted+1)
						remaining = inserted+1;
						/// (inserted = 0)
						inserted = 0;
						/// Add e to tree, updating parent and rank arrays
						neig[uLine].insert(eLine);
						neig[vLine].insert(eLine);
						/// Make sure e will stay in tree (keepChange = true)
						keepChange = true;
						/// Break
						break;
					}
				}

				/// If no improvement was encountered (keepChange == false):
				if (!keepChange) {
					/// Undo adding e to tree
					treeDegrees[u]--;
					treeDegrees[v]--; 
					/// Add e to cycleEdges
					cycleEdges.push(e);
					/// (inserted++)
					inserted++;
				}
			}
			/// (remaining--)
			remaining--;
			
		}
		/**/

		/// Save results
		minBranchVertexCount = (branchVertexCount < minBranchVertexCount) ? branchVertexCount : minBranchVertexCount;
		/// reset solution
		for(int i = 0; i < G.nVertices - 1 - obl; ++i) {	
			/// Undo nVertices - 1 - obl activations
			gS.undoActivateEdge();
		}

		uDist.param(uniform_int_distribution<int>::param_type(0, (iter>randomMargin) ? randomMargin : iter));
	}

	return minBranchVertexCount;
}

int main(int argc, char const *argv[])
{
	int nVertices, mEdges, trash;

	cin >> nVertices >> mEdges >> trash;
	//cout << nVertices << " vertices and " << mEdges << " edges."<<endl;
	cout << nVertices << "\t";

	// Create and fill Graph
	Graph G(nVertices, mEdges);
	int a, b;
	for (int e = 0; e <= mEdges; e++) {
		cin >> a >> b >> trash;
		G.addEdge(a-1, b-1);
	}

	// Calculate heuristic starting point
	int heuristicResult = MBVGrasp (G, (int)(mEdges*0.8), (int)(mEdges*0.75));
	//cout<<"Best heuristic result: "<<heuristicResult<<endl;
	int minimumSol = heuristicResult;
	
	int repetitions = 3;

	float rrs[10] = {0.01, 0.025, 0.05, 0.075, 0.1, 0.15, 0.20, 0.25, 0.35, 0.50};

	for (float rr : rrs) {
		
		auto timeAccum = chrono::duration <double, milli>(0.0).count();
		bool first = true;

		for (int repetition = 0; repetition < repetitions; repetition++) {
			// Make solution that will be used on branch and bound
			MBVSolutionUndo S(G);

			//make bridges obligatory
			for (int i = 0; i < nVertices; ++i)
			{
				if (G.vertexDegrees[i] == 1)
				{
					S.activateEdge(G.vertices[i].front());
				}
			}

			int depth = -1;
			float rr = rrs[repetition];

			int solution = minimumSol;

			auto start = chrono::steady_clock::now();
			solution = Backtrack(G, S, solution, start, depth, rr);
			auto end = chrono::steady_clock::now();
			auto dur = chrono::duration <double, milli> (end-start).count();

			if(first){
				timeAccum = dur;
				first = false;
			}
			else {
				timeAccum += dur;
			}
		}

		cout<< timeAccum/(float)(repetitions) <<"\t";
	}

	cout << endl;
	
	return 0;
}