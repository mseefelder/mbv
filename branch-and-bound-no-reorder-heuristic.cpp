#include <limits>
#include <iostream>
#include <forward_list>
#include <vector>
#include <queue>
#include <utility>
#include <algorithm>
#include <chrono>
#include <random>
#include "mbv-utils.h"

using namespace std;

int Backtrack(Graph &G, MBVSolutionUndo &sol, int &minimumSol, auto &start, long long &pruneCount)
{
	int e = 0;
	// Prune by inviability
	if (sol.getBranchVertexCount() >= minimumSol || minimumSol == 0) {
		return minimumSol;
	}
	// A valid solution was found
	else if (sol.getActiveEdgeCount() == G.nVertices-1) {
		if(sol.getBranchVertexCount() < minimumSol) {
			cout<<"Solution "<<sol.getBranchVertexCount()<<"(";
			auto end = chrono::steady_clock::now();
			cout << chrono::duration <double, milli> (end-start).count() << " ms)" << endl;
			minimumSol = sol.getBranchVertexCount();
		}
	}
	// No solution yet
	else{
		for (int i = 0; i<G.mEdges; i++) {
			e = sol.edgeIndex[i];

			if(sol.getEdgeState(e) == 0) {
				if(sol.activateEdge(e)) {
                    minimumSol = Backtrack(G, sol, minimumSol, start, pruneCount);
					sol.undoActivateEdge();
				}

				sol.prohibitEdge(e);				
				//sol.sortEdges();

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
					cout<<"Solution on Kruskal "<<sol.getBranchVertexCount()<<"(";
					auto end = chrono::steady_clock::now();
					cout << chrono::duration <double, milli> (end-start).count() << " ms)" << endl;
					minimumSol = sol.getBranchVertexCount();
				}

                //Check if we should branch
                bool doBranch = (sol.getRelaxedSolution() < minimumSol) && (sol.getActiveEdgeCount() == G.nVertices-1);
                if(!doBranch)
                    pruneCount++;

				for (int x = 0; x < opcount; ++x) {
					sol.undoActivateEdge();
				}

				//Only branch if relaxed bound is < minimumSol
                if(doBranch)
                    minimumSol = Backtrack(G, sol, minimumSol, start, pruneCount);

				sol.undoProhibitEdge(e);
				//sol.sortEdges();

				break;
			}

		}
	}

	return minimumSol;
}

int MBVGrasp (Graph &G, int maxIter, int randomMargin) {
	// Solution
	int minBranchVertexCount = numeric_limits<int>::max();

	// Random setup
	std::random_device r;
	std::default_random_engine engine(r());
	std::uniform_int_distribution<int> uDist(0, randomMargin);

	// Create solution
	MBVSolutionUndo gS(G);

	// Make edges connecting vertices with degree 1 obligatory
	// Number of obligatory edges
	int obl = 0;
	for (int i = 0; i < G.nVertices; ++i)
	{
		if (G.vertexDegrees[i] == 1)
		{
			gS.activateEdge(G.vertices[i].front());
			obl++;
		}
	}

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
		int parent[G.nVertices];
		int treeDegrees[G.nVertices];
		int rank[G.nVertices];
		for (int i = 0; i < G.nVertices; ++i)
		{
			parent[i] = -1;
			treeDegrees[i] = gS.vertexDegrees[i];
			rank[i] = 0;
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
				// Add edge to auxiliar structure for local search (in the form of a parent relation)
				if(treeDegrees[v] > treeDegrees[u]) {
					parent[u] = v;
					rank[u] = rank[v]+1;
				}
				else {
					parent[v] = u;
					rank[v] = rank[u]+1;
				}
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
		cout<<branchVertexCount<<endl;

		/**
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
				
				/// Get sequence SEQ of edges in cycle C
				
				/// For each edge e'(i,j) on SEQ:
					/// If removing e' from tree reduces the degree of (i || j) from 3 to 2
						/// Remove e' from tree
						/// Add e' to cycleEdges
						/// (remaining = inserted+1)
						/// (inserted = 0)
						/// Add e to tree, updating parent and rank arrays
						/// Make sure e will stay in tree (keepChange = true)
						/// Break

				/// If no improvement was encountered (keepChange == false):
					/// Undo adding e to tree 
					/// Add e to cycleEdges
					/// (inserted++)
			}
			/// (remaining--)
			
		}
		/**/

		/// Save results
		minBranchVertexCount = (branchVertexCount < minBranchVertexCount) ? branchVertexCount : minBranchVertexCount;
		/// reset solution
		for(int i = 0; i < G.nVertices - 1 - obl; ++i) {	
			/// Undo nVertices - 1 - obl activations
			gS.undoActivateEdge();
		}
	}

	return minBranchVertexCount;
}

int randomGreedy() {

}

int main(int argc, char const *argv[])
{
	int nVertices, mEdges, trash;

	cin >> nVertices >> mEdges >> trash;
	cout << nVertices << " vertices and " << mEdges << " edges."<<endl;

	Graph G(nVertices, mEdges);
	
	int a, b;
	for (int e = 0; e <= mEdges; e++) {
		cin >> a >> b >> trash;
		G.addEdge(a-1, b-1);
	}

	int heuristicResult = MBVGrasp (G, 30, 15);
	cout<<"Best heuristic result: "<<heuristicResult<<endl;
	
	/**
	int minimumSol = numeric_limits<int>::max();
	MBVSolutionUndo S(G);
	MBVSolutionUndo initialKruskal(G);

	//make bridges obligatory
	for (int i = 0; i < nVertices; ++i)
	{
		if (G.vertexDegrees[i] == 1)
		{
			S.activateEdge(G.vertices[i].front());
		}
	}

	//Run Kruskal to get a headstart
	int e = 0;
	for (int i = 0; i < mEdges; ++i)
	{
		e = initialKruskal.edgeIndex[i];
		bool activated = initialKruskal.activateEdge(e);

        cout<<e<<": "<<initialKruskal.edgeWeight[e]<<"(";
        cout<<G.vertexDegrees[G.edges[e].first]<<","<<G.vertexDegrees[G.edges[e].second];
        cout<<") "<<activated<<endl;

		if(initialKruskal.getActiveEdgeCount() == nVertices - 1 ){
            cout<<endl;
            cout<<"Relaxed solution "<<initialKruskal.getRelaxedSolution()<<endl;
			break;
		}
	}

	minimumSol = initialKruskal.getBranchVertexCount();
	cout<<"Initial heuristic: "<<minimumSol<<endl;
	/*
    for (int i = 0; i < nVertices-1; i++) {
        initialKruskal.undoActivateEdge();
    }
    cout<<"Relaxed solution "<<initialKruskal.getRelaxedSolution()<<endl;

    long long pruneCount = 0;

	auto start = chrono::steady_clock::now();
	minimumSol = Backtrack(G, S, minimumSol, start, pruneCount);
	auto end = chrono::steady_clock::now();
	cout << "Solution has " << minimumSol << " branch vertices! (";
	cout << chrono::duration <double, milli> (end-start).count() << " ms), ";
    cout << pruneCount << " prunes." << endl;
	*/
	return 0;
}