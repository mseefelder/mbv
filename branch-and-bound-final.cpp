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
#include "mbv-grasp.h"

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


int main(int argc, char const *argv[])
{
	int nVertices, mEdges, trash;

	cin >> nVertices >> mEdges >> trash;
	//cout << nVertices << " vertices and " << mEdges << " edges."<<endl;
	cout << nVertices << "\t" << mEdges << "\t";

	// Create and fill Graph
	Graph G(nVertices, mEdges);
	int a, b;
	for (int e = 0; e <= mEdges; e++) {
		cin >> a >> b >> trash;
		G.addEdge(a-1, b-1);
	}

	// Calculate heuristic starting point
	int heuristicResult = MBVGrasp (G, 50, 25);
	int minimumSol = heuristicResult;
	
	int repetitions = 5;
		
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
		float rr = 0.15;

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

	cout<< timeAccum/(float)(repetitions);

	cout << endl;
	
	return 0;
}