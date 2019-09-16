/*
  Heuristic search Minor embedding algorithm <zyan632@aucklanduni.ac.nz> 
  Basic version with new improvements (random selection)
  usage: embedding.exe guest 
  put host graph (host.in) in the same folder
*/

#include "Graph.h"
#include <string>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <ctime>
#include <fstream>
#include <bits/stdc++.h>

vector<int> randomizeVertexOrder(int size)
{
  vector<int> myvector;
  for (int i = 0; i < size; ++i)
  {
    myvector.push_back(i); //set value 0 1 2 3 ... size-1
  }
  random_shuffle(myvector.begin(), myvector.end()); // using built-in random generator:
  return myvector;
}

class Qubit
{
public:
  int index;
  Graph graph;
  vector<Qubit> freeNeighbors;

  Qubit(int index, Graph graph, vector<Qubit> neighbors) : index(index), graph(graph), freeNeighbors(neighbors) {}
  ~Qubit() {}
};

// sort vertex order by size
class Vertex
{
public:
  int index;
  vector<Qubit> qubits;
  int size();

  Vertex();
  Vertex(Qubit qubit);
  ~Vertex();
  void addQubit(Qubit qubit);
};

int Vertex::size()
{
  return qubits.size();
}

void Vertex::addQubit(Qubit qubit)
{
  qubits.push_back(qubit);
}

bool findMinorEmbedding(Graph &G, Graph &H)
{
  srand(unsigned(time(0))); // for randomizing
#if 1
  vector<int> vertexOrder = randomizeVertexOrder(H.order());
#else
  vector<int> vertexOrder = sortVertexOrder(H); // sort order by size
#endif
  vector<int> weight;

  vector<int> overlap(G.order(), 0);
  vector<int> newOverlap(G.order(), 0);

  vector<vector<int>> mapping = initializeMapping(H.order());
  vector<vector<int>> newMapping = initializeMapping(H.order());

  int stage = 1;
  while (getMaxOverlap(newOverlap) < getMaxOverlap(overlap) || getTotalSize(newMapping) < getTotalSize(mapping) || stage <= 2)
  {
    mapping = newMapping;
    overlap = newOverlap;
    for (int i = 0; i < vertexOrder.size(); i++)
    {
      int currentVertex = vertexOrder.at(i);
      // assign weight
      weight.clear();
      for (int j = 0; j < G.order(); j++)
      {
        bool isIn = find(newMapping.at(currentVertex).begin(), newMapping.at(currentVertex).end(), j) != newMapping.at(currentVertex).end();
        weight.push_back(getVertexWeight(j, newOverlap, isIn));
      }
      findMinimalVertexModel(G, H, weight, newMapping, newOverlap, currentVertex);
    }
    stage++;
  }

  //check the usage for each vertex in G, if all <=1, then we find the embedding. otherwise "failure"
  bool embedding = true;
  for (int i = 0; i < overlap.size(); i++)
  {
    if (overlap.at(i) > 1)
    {
      embedding = false;
      break;
    }
  }

  //print the embedding, if find
  if (embedding)
  {
    int totalQubit = 0;
    int max = 0;
    cout << "[";
    for (int i = 0; i < mapping.size(); i++)
    {
      cout << "[";
      totalQubit = totalQubit + mapping.at(i).size();
      if (mapping.at(i).size() > max)
        max = mapping.at(i).size();
      sort(mapping.at(i).begin(), mapping.at(i).end());
      for (int j = 0; j < mapping.at(i).size() - 1; j++)
      {
        cout << mapping.at(i).at(j) << ", ";
      }
      cout << mapping.at(i).at(mapping.at(i).size() - 1);
      if (i == mapping.size() - 1)
      {
        cout << "]";
      }
      else
      {
        cout << "], " << endl;
      }
    }
    cout << "]" << endl;
    cout << "Total number of Qubits used is: " << totalQubit << endl;
    cout << "Max chain length is: " << max << endl;
    return true;
  }
  else
  {
    return false;
  }
}

int main(int argc, char **argv)
{
  clock_t start = clock();
  ifstream host("host.in");
  Graph G(0);
  host >> G;
  Graph H(0);
  cin >> H;
  float ratio = std::stof(argv[1]);
  bool foundEmbedding;
  cout << "Guest graph: " << endl;
  cout << H << endl;
  for (int i = 1; i <= 2; i++)
  {
    cout << "try " << i << ":" << endl;
    foundEmbedding = findMinorEmbedding(G, H);
    if (foundEmbedding)
    {
      break;
    }
  }
  if (!foundEmbedding)
    cout << "No embedding found" << endl;

  clock_t ends = clock();
  //cout << "Running Time : " << (double)(ends - start) / CLOCKS_PER_SEC << endl;
  return 0;
}