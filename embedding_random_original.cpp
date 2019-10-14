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

// random generator function:
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

// sort vertex order by size
class Vertex
{
public:
  int index;
  int size;
};
bool sortByVertexSize(const Vertex &v1, const Vertex &v2)
{
  return v1.size > v2.size;
}

vector<int> sortVertexOrder(Graph &H)
{
  vector<Vertex> toSort;
  vector<int> myvector;
  for (int i = 0; i < H.order(); i++)
  {
    Vertex v;
    v.index = i;
    v.size = H.deg(i);
    toSort.push_back(v);
  }
  sort(toSort.begin(), toSort.end(), sortByVertexSize);
  for (int i = 0; i < H.order(); i++)
    myvector.push_back(toSort.at(i).index);
  return myvector;
}

// initialize the mapping
vector<vector<int>> initializeMapping(int size)
{
  vector<vector<int>> mapping;
  for (int i = 0; i < size; i++)
  {
    vector<int> v;
    mapping.push_back(v);
  }
  return mapping;
}

// improving 1, minimize max number of overlap
int getMaxOverlap(vector<int> overlap)
{
  int max = 0;
  for (int i = 0; i < overlap.size(); i++)
  {
    max = overlap.at(i) > max ? overlap.at(i) : max;
  }
  return max;
}

// improving 2, minimize total number of physical qubits
int getTotalSize(vector<vector<int>> mapping)
{
  int total = 0;
  for (int i = 0; i < mapping.size(); i++)
  {
    total += mapping.at(i).size();
  }
  return total;
}

// define vertex weight - weight of a vertex in G grows exponentialy with num of vertices of H represented there
int getVertexWeight(int index, vector<int> overlap, bool isIn)
{
  int diameter = 30;
  int exponential = overlap.at(index);
  if (isIn)
    exponential--;
  return pow(diameter, exponential);
}

// check all neighbors empty
bool checkAllNeighbors(Graph &H, vector<vector<int>> mapping, int currentVertex)
{
  for (int i = 0; i < H.deg(currentVertex); i++)
  {
    if (mapping.at(H.adj[currentVertex].at(i)).size() != 0)
      return false;
  }
  return true;
}

// find the weight for shortest path from A to B
void dijkstraComputePaths(int source, Graph &adjacency_list,
                          vector<int> &min_distance, vector<int> &previous, vector<int> weights)
{
  weights.push_back(0);
  int n = adjacency_list.order();
  min_distance.clear();
  min_distance.resize(n, 1 << 30);
  min_distance[source] = weights.at(source);
  ;
  previous.clear();
  previous.resize(n, -1);
  set<pair<int, int>> vertex_queue;
  vertex_queue.insert(make_pair(min_distance[source], source));

  while (!vertex_queue.empty())
  {
    int dist = vertex_queue.begin()->first;
    int u = vertex_queue.begin()->second;
    vertex_queue.erase(vertex_queue.begin());

    // Visit each edge exiting u

    for (int i = 0; i < adjacency_list.deg(u); i++)
    {
      int v = adjacency_list.adj[u].at(i);
      int weight = weights.at(v);
      int distance_through_u = dist + weight;
      if (distance_through_u < min_distance[v])
      {
        vertex_queue.erase(make_pair(min_distance[v], v));
        min_distance[v] = distance_through_u;
        previous[v] = u;
        vertex_queue.insert(make_pair(min_distance[v], v));
      }
    }
  }
}

// find the shortest
vector<int> DijkstraGetShortestPathTo(
    int vertex, const std::vector<int> &previous)
{
  vector<int> path;
  for (; vertex != -1; vertex = previous[vertex])
    path.push_back(vertex);
  return path;
}

//create a dummy graph, add a dummy vertex which adjacent to every vertex in the model
void createDummyGraph(Graph &dummyGraph, vector<int> &tempWeight, Graph &G, Graph &H, vector<vector<int>> mapping, int currentVertex, int index)
{
  dummyGraph = G;
  dummyGraph.adj.resize(G.order() + 1);
  dummyGraph.adj[G.order()].resize(0);

  for (int k = 0; k < mapping.at(H.adj[currentVertex].at(index)).size(); k++)
  {
    dummyGraph.addArc(G.order(), mapping.at(H.adj[currentVertex].at(index)).at(k));
    dummyGraph.addArc(mapping.at(H.adj[currentVertex].at(index)).at(k), G.order());
    tempWeight[mapping.at(H.adj[currentVertex].at(index)).at(k)] = 0; //excluding w(neightbor vertex model)
  }
}

//find the union paths from root g* to each vertex model
void updateModels(Graph &G, Graph &H, int root, int currentVertex, vector<int> weight, vector<vector<int>> &mapping, vector<int> &overlap)
{
  map<int, int> vertex_usage;
  int dist[G.order()];

  vector<int> path;
  for (int i = 0; i < H.deg(currentVertex); i++)
  {
    path.clear();
    if (mapping.at(H.adj[currentVertex].at(i)).size() != 0)
    {
      Graph dummyGraph(0);
      vector<int> tempWeight = weight;
      createDummyGraph(dummyGraph, tempWeight, G, H, mapping, currentVertex, i);
      vector<int> min_distance;
      vector<int> previous;
      dijkstraComputePaths(root, dummyGraph, min_distance, previous, tempWeight);
      path = DijkstraGetShortestPathTo(G.order(), previous);
    }
    for (int j = 2; j < path.size(); j++)
    { // don't want include the dummy vertex and the vertex in neighbor's map
      if (vertex_usage.find(path.at(j)) == vertex_usage.end())
      {
        vertex_usage.insert(pair<int, int>(path.at(j), 1)); // if not exist, create a new one
        dist[path.at(j)] = H.adj[currentVertex].at(i);      // destination
      }
      else
      {
        vertex_usage[path.at(j)]++; // update the usage
      }
    }
  }

  // remove old usage
  for (int j = 0; j < mapping.at(currentVertex).size(); j++)
  {
    overlap.at(mapping.at(currentVertex).at(j))--;
  }
  mapping.at(currentVertex).clear();

  //update all vertex models
  mapping.at(currentVertex).push_back(root);
  overlap.at(root)++;
  for (map<int, int>::iterator ii = vertex_usage.begin(); ii != vertex_usage.end(); ++ii)
  {
    if ((*ii).first == root)
      continue;
    if ((*ii).second > 1)
    { // appears in more than one path, assign to new model ****************************
      if (find(mapping.at(currentVertex).begin(), mapping.at(currentVertex).end(), (*ii).first) == mapping.at(currentVertex).end())
      {
        mapping.at(currentVertex).push_back((*ii).first);
        overlap.at((*ii).first)++;
      }
    }
    else
    { // appears once, update
      if (find(mapping.at(dist[(*ii).first]).begin(), mapping.at(dist[(*ii).first]).end(), (*ii).first) == mapping.at(dist[(*ii).first]).end())
      {
        mapping.at(dist[(*ii).first]).push_back((*ii).first);
        overlap.at((*ii).first)++;
      }
    }
  }
}

//find the minimal vertex model
void findMinimalVertexModel(Graph &G, Graph &H, vector<int> weight, vector<vector<int>> &mapping, vector<int> &overlap, int currentVertex)
{
  // if all neighbors are empty
  if (checkAllNeighbors(H, mapping, currentVertex))
  {
    vector<int> vertexOrder = randomizeVertexOrder(G.order());
    for (int i = 0; i < vertexOrder.size(); i++)
    {
      if (weight.at(vertexOrder.at(i)) == 1 && G.adj[vertexOrder.at(i)].size() != 0)
      { // isolated vertex, ingore
        mapping.at(currentVertex).push_back(vertexOrder.at(i));
        overlap.at(vertexOrder.at(i))++; // update num of vertices of H represented at vertex of G
        return;
      }
    }
  }

  // new improvement:  choose Candidate from first level
  vector<int> rootCandidate;
  vector<int> rootCandidateWeight;
  for (int i = 0; i < H.deg(currentVertex); i++)
  {
    for (int j = 0; j < mapping.at(H.adj[currentVertex].at(i)).size(); j++)
    {
      for (int k = 0; k < G.adj[mapping.at(H.adj[currentVertex].at(i)).at(j)].size(); k++)
      {
        int neigh = G.adj[mapping.at(H.adj[currentVertex].at(i)).at(j)].at(k);
        if (find(rootCandidate.begin(), rootCandidate.end(), neigh) == rootCandidate.end() && weight.at(neigh) <= 1)
        { // if not contain
          rootCandidate.push_back(neigh);
        }
      }
    }
  }

  int cost[rootCandidate.size()][H.deg(currentVertex)];
  int minCost = 1 << 30;
  int root = 0;
  for (int i = 0; i < rootCandidate.size(); i++)
  {
    int sumCost = 0;
    int current = rootCandidate.at(i);
    for (int j = 0; j < H.deg(currentVertex); j++)
    {
      if (mapping.at(H.adj[currentVertex].at(j)).size() == 0)
      {
        cost[i][j] = 0;
      }
      else if (find(mapping.at(H.adj[currentVertex].at(j)).begin(), mapping.at(H.adj[currentVertex].at(j)).end(), current) != mapping.at(H.adj[currentVertex].at(j)).end())
      {
        cost[i][j] = weight.at(current);
      }
      else
      {
        // add a dummy vertex which adjacent to every vertex in the model - compute the distance from a vertex to a vertex-model
        Graph dummyGraph(0);
        vector<int> tempWeight = weight;
        createDummyGraph(dummyGraph, tempWeight, G, H, mapping, currentVertex, j);
        vector<int> min_distance;
        vector<int> previous;
        dijkstraComputePaths(current, dummyGraph, min_distance, previous, tempWeight);
        cost[i][j] = min_distance.at(G.order());
      }
      sumCost += cost[i][j];
      if (sumCost > minCost)
        break;
    }
    rootCandidateWeight.push_back(sumCost);
    // find root g* with minimum cost
    if (sumCost <= minCost)
    {
      minCost = sumCost;
    }
  }

  // random get the root from the first level rootCandidate
  vector<int> candidates;
  for (int i = 0; i < rootCandidate.size(); i++)
  {
    if (minCost == rootCandidateWeight.at(i))
      candidates.push_back(rootCandidate.at(i));
  }

  // random get the root from the candidates
  if (candidates.size() == 0)
  {
    vector<int> vertexOrder = randomizeVertexOrder(G.order());
    for (int i = 0; i < G.order(); i++)
    {
      if (G.deg(vertexOrder.at(i)) == 0)
      { // isolated vertex, ingore
        continue;
      }
      if (weight.at(vertexOrder.at(i)) <= 1)
      { // if no candidates, random choose one
        root = vertexOrder.at(i);
        break;
      }
    }
  }
  else
  {
    random_shuffle(candidates.begin(), candidates.end()); // using built-in random generator:
    root = candidates.at(0);
  }

  //find the union paths from root g* to each vertex model
  updateModels(G, H, root, currentVertex, weight, mapping, overlap);
}

// main method for finding Minor embedding
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
  cout << "Running Time : " << (double)(ends - start) / CLOCKS_PER_SEC << endl;
  return 0;
}
