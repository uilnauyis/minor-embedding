/*
 *  Simple (Directed) Graph Data Structure
 *  <mjd@cs.auckland.ac.nz> -- Sept 2002
 */
#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <algorithm>

using namespace std;

typedef vector<int> NLIST;     // node (neighbor) list
typedef vector<NLIST> ADJLIST; // graph data structure

class Graph
{
private:
public:        // moved here for lucas
  ADJLIST adj; // adjacency lists.

  // creators + destroyers
  //
  Graph(int n) : adj(n) {}
  ~Graph() {}

  // stream I/0
  //
  friend ostream &operator<<(ostream &, const Graph &);
  friend istream &operator>>(istream &, Graph &);

  // mutators
  //
  bool addArc(int u, int v) // false if already exists
  {
    if (isArc(u, v))
      return false;

    adj[u].push_back(v);
    return true;
  }

  bool addEdge(int u, int v) // bi-directional arcs
  {
    return addArc(u, v) && addArc(v, u);
  }

  bool rmArc(int u, int v) // false if doesn't exist
  {
    int n = order();
    assert(u >= 0 && u < n);
    assert(v >= 0 && v < n);
    assert(u != v);

    const NLIST::iterator it = find(adj[u].begin(), adj[u].end(), v);
    if (it == adj[u].end())
      return false;
    adj[u].erase(it);
    return true;
  }

  bool rmEdge(int u, int v)
  {
    return rmArc(u, v) && rmArc(v, u);
  }

  // accessors
  //
  int order() const { return adj.size(); }

  bool isArc(int u, int v) const
  {
    int n = order();
    assert(u >= 0 && u < n);
    assert(v >= 0 && v < n);
    assert(u != v);

    NLIST::const_iterator it = find(adj[u].begin(), adj[u].end(), v);
    if (it == adj[u].end())
      return false;
    else
      return true;
  }

  bool isEdge(int u, int v) const
  {
    return isArc(u, v) && isArc(v, u);
  }

  // out-degree
  int deg(int u) const { return adj[u].size(); }

  const NLIST &N(int u) const { return (const NLIST &)adj[u]; }
};
#endif
