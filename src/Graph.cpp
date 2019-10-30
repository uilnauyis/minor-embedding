/*
 *  Simple Graph Data Structure (I/O)  -- <mjd@cs.auckland.ac.nz>
 */
#include "Graph.h"
#include <string>
#include <sstream>
#include <algorithm>

ostream &operator<<(ostream &o, const Graph &G)
{
  o << G.order() << "\n";
  for (int i = 0; i < G.order(); i++)
  {
    for (int j = 0; j < G.adj[i].size(); j++)
      o << G.adj[i][j] << " ";
    o << "\n";
  }
  return o;
}

istream &operator>>(istream &in, Graph &G)
{
  string line;
  int n;
  getline(in, line);
  istringstream lineIn(line);
  lineIn >> n;

  G.adj.resize(n);
  for (int i = 0; i < n; i++)
  {
    G.adj[i].resize(0);
  }

  for (int v1 = 0; v1 < n; v1++)
  {
    getline(in, line);
    istringstream lineIn(line);
    int v2;
    while (lineIn >> v2)
    {
      G.addArc(v1, v2);
    }
  }
  return in;
}
