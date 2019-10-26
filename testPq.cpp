/*
  This test shows that when size is about 500, 'set' performs better than 'priority queue'.
*/

#include <functional>
#include <queue>
#include <vector>
#include <iostream>
#include <utility>
#include <set>
#include <time.h>

using namespace std;

void print_queue(priority_queue<pair<int, int>> q)
{
  while (!q.empty())
  {
    std::cout << q.top().first << "," << q.top().second << " ";
    q.pop();
  }
  std::cout << '\n';
}

void print_set(set<pair<int, int>> s)
{
  while (!s.empty())
  {
    set<pair<int, int>>::iterator dist = s.begin();
    std::cout << (*dist).first << "," << (*dist).second << " ";
    s.erase(dist);
  }
  std::cout << '\n';
}


int main()
{
  clock_t t_1;
  t_1 = clock();

  priority_queue<pair<int, int>> q;

  for (int i = 0; i < 50000; i++)
  {
    q.push(make_pair(rand() % 1000, i));
  }

  print_queue(q);
  t_1 = clock() - t_1;
  printf("It took me %d clicks (%f seconds).\n", t_1, ((float)t_1) / CLOCKS_PER_SEC);

  t_1 = clock();
  set<pair<int, int>> s;

  for (int i = 0; i < 50000; i++)
  {
    s.insert(make_pair(rand() % 1000, i));
  }

  print_set(s);
  t_1 = clock() - t_1;
  printf("It took me %d clicks (%f seconds).\n", t_1, ((float)t_1) / CLOCKS_PER_SEC);
}