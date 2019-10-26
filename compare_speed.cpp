/*
  Compare speed between set and array. Set faster WHEN SIZE IS AROUND 100
*/

#include <iostream>
#include <queue>
#include <set>
#include <stdio.h> /* printf */
#include <time.h>  /* clock_t, clock, CLOCKS_PER_SEC */
#include <math.h>
#include <algorithm>

using namespace std;
int ITEM_NUM = 100;
int main()
{
  clock_t t_1;
  t_1 = clock();
  set<int> s;
  for (int i = 0; i < ITEM_NUM; i++)
  {
    s.insert(rand() % 10000);
  }
  for (int i = 0; i < ITEM_NUM; i++)
  {
    cout << (s.find(rand() % 10000) != s.end());
  }

  cout << "\n";
  t_1 = clock() - t_1;
  printf("It took me %d clicks (%f seconds).\n", t_1, ((float)t_1) / CLOCKS_PER_SEC);

  clock_t t_2;
  t_2 = clock();
  vector<int> v;
  for (int i = 0; i < ITEM_NUM; i++)
  {
    v.push_back(rand() % 10000);
  }
  for (int i = 0; i < ITEM_NUM; i++)
  {
    cout << (find(v.begin(), v.end(), rand() % 10000) != v.end());
  }
  cout << "\n";
  t_2 = clock() - t_2;
  printf("It took me %d clicks (%f seconds).\n", t_2, ((float)t_2) / CLOCKS_PER_SEC);
}