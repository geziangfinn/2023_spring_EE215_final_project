// core_router.cpp : Defines the entry point for the console application.
//
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <vector>
#include <unordered_map>

#include "omp.h"

#include "assert.h"
using namespace std;
///////////////////////////////////////////////////////////////////////////////////////////

struct COORDINATE
{
  size_t layer;
  size_t x;
  size_t y;

  bool operator<(const COORDINATE &coo) const
  {
    if (layer != coo.layer)
      return layer < coo.layer;

    if (x != coo.x)
      return x < coo.x;

    return y < coo.y;
  }

  bool operator==(const COORDINATE &coo) const
  {
    return layer == coo.layer && x == coo.x && y == coo.y;
  }
};

std::ostream &operator<<(std::ostream &o, const COORDINATE &c)
{
  o << '[' << c.layer << "," << c.x << "," << c.y << "]";
  return o;
}

/////////////////////////////////////////////////////////////////

struct NETCOORDINATES
{
  size_t netid;
  COORDINATE coo1, coo2;
};

typedef std::vector<NETCOORDINATES> NETLIST;

/////////////////////////////////////////////////////////////////

enum DIRECTION
{
  UNDEFINED = -1,
  NONE,
  N,
  S,
  E,
  W,
  U,
  D
};
vector<DIRECTION> allowed_expansion = {N, S, E, W, U, D};

// set<DIRECTION> allowed_expansion(v.begin(), v.end());

bool is_bend(DIRECTION d1, DIRECTION d2)
{
  if ((d1 == N || d1 == S) && (d2 == E || d2 == W))
    return true;

  if ((d1 == E || d1 == W) && (d2 == N || d2 == S))
    return true;

  return false;
}

bool is_via(DIRECTION d) { return (d == U || d == D); }

///////////////////////////////////////////////////////

int const COST_UNDEFINED = -2;
int const COST_BLOCKED = -1;

struct WAVE_CELL
{
  COORDINATE coord;
  int cost;
  DIRECTION pred;

  WAVE_CELL(void)
  {
    coord.layer = -1;
    coord.x = -1;
    coord.y = -1;
    cost = COST_UNDEFINED;
    pred = UNDEFINED;
  }
  bool operator>(const WAVE_CELL &cell) const { return cost > cell.cost; }
};

bool CompareWaveCell(const WAVE_CELL &a, const WAVE_CELL &b)
{
  return a > b; //? < or >?
}

// typedef std::map<COORDINATE, WAVE_CELL> WAVEFRONT;
typedef std::vector<COORDINATE> PATH;
typedef std::set<COORDINATE> SETCOO;
typedef std::map<size_t, PATH> ROUTING;
typedef std::map<COORDINATE, int> COSTS;
// typedef priority_queue<WAVE_CELL, vector<WAVE_CELL>, greater<WAVE_CELL>> WAVEFRONT;

std::vector<WAVE_CELL>::iterator findWaveCellInWaveFront(COORDINATE coord, vector<WAVE_CELL> &WAVEFRONT)
{
  std::vector<WAVE_CELL>::iterator iter = WAVEFRONT.begin();
  int size = WAVEFRONT.size();
  int iternumber = size;

  for (int i = 0; i < size; i++)
  {
    if (WAVEFRONT[i].coord == coord)
    {
      iternumber = i;
    }
  }
  iter = WAVEFRONT.begin() + iternumber;
  return iter;
}

int get_cost(const COSTS &costs, const COORDINATE &coo)
{
  COSTS::const_iterator it = costs.find(coo);
  assert(it != costs.end());

  if (it == costs.end())
    return 0; //?

  return it->second;
}

bool is_blocked(const COSTS &costs, const COORDINATE &coo)
{
  return get_cost(costs, coo) == COST_BLOCKED;
}

void block(COSTS &costs, const COORDINATE &coo) { costs[coo] = COST_BLOCKED; }

void unblock(COSTS &costs, const COORDINATE &coo)
{
  COSTS::iterator it = costs.find(coo);
  if (it != costs.end() && it->second != COST_BLOCKED)
    return;

  it->second = 1;
}

/////////////////////////////////////////////////////////////////////////////////////////////

struct PENALTY
{
  size_t via;
  size_t bend;
  size_t layer1_vertical;
  size_t layer1_horizontal;
  size_t layer2_vertical;
  size_t layer2_horizontal;
};

void read_specs(const string &path, const string &name, NETLIST &netlist,
                COORDINATE &dims, COSTS &costs, PENALTY &penalty)
{
  string p = path;

  string path_netlist = p + "/" + (name + ".nl");
  string path_grid = p + "/" + (name + ".grid");

  ifstream in_netlist(path_netlist);

  if (in_netlist)
  {
    size_t nof_nets;
    in_netlist >> nof_nets;

    netlist.reserve(nof_nets);

    for (size_t n = 0; n < nof_nets; n++)
    {
      NETCOORDINATES pcs;

      in_netlist >> pcs.netid;

      in_netlist >> pcs.coo1.layer;
      pcs.coo1.layer--;

      in_netlist >> pcs.coo1.x;
      in_netlist >> pcs.coo1.y;

      in_netlist >> pcs.coo2.layer;
      pcs.coo2.layer--;

      in_netlist >> pcs.coo2.x;
      in_netlist >> pcs.coo2.y;

      netlist.push_back(pcs);
    }
  }

  ifstream in_grid(path_grid);

  if (in_grid)
  {
    size_t size_x, size_y;
    in_grid >> size_x >> size_y >> penalty.bend >> penalty.via;

    dims.layer = 2;
    dims.x = size_x;
    dims.y = size_y;

    for (size_t layer = 0; layer < 2; layer++)
      for (size_t y = 0; y < size_y; y++)
        for (size_t x = 0; x < size_x; x++)
        {
          int cost;
          in_grid >> cost;

          COORDINATE coo;
          coo.layer = layer;
          coo.x = x;
          coo.y = y;
          costs[coo] = cost;
        }

    std::cout << "done" << std::endl;
  }
}

void out_results(ofstream &out, ROUTING &r)
{
  out << r.size() << std::endl;
  for (ROUTING::value_type v : r)
  {
    out << v.first << std::endl;

    for (COORDINATE coo : v.second)
    {
      out << coo.layer + 1;
      out << " ";
      out << coo.x;
      out << " ";
      out << coo.y;
      out << std::endl;
    }

    out << 0 << std::endl;
  }
}

//////////////////////////////////////////////////////////////////////////////////////

bool is_valid(const COORDINATE &dims, const COORDINATE &coo)
{
  return (coo.layer < dims.layer && coo.x < dims.x && coo.y < dims.y);
}

COORDINATE back_trace(DIRECTION d, const COORDINATE &coo)
{
  COORDINATE coo_prev = coo;

  switch (d)
  {
  case N:
    coo_prev.y--;
    break;

  case S:
    coo_prev.y++;
    break;

  case E:
    coo_prev.x--;
    break;

  case W:
    coo_prev.x++;
    break;

  case U:
    coo_prev.layer--;
    break;

  case D:
    coo_prev.layer++;
    break;
  }

  return coo_prev;
}

bool next_coo(DIRECTION d, const COORDINATE &dims, const COORDINATE coo,
              COORDINATE &coo_next)
{
  coo_next = coo;

  switch (d)
  {
  case N:
    if (coo.y + 1 >= dims.y)
      return false;

    coo_next.y++;
    break;

  case S:
    if (coo.y == 0)
      return false;

    coo_next.y--;
    break;

  case E:
    if (coo.x + 1 >= dims.x)
      return false;

    coo_next.x++;
    break;

  case W:
    if (coo.x == 0)
      return false;

    coo_next.x--;
    break;

  case U:
    if (coo.layer + 1 >= dims.layer)
      return false;

    coo_next.layer++;
    break;

  case D:
    if (coo.layer == 0)
      return false;

    coo_next.layer--;
    break;
  }

  return true;
}

////////////////////////////////////////////////////////////////////////////////////////////

void route_single_route(const COORDINATE &dims, const COSTS &grid,
                        const NETCOORDINATES &pcs, const PENALTY &penalty,
                        PATH &path, long long &cost, bool lee)
{
  // std::queue<COORDINATE> qq;
  // WAVEFRONT wf;
  // set<COORDINATE> EXPANDED;

  map<COORDINATE, DIRECTION> PRED;

  vector<vector<vector<bool>>> EXPANDED(dims.layer, vector<vector<bool>>(dims.x, vector<bool>(dims.y, false)));

  vector<vector<vector<bool>>> inWavefront(dims.layer, vector<vector<bool>>(dims.x, vector<bool>(dims.y, false)));

  COORDINATE const coo_start = pcs.coo1;
  COORDINATE const coo_end = pcs.coo2;

  // qq.push(coo_start);  // push in start

  WAVE_CELL wc;

  COSTS::const_iterator itC = grid.find(coo_start);
  assert(itC != grid.end());

  wc.cost = itC->second;
  wc.pred = NONE;
  wc.coord = coo_start;

  vector<WAVE_CELL> WAVEFRONT; //! minheap of wavefront
  WAVEFRONT.emplace_back(wc);

  inWavefront[wc.coord.layer][wc.coord.x][wc.coord.y] = true;

  bool reach = false;

  while (!reach)
  {
    if (WAVEFRONT.empty())
    {
      break;
    }

    // make_heap(WAVEFRONT.begin(), WAVEFRONT.end(), CompareWaveCell);
    WAVE_CELL cellToBeExpand = WAVEFRONT.front();
    WAVEFRONT.erase(WAVEFRONT.begin());

    inWavefront[cellToBeExpand.coord.layer][cellToBeExpand.coord.x][cellToBeExpand.coord.y] = false;

    EXPANDED[cellToBeExpand.coord.layer][cellToBeExpand.coord.x][cellToBeExpand.coord.y] = true;

    if (cellToBeExpand.coord == coo_end)
    {
      reach = true; // no need for reach?
      cout << WAVEFRONT.size() << endl;
      cost += cellToBeExpand.cost;
      break;
    }

    for (DIRECTION d : allowed_expansion)
    {
      COORDINATE coo_next;
      if (next_coo(d, dims, cellToBeExpand.coord, coo_next))
      {
        assert(is_valid(dims, coo_next));
        if (EXPANDED[coo_next.layer][coo_next.x][coo_next.y]) //! UNEXPANDED NEIGHBOR
        {
          continue;
        }
        if (is_blocked(grid, coo_next))
        {
          continue;
        }
        size_t cost_to_move;

        cost_to_move = get_cost(grid, coo_next);

        if (is_bend(cellToBeExpand.pred, d))
          cost_to_move += penalty.bend;
        if (is_via(d))
          cost_to_move += penalty.via;

        if (d == N || d == S)
        {
          if (coo_next.layer == 0)
          {
            cost_to_move += penalty.layer1_vertical;
          }
          else if (coo_next.layer == 1)
          {
            cost_to_move += penalty.layer2_vertical;
          }
        }
        if (d == E || d == W)
        {
          if (coo_next.layer == 0)
          {
            cost_to_move += penalty.layer1_horizontal;
          }
          else if (coo_next.layer == 1)
          {
            cost_to_move += penalty.layer2_horizontal;
          }
        }

        size_t path_cost_new = cellToBeExpand.cost + cost_to_move;

        if (!inWavefront[coo_next.layer][coo_next.x][coo_next.y])
        {
          WAVE_CELL cellToBeAdd;
          cellToBeAdd.coord = coo_next;
          cellToBeAdd.pred = d;
          cellToBeAdd.cost = path_cost_new;
          PRED[cellToBeAdd.coord] = d;
          WAVEFRONT.emplace_back(cellToBeAdd);
          if (!lee)
          {
            push_heap(WAVEFRONT.begin(), WAVEFRONT.end(), CompareWaveCell);
          }
          inWavefront[coo_next.layer][coo_next.x][coo_next.y] = true;
        }
        else
        {
          vector<WAVE_CELL>::iterator iter = findWaveCellInWaveFront(coo_next, WAVEFRONT);

          assert(iter != WAVEFRONT.end());
          if (iter->cost > path_cost_new)
          {
            iter->pred = d;
            iter->cost = path_cost_new;
            PRED[iter->coord] = d;
            if (!lee)
            {
              make_heap(WAVEFRONT.begin(), WAVEFRONT.end(), CompareWaveCell);
            }
          }
        }
      }
    }
  }

  assert(path.empty());
  COORDINATE coo = coo_end;

  if (!EXPANDED[coo_end.layer][coo_end.x][coo_end.y])
  {
    cout << "... No route ";
    return;
  }

  while (!(coo == coo_start))
  {
    path.push_back(coo);

    DIRECTION predforpath = PRED[coo];

    if (is_via(predforpath))
    {
      COORDINATE coo_via = coo;
      coo_via.layer = 2;
      path.push_back(coo_via);
    }
    coo = back_trace(predforpath, coo);

    assert(is_valid(dims, coo));
  }

  path.push_back(coo);
}

/////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char const *argv[])
{
  if (argc < 3)
  {
    cout << "core_router path specsfile" << std::endl;
    return -1;
  }

  NETLIST netlist;
  COORDINATE dims;
  COSTS costs;
  PENALTY penalty;

  bool lee = false;

  cout << argv[1] << " -- " << argv[2] << std::endl;

  read_specs(argv[1], argv[2], netlist, dims, costs, penalty);

  ROUTING routing;
  penalty.layer1_horizontal = 0;
  penalty.layer1_vertical = 0;
  penalty.layer2_horizontal = 0;
  penalty.layer2_vertical = 0;

  if (netlist.size() > 120) // bench5, fract2, primary1, industry1
  {
    penalty.layer1_vertical = penalty.bend;
    penalty.layer2_horizontal = penalty.bend;
  }

  if (netlist.size() > 800)//primary1, industry1
  {
    lee = true;
  }

  long long totalcost = 0;
  clock_t start_time = clock();

  for (NETCOORDINATES &pcs : netlist)
  {
    cout << "routing net " << pcs.netid << "... ";
    PATH path;

    unblock(costs, pcs.coo1);
    unblock(costs, pcs.coo2);

    route_single_route(dims, costs, pcs, penalty, path, totalcost, lee);

    for (COORDINATE coo : path)
    {
      block(costs, coo);
    }
    reverse(path.begin(), path.end());
    routing[pcs.netid] = path;

    cout << "done" << std::endl;
  }

  clock_t end_time = clock();

  cout << "Total cost: " << totalcost << endl
       << "Running time: " << (double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << endl;

  string path_result =
      (std::string(argv[3]) + "/" + std::string(argv[2]) + ".route");
  ofstream out_result(path_result);
  out_results(out_result, routing);

  return 0;
}
