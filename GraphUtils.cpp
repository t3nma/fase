/* -------------------------------------------------

//                                                 
//  88888888888           ad88888ba   88888888888  
//  88                   d8"     "8b  88           
//  88                   Y8,          88           
//  88aaaaa  ,adPPYYba,  `Y8aaaaa,    88aaaaa      
//  88"""""  ""     `Y8    `"""""8b,  88"""""      
//  88       ,adPPPPP88          `8b  88           
//  88       88,    ,88  Y8a     a8P  88           
//  88       `"8bbdP"Y8   "Y88888P"   88888888888  
//                                                 
//

Pedro {Paredes, Ribeiro} - DCC/FCUP
----------------------------------------------------
Graph Utilities

Adapted from gtrieScanner - http://www.dcc.fc.up.pt/gtries/

---------------------------------------------------- */

#include "GraphUtils.h"

int *GraphUtils::_degree;
int **GraphUtils::_ds;
int *GraphUtils::_neighbours;

void GraphUtils::readFileTxt(Graph *g, const char *s, bool dir, bool weight, int base) {

  FILE *f = fopen(s, "r");
  if (!f)
    exit(EXIT_FAILURE);

  int i, a, b, c, size, max;
  vector<int> va, vb;

  size = max = 0;
  while (fscanf(f,"%d %d", &a, &b)==2) {
    if (weight) i=fscanf(f,"%d", &c);
    a += 1 - base; b += 1 - base;
    va.push_back(a);
    vb.push_back(b);
    if (a>max) max=a;
    if (b>max) max=b;
    size++;
  }
    
  fclose(f);

  if (dir) g->createGraph(max, DIRECTED);
  else     g->createGraph(max, UNDIRECTED);

  for (i=0; i<size; i++)
    if (va[i] != vb[i])
      g->addEdge(va[i]-1, vb[i]-1);

  g->preProcess();
  va.clear();
  vb.clear();
}

void GraphUtils::strToGraph(Graph *g, const char *s, int size, bool dir) {
  int i,j;

  if (dir)  g->createGraph(size, DIRECTED);
  else      g->createGraph(size, UNDIRECTED);

  for (i=0; i<size; i++)
    for (j=0; j<size; j++)
      if (s[i*size+j]=='1') g->addEdge(i, j);  
}


int GraphUtils::int_compare(const void *a, const void *b) {
  return (*((int *)b)) - (*((int *)a));
}

int GraphUtils::_compare_int(const void *a, const void *b) {
  return *((int *)a) - *((int *)b);
}

int GraphUtils::_compare_int_descending(const void *a, const void *b) {
  return *((int *)b) - *((int *)a);
}

int GraphUtils::_compare_degree(const void *a, const void *b) {
  int n1 = *((int *)a);
  int n2 = *((int *)b);

  if (_neighbours[n1] < _neighbours[n2]) return -1;
  if (_neighbours[n1] > _neighbours[n2]) return +1;

  return 0;  
}


// Order graph by increasing degree, then by increasing neighbour degree sequence
void GraphUtils::orderGraph(Graph *old, Graph *g) {
  int i, j, aux;
  int size= old->numNodes();
  int degree[size];
  int *ds[size];
  int neighbours[size];
  vector<int> *v;

  for (i=0; i<size; i++) {
    degree[i] = old->nodeOutEdges(i) + old->nodeInEdges(i);
    neighbours[i] = old->numNeighbours(i);
    ds[i] = new int[neighbours[i]];        
    v = old->neighbours(i);
    for (j=0; j<neighbours[i]; j++)
      ds[i][j] = old->nodeOutEdges((*v)[j]) + old->nodeOutEdges((*v)[j]);
    qsort(ds[i], old->numNeighbours(i), sizeof(int), _compare_int);
  }
  
  _degree =  degree;
  _ds     = ds;
  _neighbours = neighbours;
  int n[size], r[size];
  for (i=0; i<size; i++) n[i]=i;
  qsort(n, size, sizeof(int), _compare_degree);
  for (i=0; i<size; i++) r[n[i]]=i;

  for (i=0; i<size; i++)
    delete[] ds[i];

  g->createGraph(size, old->type());
  for (i=0; i<size; i++) {
    v = old->outEdges(i);
    aux = v->size();
    for (j=0; j<aux; j++)
      g->addEdge(r[i], r[(*v)[j]]);
    v = old->inEdges(i);
    aux = v->size();
    for (j=0; j<aux; j++)
      g->addEdge(r[(*v)[j]], r[i]);
  }
  
}
