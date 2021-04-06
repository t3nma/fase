#ifndef _MYFASE_
#define _MYFASE_

#include "Common.h"
#include "Graph.h"
#include "Random.h"
#include "Label.h"
#include "MyIGtrie.h"
#include "Isomorphism.h"

class MyFase
{
 private:
  bool directed;
  bool sampling;
  Graph *graph;
  int K;
  int motifCount;
  MyIGtrie igtrie;
  map<string, int> canonicalTypes;

  int** vext;
  int* vextSz;
  int* vsub;
  double* sampProb;
  char sadjM[MAXMOTIF * MAXMOTIF + 1];
  char nauty_s[MAXMOTIF * MAXMOTIF + 1];

  void reduceCanonicalTypes();
  void expandEnumeration(int depth, int labelNode, long long int label);
  void getSubgraphFrequency(pair<long long int, int> element, Isomorphism* iso);

  void genQueryVersions(int, vector<int>&, int*, int**, int*, Graph*, vector< vector<int> >&);
  void expandQuery(Graph*, int*);

 public:
  MyFase(Graph*, bool, int);
  ~MyFase();

  int getTypes();
  void runCensus();
  void initSampling(int sz, double* _sampProb);
  int getMotifCount() {return motifCount;}
  vector<pair<int, string> > subgraphCount();

  void setQuery(Graph*);
  void igtrieEnum();
};

#endif
