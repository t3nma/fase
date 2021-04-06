#ifndef _FASE_
#define _FASE_

#include "Common.h"
#include "Graph.h"
#include "Random.h"
#include "Label.h"
#include "IGtrie.h"
#include "Isomorphism.h"

class Fase
{
 private:
  bool directed;
  bool sampling;
  Graph *graph;
  int K;
  int motifCount;
  IGtrie igtrie;
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
  void genQueryVersions(int depth, vector<int>& vsub, int* used, int** vext, int* vextSz,
                        Graph* g, vector< vector<int> >& acc);
  void expandQuery(Graph* g, int* vsub);

 public:
  Fase(Graph*, bool, int);
  ~Fase();

  int getTypes();
  void runCensus();
  void initSampling(int sz, double* _sampProb);
  int getMotifCount() {return motifCount;}
  vector<pair<int, string> > subgraphCount();
  void setQuery(Graph* g);
};

#endif
