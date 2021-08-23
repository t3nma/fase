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
  long long int motifCount;
  IGtrie igtrie;
  map<int, string> canonicalTypes;
  map<string, int> canonicalCounts;

  int** vext;
  int* vextSz;
  int* vsub;
  double* sampProb;
  char sadjM[MAXMOTIF * MAXMOTIF + 1];
  char nauty_s[MAXMOTIF * MAXMOTIF + 1];

  void reduceCanonicalTypes();
  void expandEnumeration(int depth, int labelNode);
  void expandQueryEnumeration(int, int, Graph*);
  void getSubgraphFrequency(pair< pair<int, long long int>, pair<int, int> > element, Isomorphism* iso);
  void dfsUpdate(int, bool, int, int);
  void dfsUpdateM(int, int);
  void dfsUpdateM2(int, int, bool, int);
  void incrementCount(int, int);

 public:
  Fase(Graph*, bool, int);
  ~Fase();

  int getTypes();
  void runCensus();
  void updateCensus(int, int, bool);
  void monitor(int, int, bool);
  void monitor2(int, int, bool);
  void initSampling(int sz, double* _sampProb);
  long long int getMotifCount() {return motifCount;}
  vector< pair<int,string> > subgraphCount();
  void setQuery(Graph*);
  void setQuery2(Graph*);
  void setup();
};

#endif
