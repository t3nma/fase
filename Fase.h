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
  map<long long int, string> labelCanonicalType;
  map<string, int> canonicalTypes;

  int** vext;
  int* vextSz;
  int* vsub;
  double* sampProb;
  char sadjM[MAXMOTIF * MAXMOTIF + 1];
  char nauty_s[MAXMOTIF * MAXMOTIF + 1];

  void reduceCanonicalTypes();
  void expandEnumeration(int depth, int labelNode, long long int label);
  void expandQueryEnumeration(int, int, Graph*);
  void getSubgraphFrequency(pair<pair<long long int, int>, int> element, Isomorphism* iso);
  void dfsUpdate(int, bool, int, long long int, int, long long int);
  void dfsUpdateM(int, int, long long int);
  void dfsUpdateM2(int, int, bool, int, long long int);

 public:
  Fase(Graph*, bool, int);
  ~Fase();

  int getTypes();
  void runCensus();
  void setupMonitor();
  void updateCensus(int, int, bool);
  void monitor(int, int, bool);
  void monitor2(int, int, bool);
  void initSampling(int sz, double* _sampProb);
  int getMotifCount() {return motifCount;}
  vector< pair<string,int> > subgraphCount(bool);
  void setQuery(Graph*);
  void setQuery2(Graph*);
};

#endif
