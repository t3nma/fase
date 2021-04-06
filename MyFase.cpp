#include "MyFase.h"

MyFase::MyFase(Graph* _g, bool _directed, int _K)
{
  directed = _directed;
  graph = _g;
  sampling = false;
  K = _K;

  vext = new int*[MAXMOTIF];
  for (int i = 1; i < MAXMOTIF; i++)
    vext[i] = new int[graph->numNodes()];

  vextSz = new int[MAXMOTIF];
  vsub = new int[MAXMOTIF];
  sampProb = new double[MAXMOTIF];

  igtrie.init(K);
}

MyFase::~MyFase()
{
  for (int i = 1; i < MAXMOTIF; i++)
    delete[] vext[i];
  delete[] vext;
  delete[] vextSz;
  delete[] vsub;
  delete[] sampProb;

  igtrie.destroy();
}

void MyFase::initSampling(int sz, double* _sampProb)
{
  int i;
  for (i = 0; i < sz; i++)
    sampProb[i] = _sampProb[i];

  sampling = true;
}

void MyFase::runCensus()
{
  motifCount = 0;

  Label::init(graph, directed);

  for (int i = 0; i < graph->numNodes(); i++)
    if (!sampling || Random::testProb(sampProb[0]))
    {
      vsub[0] = i;
      int *nei = graph->arrayNeighbours(i);
      int neiNum = graph->numNeighbours(i);

      vextSz[1] = 0;
      for (int j = 0; j < neiNum; j++)
        if (nei[j] > i)
          vext[1][vextSz[1]++] = nei[j];

      expandEnumeration(1, 0, 0LL);
    }
}

void MyFase::expandEnumeration(int depth, int labelNode, long long int label)
{
  if (depth == K - 1)
  {
    while (vextSz[depth])
    {
      int currentVertex = vext[depth][--vextSz[depth]];

      if (!sampling || Random::testProb(sampProb[depth]))
      {
        long long int clabel = Label::updateLabel(vsub, currentVertex, depth);
        int clabelNode = igtrie.insertLabel(labelNode, clabel, Label::repDigits(depth), false);

        if (clabelNode != -1)
        {
          igtrie.incrementLabel(clabelNode, 1);
          motifCount++;
        }
      }
    }

    return;
  }

  int i, j;
  long long int clabel = label;
  int clabelNode = labelNode;

  for (i = 0; i < vextSz[depth]; i++)
    vext[depth + 1][i] = vext[depth][i];

  while (vextSz[depth])
  {
    int currentVertex = vext[depth][--vextSz[depth]];

    if (!sampling || Random::testProb(sampProb[depth]))
    {
      vextSz[depth + 1] = vextSz[depth];
      vsub[depth] = currentVertex;

      int *eExcl = graph->arrayNeighbours(currentVertex);
      int eExclNum = graph->numNeighbours(currentVertex);

      for (i = 0; i < eExclNum; i++)
      {
        if (eExcl[i] <= vsub[0])
          continue;

        for (j = 0; j < depth; j++)
          if (eExcl[i] == vsub[j] || graph->isConnected(eExcl[i], vsub[j]))
            break;

        if (j == depth)
          vext[depth + 1][vextSz[depth + 1]++] = eExcl[i];
      }

      if (depth >= 1)
      {
        clabel = Label::updateLabel(vsub, currentVertex, depth);
        clabelNode = igtrie.insertLabel(labelNode, clabel, Label::repDigits(depth), false);
      }

      if (clabelNode != -1)
        expandEnumeration(depth + 1, clabelNode, clabel);
    }
  }
}

void MyFase::getSubgraphFrequency(pair<long long int, int> element, Isomorphism* iso)
{
  Label::fillNautyMatrix(sadjM, K, element.first);

  nauty_s[0] = '\0';
  iso->canonicalStrNauty(sadjM, nauty_s);
  string str = string(nauty_s);
  canonicalTypes[str] += element.second;
}

void MyFase::genQueryVersions(int depth, vector<int>& vsub, int* used, int** vext, int* vextSz,
                               Graph* g, vector< vector<int> >& acc) {
  if (depth == K-1) {
    while (vextSz[depth]) {
      vsub[depth] = vext[depth][--vextSz[depth]];
      acc.push_back(vsub);
    }
    return;
  }

  for (int i=0; i!=vextSz[depth]; ++i) {
    vext[depth+1][i] = vext[depth][i];
  }

  while (vextSz[depth]) {
    int curNode = vext[depth][--vextSz[depth]];

    vextSz[depth+1] = vextSz[depth];
    vsub[depth] = curNode;
    used[curNode] = 1;

    for (int v: *g->neighbours(curNode)) {
      if (used[v]) {
        continue;
      }

      int i;
      for (i=0; i!=vextSz[depth]; ++i) {
        if (v == vext[depth][i]) {
          break;
        }
      }

      if (i == vextSz[depth]) {
        vext[depth+1][vextSz[depth+1]++] = v;
      }
    }

    genQueryVersions(depth+1, vsub, used, vext, vextSz, g, acc);
    used[curNode] = 0;
  }
}

void MyFase::expandQuery(Graph *g, int *vsub) {
  int depth = 1;
  int nodeLabel = 0;
  long long int label = 0;

  while (depth != K) {
    label = Label::updateLabel(vsub, vsub[depth], depth);
    nodeLabel = igtrie.insertLabel(nodeLabel, label, Label::repDigits(depth));
    depth++;
  }
}

void MyFase::reduceCanonicalTypes()
{
  if (!canonicalTypes.empty())
    return;

  Isomorphism *iso = new Isomorphism();
  iso->initNauty(K, directed);
  for (auto element : igtrie.enumerate(K))
    getSubgraphFrequency(element, iso);
  iso->finishNauty();
}

int MyFase::getTypes()
{
  reduceCanonicalTypes();
  return (int)canonicalTypes.size();
}

vector<pair<int, string> > MyFase::subgraphCount()
{
  reduceCanonicalTypes();

  vector<pair<int, string> > subgraphVector;
  for (auto element : canonicalTypes)
    subgraphVector.push_back(make_pair(element.second, element.first));

  sort(subgraphVector.begin(), subgraphVector.end());
  reverse(subgraphVector.begin(), subgraphVector.end());

  return subgraphVector;
}

void MyFase::setQuery(Graph *g) {
  vector<int> vsub(K);
  vector< vector<int> > perm;
  int used[K] = {0};
  int **vext = new int*[K];
  int vextSz[K];

  Label::init(g, directed);

  for (int i=1; i!=K; ++i) {
    vext[i] = new int[K];
  }

  for (int i=0; i!=K; ++i) {
    vsub[0] = i;
    used[i] = 1;

    vextSz[1] = 0;
    for (int v: *g->neighbours(i)) {
      vext[1][vextSz[1]++] = v;
    }

    genQueryVersions(1, vsub, used, vext, vextSz, g, perm);
    used[i] = 0;
  }

  for (auto p: perm) {
    expandQuery(g, &p[0]);
  }

  for (int i=1; i!=K; ++i) {
    delete vext[i];
  }
  delete vext;

}

void MyFase::igtrieEnum() {
  for (auto p: igtrie.enumerate(K)) {
    cout << p.first << ": " << p.second << "\n";
  }
}
