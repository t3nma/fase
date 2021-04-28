#include "Fase.h"

Fase::Fase(Graph* _g, bool _directed, bool _monitor, int _K)
{
  directed = _directed;
  monitor = _monitor;
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

Fase::~Fase()
{
  for (int i = 1; i < MAXMOTIF; i++)
    delete[] vext[i];
  delete[] vext;
  delete[] vextSz;
  delete[] vsub;
  delete[] sampProb;

  igtrie.destroy();
}

void Fase::initSampling(int sz, double* _sampProb)
{
  int i;
  for (i = 0; i < sz; i++)
    sampProb[i] = _sampProb[i];

  sampling = true;
}

void Fase::runCensus()
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

/*
 * Stream update handler.
 * Sort-of-StreamFaSE but taking advantage of
 * the previously known igtrie information.
 */
void Fase::updateCensus(char op, int u, int v)
{
  vector<int> vsub(K), vext;
  bool connected = false;
  map< int, vector<int> > origin;
  int nodeInLabel, nodeExLabel;
  long long int inLabel, exLabel;

  if (op == 'A')
  {
    graph->addEdge(u, v);
    if (!directed)
      graph->addEdge(v, u);
  }

  vsub[0] = u;
  vsub[1] = v;

  for (int n: *graph->outEdges(u))
  {
    if (n == v)
      continue;

    vext.push_back(n);
    origin[n].push_back(u);
  }

  for (int n: *graph->outEdges(v))
  {
    if (n == u)
      continue;

    if (find(vext.begin(), vext.end(), n) == vext.end())
      vext.push_back(n);
    origin[n].push_back(v);
  }

  if (directed && graph->hasEdge(v, u))
    connected = true;

  inLabel = Label::updateLabel(&vsub[0], vsub[1], 1);
  nodeInLabel = igtrie.insertLabel(0, inLabel, Label::repDigits(1), false);
  exLabel = (!directed) ? 0 : 2 * graph->hasEdge(v, u);
  nodeExLabel = igtrie.insertLabel(0, exLabel, Label::repDigits(1), false);

  if (monitor)
  {
    reduceCanonicalTypes();
    char _op = (op == 'A') ? '+' : '-';
    cout << _op << "(" << u+1 << "," << v+1 << "):\n";
  }

  dfsUpdate(vsub, vext, 2, op, connected, origin, nodeInLabel, inLabel, nodeExLabel, exLabel, inLabel, exLabel);

  if (op == 'R')
  {
    graph->rmEdge(u, v);
    if (!directed)
      graph->rmEdge(v, u);
  }
}

void Fase::expandEnumeration(int depth, int labelNode, long long int label)
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

      // keep expanding unless we're not interested
      // in the current igtrie path
      if (clabelNode != -1)
        expandEnumeration(depth + 1, clabelNode, clabel);
    }
  }
}

void Fase::getSubgraphFrequency(pair<long long int, int> element, Isomorphism* iso)
{
  Label::fillNautyMatrix(sadjM, K, element.first);

  nauty_s[0] = '\0';
  iso->canonicalStrNauty(sadjM, nauty_s);
  string str = string(nauty_s);
  labelCanonicalType[element.first] = str;
  canonicalTypes[str] += element.second;
}

/*
 * Recursive runner for setQuery().
 * Starting with 2 fixed nodes it simulates a sort-of-FaSE
 * expansion in order to generate paths in the igtrie.
 */
void Fase::expandQueryEnumeration(vector<int>& vsub, int* vext, int numVext, int depth, int nodeLabel, long long int label, Graph* g)
{
  if (depth == K-1)
  {
    for (int i=0; i!=numVext; ++i)
    {
      label = Label::updateLabel(&vsub[0], vext[i], depth);
      igtrie.insertLabel(nodeLabel, label, Label::repDigits(depth));
    }
    return;
  }

  while (numVext)
  {
    int next = vext[--numVext];
    vsub[depth] = next;

    for (int w: *graph->neighbours(next))
    {
      if (find(vsub.begin(), vsub.begin()+depth, w) != vsub.begin()+depth ||
          find(vext, vext+numVext, w) != vext+numVext)
        continue;

      vext[numVext++] = w;
    }

    long long int cLabel = Label::updateLabel(&vsub[0], next, depth);
    int cNodeLabel = igtrie.insertLabel(nodeLabel, label, Label::repDigits(depth));

    expandQueryEnumeration(vsub, vext, numVext, depth+1, cNodeLabel, cLabel, g);
  }
}

/*
 *
 */
void Fase::dfsUpdate(vector<int>& vsub, vector<int>& vext, int depth, char op, bool connected, map<int, vector<int> >& origin,
                     int nodeInLabel, long long int inLabel, int nodeExLabel, long long int exLabel, long long int inPath, long long int exPath)
{
  if (depth == K)
  {
    if (op == 'A')
    {
      if (nodeInLabel != -1)
      {
        igtrie.incrementLabel(nodeInLabel, 1);
        if (monitor)
          // TODO: check inPath against values of K > LB_WORD_LEN
          cout << "new occurrence of " << labelCanonicalType[inPath] << "\n";
      }
      if (connected && nodeExLabel != -1)
        igtrie.incrementLabel(nodeExLabel, -1);
    }
    else
    {
      if (nodeInLabel != -1)
        igtrie.incrementLabel(nodeInLabel, -1);
      if (connected && nodeExLabel != -1)
      {
        igtrie.incrementLabel(nodeExLabel, 1);
        if (monitor)
          // TODO: check exPath against values of K > LB_WORD_LEN
          cout << "new occurrence of " << labelCanonicalType[exPath] << "\n";
      }
    }

    return;
  }

  while (!vext.empty())
  {
    int next = vext.back();
    vsub[depth] = next;
    vext.pop_back();

    int search = 0;
    for (int w: origin[next])
    {
      if (w == vsub[0])
        search |= (1 << 0);
      if (w == vsub[1])
        search |= (1 << 1);
    }

    bool _connected = (search == (1 << 2) - 1);

    long long int cInLabel = inLabel;
    int cNodeInLabel = nodeInLabel;
    long long int cInPath = inPath;
    long long int cExLabel = exLabel;
    int cNodeExLabel = nodeExLabel;
    long long int cExPath = exPath;

    if (nodeInLabel != -1)
    {
      cInLabel = Label::updateLabel(&vsub[0], vsub[depth], depth);
      cNodeInLabel = igtrie.insertLabel(nodeInLabel, cInLabel, Label::repDigits(depth), false);
      cInPath = (inPath << depth) | cInLabel;
    }

    if (nodeExLabel != -1)
    {
      cExLabel = Label::updateLabel(&vsub[0], vsub[depth], depth);
      cNodeExLabel = igtrie.insertLabel(nodeExLabel, cExLabel, Label::repDigits(depth), false);
      cExPath = (exPath << depth) | cExLabel;
    }

    // we skip the current iteration unless:
    // . we are interested in the connected path
    // OR
    // . we're in a connected enumeration AND we are interested
    //   in the disconnected path
    if (cNodeInLabel == -1 && (!_connected || cNodeExLabel == -1))
      continue;

    for (int w: *graph->neighbours(next))
    {
      int i;
      for (i=0; i!=depth; ++i)
        if (w == vsub[i] || graph->isConnected(w, vsub[i]))
          break;

      if (i == depth)
      {
        vext.push_back(w);

        vector<int> orig { next };
        if (search & (1 << 0))
          orig.push_back(vsub[0]);
        if (search & (1 << 1))
          orig.push_back(vsub[1]);

        set_union(orig.begin(), orig.end(),
                  origin[next].begin(), origin[next].end(),
                  origin[w].begin());
      }
    }

    dfsUpdate(vsub, vext, depth+1, op, _connected, origin, cNodeInLabel, cInLabel, cNodeExLabel, cExLabel, cInPath, cExPath);
  }
}

void Fase::reduceCanonicalTypes()
{
  if (!canonicalTypes.empty())
    return;

  Isomorphism *iso = new Isomorphism();
  iso->initNauty(K, directed);
  for (auto element : igtrie.enumerate(K))
    getSubgraphFrequency(element, iso);
  iso->finishNauty();
}

int Fase::getTypes()
{
  reduceCanonicalTypes();
  return (int)canonicalTypes.size();
}

vector<pair<int, string> > Fase::subgraphCount()
{
  if (canonicalTypes.empty())
    reduceCanonicalTypes();
  else
  {
    for (auto it=canonicalTypes.begin(); it!=canonicalTypes.end(); ++it)
      it->second = 0;

    for (auto element : igtrie.enumerate(K))
      canonicalTypes[labelCanonicalType[element.first]] += element.second;
  }

  vector<pair<int, string> > subgraphVector;
  for (auto element : canonicalTypes)
    subgraphVector.push_back(make_pair(element.second, element.first));

  sort(subgraphVector.begin(), subgraphVector.end());
  reverse(subgraphVector.begin(), subgraphVector.end());

  return subgraphVector;
}

/*
 * Insert a given subgraph in the igtrie by
 * enumerating all possible paths starting with
 * any pair of vertex.
 */
void Fase::setQuery(Graph *g)
{
  vector<int> vsub(K);
  int nodeLabel, *vext, numVext;
  long long int label;

  Label::init(g, directed);

  vext = new int[K];
  numVext = 0;

  for (int u=0; u!=K; ++u)
    for (int v=0; v!=K; ++v)
    {
      if (u == v)
        continue;

      vsub[0] = u;
      vsub[1] = v;

      for (int w: *g->neighbours(u))
        if (w != v)
          vext[numVext++] = w;

      for (int w: *g->neighbours(v))
        if (w != u && find(vext, vext+numVext, w) == vext+numVext)
          vext[numVext++] = w;

      label = Label::updateLabel(&vsub[0], vsub[1], 1);
      nodeLabel = igtrie.insertLabel(0, label, Label::repDigits(1));

      expandQueryEnumeration(vsub, vext, numVext, 2, nodeLabel, label, g);

      numVext = 0;
    }

  delete[] vext;
}
