#include "Fase.h"

Fase::Fase(Graph* _g, bool _directed, int _K)
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

void Fase::setup()
{
  motifCount = 0;
  reduceCanonicalTypes();
  Label::init(graph, directed);
}

void Fase::runCensus()
{
  for (int i = 0; i < graph->numNodes(); i++)
    if (!sampling || Random::testProb(sampProb[0]))
    {
      vsub[0] = i;

      vextSz[1] = 0;
      for (int w: *graph->neighbours(i))
        if (w > i)
          vext[1][vextSz[1]++] = w;

      expandEnumeration(1, 0);
    }
}

/*
 * Stream update handler.
 * Sort-of-StreamFaSE but taking advantage of
 * the previously known igtrie information.
 */
void Fase::updateCensus(int u, int v, bool increment)
{
  int nodeInLabel, nodeExLabel;
  long long int inLabel, exLabel;

  vsub[0] = u;
  vsub[1] = v;

  if (increment)
  {
    exLabel = Label::updateLabel(vsub, v, 1);

    graph->addEdge(u, v);
    if (!directed)
      graph->addEdge(v, u);

    inLabel = Label::updateLabel(vsub, v, 1);
  }
  else
  {
    inLabel = Label::updateLabel(vsub, v, 1);

    graph->rmEdge(u, v);
    if (!directed)
      graph->rmEdge(v, u);

    exLabel = Label::updateLabel(vsub, v, 1);
  }

  nodeInLabel = igtrie.insertLabel(0, inLabel, Label::repDigits(1), false);
  nodeExLabel = igtrie.insertLabel(0, exLabel, Label::repDigits(1), false);

  vextSz[2] = 0;

  for (int w: *graph->neighbours(u))
    if (w != v)
      vext[2][vextSz[2]++] = w;

  for (int w: *graph->neighbours(v))
  {
    if (w == u)
      continue;

    if (find(vext[2], vext[2]+vextSz[2], w) == vext[2]+vextSz[2])
      vext[2][vextSz[2]++] = w;
  }

  dfsUpdate(2, increment, nodeInLabel, nodeExLabel);
}

/*
 * Stream update handler for monitor mode.
 */
void Fase::monitor(int u, int v, bool increment)
{
  int nodeLabel;
  long long int label;

  // reset subgraph counters
  motifCount = 0;
  for (auto it=canonicalCounts.begin(); it!=canonicalCounts.end(); ++it)
    it->second = 0;

  if (increment)
  {
    graph->addEdge(u, v);
    if (!directed)
      graph->addEdge(v, u);
  }
  else
  {
    graph->rmEdge(u, v);
    if (!directed)
      graph->rmEdge(v, u);
  }

  vsub[0] = u;
  vsub[1] = v;

  vextSz[2] = 0;

  for (int w: *graph->neighbours(u))
    if (w != v)
      vext[2][vextSz[2]++] = w;

  for (int w: *graph->neighbours(v))
  {
    if (w == u)
      continue;

    if (find(vext[2], vext[2]+vextSz[2], w) == vext[2]+vextSz[2])
      vext[2][vextSz[2]++] = w;
  }

  label = Label::updateLabel(vsub, v, 1);
  nodeLabel = igtrie.insertLabel(0, label, Label::repDigits(1), false);

  dfsUpdateM(2, nodeLabel);
}

void Fase::monitor2(int u, int v, bool increment)
{
  if (increment)
  {
    monitor(u, v, increment);
    return;
  }

  // reset subgraph counters
  motifCount = 0;
  for (auto it=canonicalCounts.begin(); it!=canonicalCounts.end(); ++it)
    it->second = 0;

  graph->rmEdge(u, v);
  if (!directed)
    graph->rmEdge(v, u);

  if (graph->neighbours(u)->size() <= graph->neighbours(v)->size())
  {
    vsub[0] = u;
    vsub[K-1] = v;
  }
  else
  {
    vsub[0] = v;
    vsub[K-1] = u;
  }

  vextSz[1] = 0;
  for (int w: *graph->neighbours(vsub[0]))
    if (w != vsub[K-1])
      vext[1][vextSz[1]++] = w;

  dfsUpdateM2(1, vsub[K-1], false, 0);
}

void Fase::incrementCount(int labelNode, int value)
{
  canonicalCounts[ canonicalTypes[labelNode] ] += value;
  motifCount += value;
}

void Fase::expandEnumeration(int depth, int labelNode)
{
  if (igtrie.isFinal(labelNode))
    incrementCount(labelNode, 1);

  if (depth == K - 1)
  {
    long long int clabel;
    int currentVertex, clabelNode;

    while (vextSz[depth])
    {
      currentVertex = vext[depth][--vextSz[depth]];

      clabel = Label::updateLabel(vsub, currentVertex, depth);
      clabelNode = igtrie.insertLabel(labelNode, clabel, Label::repDigits(depth), false);

      if (clabelNode != -1)
        incrementCount(clabelNode, 1);
    }

    return;
  }

  int i, j, currentVertex, clabelNode = -1;
  long long int clabel;

  for (i = 0; i < vextSz[depth]; i++)
    vext[depth + 1][i] = vext[depth][i];

  while (vextSz[depth])
  {
    currentVertex = vext[depth][--vextSz[depth]];

    vextSz[depth + 1] = vextSz[depth];
    vsub[depth] = currentVertex;

    clabel = Label::updateLabel(vsub, currentVertex, depth);
    clabelNode = igtrie.insertLabel(labelNode, clabel, Label::repDigits(depth), false);

    // keep expanding unless we're not interested
    // in the current path
    if (clabelNode == -1)
      continue;

    for (int w: *graph->neighbours(currentVertex))
    {
      if (w <= vsub[0])
        continue;

      for (j = 0; j < depth; j++)
        if (w == vsub[j] || graph->isConnected(w, vsub[j]))
          break;

      if (j == depth)
        vext[depth + 1][vextSz[depth + 1]++] = w;
    }

    expandEnumeration(depth + 1, clabelNode);
  }
}

void Fase::getSubgraphFrequency(pair< pair<int, long long int>, pair<int, int> > element, Isomorphism* iso)
{
  Label::fillNautyMatrix(sadjM, element.second.first, element.first.second);

  nauty_s[0] = '\0';
  iso->canonicalStrNauty(sadjM, nauty_s);
  string str = string(nauty_s);
  canonicalTypes[element.first.first] = str;
  canonicalCounts[str] = 0;
}

/*
 * Recursive runner for setQuery().
 * Starting with 2 fixed nodes it simulates a sort-of-FaSE
 * expansion in order to generate paths in the igtrie.
 */
void Fase::expandQueryEnumeration(int depth, int nodeLabel, Graph* g)
{
  if (depth == g->numNodes()-1)
  {
    int next, cNodeLabel;
    long long int label;

    while (vextSz[depth])
    {
      next = vext[depth][--vextSz[depth]];
      vsub[depth] = next;

      if (next == -1)
        continue;

      label = Label::updateLabel(vsub, next, depth);
      cNodeLabel = igtrie.insertLabel(nodeLabel, label, Label::repDigits(depth));

      igtrie.setFinal(cNodeLabel);
    }
    return;
  }

  int next, cNodeLabel;
  long long int label;

  for (int i=0; i!=vextSz[depth]; ++i)
    vext[depth+1][i] = vext[depth][i];
  vextSz[depth+1] = vextSz[depth];

  for (int i=0; i!=vextSz[depth]; ++i)
  {
    next = vext[depth][i];

    if (next == -1)
      continue;

    vsub[depth] = next;
    vext[depth+1][i] = -1;

    for (int w: *g->neighbours(next))
    {
      int j;
      for (j=0; j!=depth; ++j)
        if (w == vsub[j] || g->isConnected(w, vsub[j]))
          break;

      if (j == depth)
        vext[depth+1][vextSz[depth+1]++] = w;
    }

    label = Label::updateLabel(vsub, next, depth);
    cNodeLabel = igtrie.insertLabel(nodeLabel, label, Label::repDigits(depth));

    expandQueryEnumeration(depth+1, cNodeLabel, g);

    vext[depth+1][i] = next;
    vextSz[depth+1] = vextSz[depth];
  }
}

void Fase::dfsUpdate(int depth, bool increment, int nodeInLabel, int nodeExLabel)
{
  if (nodeInLabel != -1 && igtrie.isFinal(nodeInLabel))
    incrementCount(nodeInLabel, 1 - 2*!increment);

  if (nodeExLabel != -1 && igtrie.isFinal(nodeExLabel))
    incrementCount(nodeExLabel, 1 - 2*increment);

  if (depth == K-1)
  {
    int next, cNodeInLabel, cNodeExLabel;
    long long int cLabel;

    while (vextSz[depth])
    {
      next = vext[depth][--vextSz[depth]];
      cLabel = Label::updateLabel(vsub, next, depth);

      if (nodeInLabel != -1)
      {
        cNodeInLabel = igtrie.insertLabel(nodeInLabel, cLabel, Label::repDigits(depth), false);

        if (cNodeInLabel != -1)
          incrementCount(cNodeInLabel, 1 - 2*!increment);
      }

      if (nodeExLabel != -1)
      {
        cNodeExLabel = igtrie.insertLabel(nodeExLabel, cLabel, Label::repDigits(depth), false);

        if (cNodeExLabel != -1)
          incrementCount(cNodeExLabel, 1 - 2*increment);
      }
    }

    return;
  }

  int next, cNodeInLabel = nodeInLabel, cNodeExLabel = nodeExLabel;
  long long int cLabel;

  for (int i=0; i!=vextSz[depth]; ++i)
    vext[depth+1][i] = vext[depth][i];

  while (vextSz[depth])
  {
    next = vext[depth][--vextSz[depth]];
    vsub[depth] = next;
    vextSz[depth+1] = vextSz[depth];

    cLabel = Label::updateLabel(vsub, next, depth);

    if (nodeInLabel != -1)
      cNodeInLabel = igtrie.insertLabel(nodeInLabel, cLabel, Label::repDigits(depth), false);

    if (nodeExLabel != -1)
      cNodeExLabel = igtrie.insertLabel(nodeExLabel, cLabel, Label::repDigits(depth), false);

    // we skip the current iteration unless we are interested
    // in at least one of the enumeration paths
    if (cNodeInLabel == -1 && cNodeExLabel == -1)
      continue;

    for (int w: *graph->neighbours(next))
    {
      int i;
      for (i=0; i!=depth; ++i)
        if (w == vsub[i] || graph->isConnected(w, vsub[i]))
          break;

      if (i == depth)
        vext[depth+1][vextSz[depth+1]++] = w;
    }

    dfsUpdate(depth+1, increment, cNodeInLabel, cNodeExLabel);
  }
}

void Fase::dfsUpdateM(int depth, int nodeLabel)
{
  if (igtrie.isFinal(nodeLabel))
    incrementCount(nodeLabel, 1);

  if (depth == K-1)
  {
    int next, cNodeLabel;
    long long int cLabel;

    while (vextSz[depth])
    {
      next = vext[depth][--vextSz[depth]];

      cLabel = Label::updateLabel(vsub, next, depth);
      cNodeLabel = igtrie.insertLabel(nodeLabel, cLabel, Label::repDigits(depth), false);

      if (cNodeLabel != -1)
        incrementCount(cNodeLabel, 1);
    }

    return;
  }

  int next, cNodeLabel = nodeLabel;
  long long int cLabel;

  for (int i=0; i!=vextSz[depth]; ++i)
    vext[depth+1][i] = vext[depth][i];

  while (vextSz[depth])
  {
    next = vext[depth][--vextSz[depth]];
    vsub[depth] = next;
    vextSz[depth+1] = vextSz[depth];

    cLabel = Label::updateLabel(vsub, next, depth);

    if (nodeLabel != -1)
      cNodeLabel = igtrie.insertLabel(nodeLabel, cLabel, Label::repDigits(depth), false);

    // we skip the current iteration unless we are
    // interested in the enumeration path
    if (cNodeLabel == -1)
      continue;

    for (int w: *graph->neighbours(next))
    {
      int i;
      for (i=0; i!=depth; ++i)
        if (w == vsub[i] || graph->isConnected(w, vsub[i]))
          break;

      if (i == depth)
        vext[depth+1][vextSz[depth+1]++] = w;
    }

    dfsUpdateM(depth+1, cNodeLabel);
  }
}

void Fase::dfsUpdateM2(int depth, int searchNode, bool connected, int nodeLabel)
{
  if (connected && igtrie.isFinal(nodeLabel))
    incrementCount(nodeLabel, 1);

  if (depth == K-1)
  {
    int next, cNodeLabel;
    long long int cLabel;

    while (vextSz[depth])
    {
      next = vext[depth][--vextSz[depth]];

      // proceed only if the subgraph is connected or
      // the next node is the one that connects it
      if (!connected && next != searchNode)
        continue;

      cLabel = Label::updateLabel(vsub, next, depth);
      cNodeLabel = igtrie.insertLabel(nodeLabel, cLabel, Label::repDigits(depth), false);

      if (cNodeLabel != -1)
        incrementCount(cNodeLabel, 1);

      if (!connected)
        break;
    }

    return;
  }

  int next, cNodeLabel;
  long long int cLabel;
  bool _connected;

  for (int i=0; i!=vextSz[depth]; ++i)
    vext[depth+1][i] = vext[depth][i];

  while (vextSz[depth])
  {
    next = vext[depth][--vextSz[depth]];

    vsub[depth] = next;
    vextSz[depth+1] = vextSz[depth];

    _connected = connected || (next == searchNode);

    cLabel = Label::updateLabel(vsub, next, depth);
    cNodeLabel = igtrie.insertLabel(nodeLabel, cLabel, Label::repDigits(depth), false);

    // we skip the current iteration unless we are
    // interested in the enumeration path
    if (cNodeLabel == -1)
      continue;

    for (int w: *graph->neighbours(next))
    {
      int i;
      for (i=0; i!=depth; ++i)
        if (w == vsub[i] || graph->isConnected(w, vsub[i]))
          break;

      if (i == depth)
        vext[depth+1][vextSz[depth+1]++] = w;
    }

    dfsUpdateM2(depth+1, searchNode, _connected, cNodeLabel);
  }
}

void Fase::reduceCanonicalTypes()
{
  if (!canonicalTypes.empty())
    return;

  Isomorphism *iso = new Isomorphism();

  for (auto element : igtrie.enumerate(K))
  {
    iso->initNauty(element.second.first, directed);
    getSubgraphFrequency(element, iso);
    iso->finishNauty();
  }
}

int Fase::getTypes()
{
  reduceCanonicalTypes();
  return (int)canonicalCounts.size();
}

vector< pair<int, string> > Fase::subgraphCount()
{
  vector< pair<int, string> > ret;

  for (auto it = canonicalCounts.begin(); it != canonicalCounts.end(); ++it) {
    if (it->second > 0)
      ret.push_back(make_pair(it->second, it->first));
  }

  sort(ret.begin(), ret.end());
  reverse(ret.begin(), ret.end());

  return ret;
}

/*
 * Insert a given subgraph in the igtrie by
 * enumerating all possible paths starting with
 * any pair of vertex.
 */
void Fase::setQuery(Graph *g)
{
  int numNodes, nodeLabel;
  long long int label;

  Label::init(g, directed);

  numNodes = g->numNodes();
  for (int u=0; u!=numNodes; ++u)
    for (int v=0; v!=numNodes; ++v)
    {
      if (u == v)
        continue;

      vsub[0] = u;
      vsub[1] = v;

      vextSz[2] = 0;

      for (int w: *g->neighbours(u))
        if (w != v)
          vext[2][vextSz[2]++] = w;

      for (int w: *g->neighbours(v))
      {
        if (w == u)
          continue;

        if (find(vext[2], vext[2]+vextSz[2], w) == vext[2]+vextSz[2])
          vext[2][vextSz[2]++] = w;
      }

      label = Label::updateLabel(vsub, v, 1);
      nodeLabel = igtrie.insertLabel(0, label, Label::repDigits(1));

      expandQueryEnumeration(2, nodeLabel, g);
    }
}

void Fase::setQuery2(Graph *g)
{
  Label::init(g, directed);

  for (int u=0; u!=g->numNodes(); ++u)
  {
    vsub[0] = u;

    vextSz[1] = 0;
    for (int w: *g->neighbours(u))
      vext[1][vextSz[1]++] = w;

    expandQueryEnumeration(1, 0, g);
  }
}
