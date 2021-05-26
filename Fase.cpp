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

void Fase::runCensus()
{
  motifCount = 0;
  reduceCanonicalTypes();
  Label::init(graph, directed);

  for (int i = 0; i < graph->numNodes(); i++)
    if (!sampling || Random::testProb(sampProb[0]))
    {
      vsub[0] = i;

      vextSz[1] = 0;
      for (int w: *graph->neighbours(i))
        if (w > i)
          vext[1][vextSz[1]++] = w;

      expandEnumeration(1, 0, 0LL);
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

  if (increment)
  {
    graph->addEdge(u, v);
    if (!directed)
      graph->addEdge(v, u);
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

  inLabel = Label::updateLabel(vsub, v, 1);
  nodeInLabel = igtrie.insertLabel(0, inLabel, Label::repDigits(1), false);
  exLabel = (!directed) ? 0 : 2 * graph->hasEdge(v, u);
  nodeExLabel = igtrie.insertLabel(0, exLabel, Label::repDigits(1), false);

  dfsUpdate(2, increment, nodeInLabel, inLabel, nodeExLabel, exLabel);

  if (!increment)
  {
    graph->rmEdge(u, v);
    if (!directed)
      graph->rmEdge(v, u);
  }
}

/*
 * Stream update handler for monitor mode.
 */
void Fase::monitor(int u, int v, bool increment)
{
  int nodeLabel;
  long long int label;

  reduceCanonicalTypes();
  for (auto& elem: canonicalTypes)
    elem.second = 0;

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

  dfsUpdateM(2, nodeLabel, label);
}

void Fase::monitor2(int u, int v, bool increment)
{
  reduceCanonicalTypes();
  for (auto& elem: canonicalTypes)
    elem.second = 0;

  if (increment)
  {
    monitor(u, v, increment);
    return;
  }

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

  dfsUpdateM2(1, vsub[K-1], false, 0, 0LL);
}

void Fase::expandEnumeration(int depth, int labelNode, long long int label)
{
  if (depth == K - 1)
  {
    long long int clabel;
    int clabelNode;

    while (vextSz[depth])
    {
      int currentVertex = vext[depth][--vextSz[depth]];

      if (!sampling || Random::testProb(sampProb[depth]))
      {
        clabel = Label::updateLabel(vsub, currentVertex, depth);
        clabelNode = igtrie.insertLabel(labelNode, clabel, Label::repDigits(depth), false);

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

      if (depth >= 1)
      {
        clabel = Label::updateLabel(vsub, currentVertex, depth);
        clabelNode = igtrie.insertLabel(labelNode, clabel, Label::repDigits(depth), false);
      }

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
void Fase::expandQueryEnumeration(int depth, int nodeLabel, Graph* g)
{
  if (depth == K-1)
  {
    int next;
    long long int label;

    while (vextSz[depth])
    {
      next = vext[depth][--vextSz[depth]];

      if (next == -1)
        continue;

      label = Label::updateLabel(vsub, next, depth);
      igtrie.insertLabel(nodeLabel, label, Label::repDigits(depth));
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

void Fase::dfsUpdate(int depth, bool increment, int nodeInLabel, long long int inLabel, int nodeExLabel, long long int exLabel)
{
  if (depth == K-1)
  {
    int next, cNodeInLabel, cNodeExLabel;
    long long int cLabel;

    while (vextSz[depth])
    {
      next = vext[depth][--vextSz[depth]];
      vsub[depth] = next;

      cLabel = Label::updateLabel(vsub, next, depth);

      if (nodeInLabel != -1)
      {
        cNodeInLabel = igtrie.insertLabel(nodeInLabel, cLabel, Label::repDigits(depth), false);

        if (cNodeInLabel != -1)
        {
          igtrie.incrementLabel(cNodeInLabel, 1 - 2*!increment);
          motifCount += 1 - 2*!increment;
        }
      }

      if (nodeExLabel != -1)
      {
        cNodeExLabel = igtrie.insertLabel(nodeExLabel, cLabel, Label::repDigits(depth), false);

        if (cNodeExLabel != -1)
        {
          igtrie.incrementLabel(cNodeExLabel, 1 - 2*increment);
          motifCount += 1 - 2*increment;
        }
      }
    }

    return;
  }

  int next, cNodeInLabel, cNodeExLabel;
  long long int cLabel, cInPath, cExPath;

  for (int i=0; i!=vextSz[depth]; ++i)
    vext[depth+1][i] = vext[depth][i];

  while (vextSz[depth])
  {
    next = vext[depth][--vextSz[depth]];
    vsub[depth] = next;
    vextSz[depth+1] = vextSz[depth];

    cLabel = Label::updateLabel(vsub, next, depth);

    if (nodeInLabel != -1)
    {
      cNodeInLabel = igtrie.insertLabel(nodeInLabel, cLabel, Label::repDigits(depth), false);
      cInPath = (inLabel << Label::repDigits(depth)) | cLabel;
    }

    if (nodeExLabel != -1)
    {
      cNodeExLabel = igtrie.insertLabel(nodeExLabel, cLabel, Label::repDigits(depth), false);
      cExPath = (exLabel << Label::repDigits(depth)) | cLabel;
    }

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

    dfsUpdate(depth+1, increment, cNodeInLabel, cInPath, cNodeExLabel, cExPath);
  }
}

void Fase::dfsUpdateM(int depth, int nodeLabel, long long int label)
{
  if (depth == K-1)
  {
    int next, cNodeLabel;
    long long int cLabel, cPath;

    while (vextSz[depth])
    {
      next = vext[depth][--vextSz[depth]];
      cLabel = Label::updateLabel(vsub, next, depth);
      cNodeLabel = igtrie.insertLabel(nodeLabel, cLabel, Label::repDigits(depth), false);

      if (cNodeLabel != -1)
      {
        cPath = (label << Label::repDigits(depth)) | cLabel;
        canonicalTypes[labelCanonicalType[cPath]]++;
      }
    }

    return;
  }

  int next, cNodeLabel;
  long long int cLabel, cPath;

  for (int i=0; i!=vextSz[depth]; ++i)
    vext[depth+1][i] = vext[depth][i];

  while (vextSz[depth])
  {
    next = vext[depth][--vextSz[depth]];
    vsub[depth] = next;
    vextSz[depth+1] = vextSz[depth];

    cLabel = Label::updateLabel(vsub, next, depth);
    cNodeLabel = igtrie.insertLabel(nodeLabel, cLabel, Label::repDigits(depth), false);
    cPath = (label << Label::repDigits(depth)) | cLabel;

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

    dfsUpdateM(depth+1, cNodeLabel, cPath);
  }
}

void Fase::dfsUpdateM2(int depth, int searchNode, bool connected, int nodeLabel, long long int label)
{
  if (depth == K-1)
  {
    int next, cNodeLabel;
    long long int cLabel, cPath;

    while (vextSz[depth])
    {
      next = vext[depth][--vextSz[depth]];

      if (!connected && next != searchNode)
        continue;

      cLabel = Label::updateLabel(vsub, next, depth);
      cNodeLabel = igtrie.insertLabel(nodeLabel, cLabel, Label::repDigits(depth), false);

      if (cNodeLabel != -1)
      {
        cPath = (label << Label::repDigits(depth)) | cLabel;
        canonicalTypes[labelCanonicalType[cPath]]++;
      }

      if (!connected)
        break;
    }

    return;
  }

  int next, cNodeLabel;
  long long int cLabel, cPath;
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
    cPath = (label << Label::repDigits(depth)) | cLabel;

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

    dfsUpdateM2(depth+1, searchNode, _connected, cNodeLabel, cPath);
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

vector< pair<string, int> > Fase::subgraphCount(bool monitor)
{
  // reset and update counts from igtrie if not in monitor mode
  if (!monitor)
  {
    for (auto& elem: canonicalTypes)
      elem.second = 0;

    for (auto elem : igtrie.enumerate(K))
      canonicalTypes[labelCanonicalType[elem.first]] += elem.second;
  }

  /*
  vector<pair<int, string> > subgraphVector;
  for (auto element : canonicalTypes)
    subgraphVector.push_back(make_pair(element.second, element.first));

  sort(subgraphVector.begin(), subgraphVector.end());
  reverse(subgraphVector.begin(), subgraphVector.end());

  return subgraphVector;
  */

  vector< pair<string, int> > ret;
  for (auto elem: canonicalTypes)
    if (!monitor || elem.second > 0)
      ret.push_back({elem.first, elem.second});

  return ret;
}

/*
 * Insert a given subgraph in the igtrie by
 * enumerating all possible paths starting with
 * any pair of vertex.
 */
void Fase::setQuery(Graph *g)
{
  int nodeLabel;
  long long int label;

  Label::init(g, directed);

  for (int u=0; u!=K; ++u)
    for (int v=0; v!=K; ++v)
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

  for (int u=0; u!=K; ++u)
  {
    vsub[0] = u;

    vextSz[1] = 0;
    for (int w: *g->neighbours(u))
      vext[1][vextSz[1]++] = w;

    expandQueryEnumeration(1, 0, g);
  }
}
