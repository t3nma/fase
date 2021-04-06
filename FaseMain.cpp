#include "Common.h"
#include "Fase.h"
#include "MyFase.h"
#include "DynamicGraph.h"
#include "GraphMatrix.h"
#include "GraphUtils.h"
#include "Timer.h"
#include "Graph.h"
#include "Random.h"

using namespace std;

Graph *G;
int K = 0;
double sampProb[MAXMOTIF], prob;
bool dir = false, detailed = false, draw = false, samp = false, largeScale = false;
char ifilename [200];
char ofilename [200];
FILE *outFile;
time_t t_start, t_end;

void init()
{
  printf("  88888888888           ad88888ba   88888888888  \n"
         "  88                   d8\"     \"8b  88           \n"
         "  88                   Y8,          88           \n"
         "  88aaaaa  ,adPPYYba,  `Y8aaaaa,    88aaaaa      \n"
         "  88\"\"\"\"\"  \"\"     `Y8    `\"\"\"\"\"8b,  88\"\"\"\"\"      \n"
         "  88       ,adPPPPP88          `8b  88           \n"
         "  88       88,    ,88  Y8a     a8P  88           \n"
         "  88       `\"8bbdP\"Y8   \"Y88888P\"   88888888888  \n\n"
         "\tVersion: 1.0\n"
         "FaSE - Fast Subgraph Enumeration (with Sampling)\n"
         "\n\n\tPedro {Paredes, Ribeiro} - DCC/FCUP\n\n\n\n\n");
  t_start = time(0);
}

void displayHelp()
{
  printf("------------ FaSE Usage --------------\nMain Settings: ./FASE -s <Subgraph Size> -i <input file> [arguments...]\n\n\tAll commands:\n-h : Displays this help information\n-s <Integer> : Subgraph Size\n-i <Filename> : Name of input file (Format in Readme.txt)\n-d : Directed Subgraph (Default undirected)\n-o : Name of output file (Default is stdout)\n-dt : Detailed Result (Displays all subgraph types and occurrences)\n-ls : Use a large scale representation (default is adjacency matrix)\n-z : Use 0-based input (Suitable for input files starting at node 0)\n-p <P1> <P2> ... <Ps> : Sets the sampling probabilities by depth (note that -s must have been selected first)\n-q : Ignore arguments and prompt input\n--------------------------------------\n");
}

void read(int argc, char **argv)
{
  int E, V, i, check = 0, itera = 0;
  int zeroBased = 1;
  ofilename[0] = '0';
  ofilename[1] = '\0';

  for (i = 1; i < argc; i++)
  {
    if (argv[i][0] != '-')
      continue;
    if (argv[i][1] == 'h')
    {
      displayHelp();
      K = 0;
      return;
    }

    if (argv[i][1] == 'd' && argv[i][2] == 't')
      detailed = true;
    else if (argv[i][1] == 'd' && argv[i][2] == 'r')
      draw = true;
    else if (argv[i][1] == 'd')
      dir = true;

    if (argv[i][1] == 'l' && argv[i][2] == 's')
      largeScale = true;

    if (argv[i][1] == 'z')
      zeroBased = 0;

    if (argv[i][1] == 'i')
    {
      i++;
      strcpy(ifilename, argv[i]);
      check |= (1 << 0);
      continue;
    }

    if (argv[i][1] == 's')
    {
      i++;
      K = argv[i][0] - '0';
      int j = 1;
      while (argv[i][j] != '\0')
      {
        K *= 10;
        K += argv[i][j] - '0';
        j++;
      }
      check |= (1 << 1);
      continue;
    }

    if (argv[i][1] == 'o')
    {
      i++;
      strcpy(ofilename, argv[i]);
      continue;
    }

    if (argv[i][1] == 'p')
    {
      int j;
      for (j = 0, i++; j < K; j++, i++)
        sampProb[j] = atof(argv[i]);
      samp = true;
      continue;
    }

    if (argv[i][1] == 'q')
    {
      itera = 1;
      break;
    }
  }

  if (largeScale)
    G = new DynamicGraph();
  else
    G = new GraphMatrix();

  if (!itera)
  {
    if (check != (1 << 2) - 1)
    {
      K = 0;
      if (check != 0)
        printf("Warning: Incorrect number of necessary arguments provided\n");
      displayHelp();
      return;
    }
    GraphUtils::readFileTxt(G, ifilename, dir, false, zeroBased);
    G->sortNeighbours();
    G->makeArrayNeighbours();
    if (ofilename[0] == '0' && ofilename[1] == '\0')
      outFile = stdout;
    else
      outFile = fopen(ofilename, "w");
    return;
  }

  // Direction
  printf("Directed? (Y/n) ");
  char chdir;
  scanf(" %c", &chdir);
  if (chdir == 'n' || chdir == 'N')
    dir = false;

  // Initial
  printf("Input 0 or 1 based: ");
  scanf("%d", &zeroBased);

  // Input filename
  printf("Insert input file name: ");
  scanf(" %s", ifilename);
  GraphUtils::readFileTxt(G, ifilename, dir, false, zeroBased);
  G->sortNeighbours();
  G->makeArrayNeighbours();

  // Subgraph Size
  printf("Input the value K of the subgraph search: ");
  scanf("%d", &K);

  // Input filename
  printf("Insert output file name or 0 to stdout: ");
  scanf(" %s", ofilename);
  if (ofilename[0] == '0' && ofilename[1] == '\0')
    outFile = stdout;
  else
    outFile = fopen(ofilename, "w");

  // Default Sampling probabilities
  if (!samp)
    for (i = 0; i < K; i++)
      sampProb[i] = 1.0;
}

void initSamplingProbabilities(Fase* fase)
{
  int i;
  prob = 1.0;
  for (i = 0; i < K; i++)
    prob *= sampProb[i];

  if (samp && fabs(prob) > 10e-7)
    fase->initSampling(K, sampProb);
  else
    prob = 1.0;
}

void output(Fase* fase)
{
  printf("Finished Calculating\n");
  FILE *f = outFile;
  fprintf(f, "\tOutput:\n");
  fprintf(f, "Network: %s\n", ifilename);
  fprintf(f, "Directed: %s\n", dir ? "Yes" : "No");
  fprintf(f, "Nodes: %d\n", G->numNodes());
  fprintf(f, "Edges: %d\n", G->numEdges() / (dir ? 1 : 2));
  fprintf(f, "Subgraph Size: %d\n", K);
  if (largeScale)
    fprintf(f, "Graph Representation: Large Scale\n");
  else
    fprintf(f, "Graph Representation: Adjacency Matrix\n");

  t_end = time(0);
  struct tm *tm_start = localtime(&t_start);
  fprintf(f, "Start of Computation: %02dh%02dm%02ds %02d/%02d/%02d\n\
", tm_start->tm_hour, tm_start->tm_min, tm_start->tm_sec, tm_start->tm_mday, tm_start->tm_mon + 1, 1900 + tm_start->tm_year);
  struct tm *tm_end   = localtime(&t_end);
  fprintf(f, "End of Computation: %02dh%02dm%02ds %02d/%02d/%02d\n", tm_end->tm_hour, tm_end->tm_min, tm_end->tm_sec, tm_end->tm_mday, tm_end->tm_mon + 1, 1900 + tm_end->tm_year);

  fprintf(f, "\n\n\tResults:\n");
  fprintf(f, "Subgraph Occurrences: %lld\n", (long long int)(fase->getMotifCount() / prob));
  fprintf(f, "Subgraph Types: %d\n", fase->getTypes());
  fprintf(f, "Computation Time (ms): %0.4lf\n", Timer::elapsed());

  if (fabs(prob - 1.0) <= 10e-7)
    fprintf(f, "\nExact Enumeration, no Sampling done\n");
  else
  {
    fprintf(f, "\n\tSampling Information:\n");
    fprintf(f, "Percentage of Sampled Subgraphs: %0.2lf\%\n", 100 * prob);
    fprintf(f, "Percentage by depth:\n");
    int i;
    for (i = 0; i < K; i++)
      fprintf(f, "P[%d]: %0.3lf\%\n", i, 100 * sampProb[i]);
  }

  if (detailed)
  {
    fprintf(f, "\n\tDetailed Output:\n");
    for (auto element : fase->subgraphCount())
      if (samp && fabs(prob) > 10e-7)
        fprintf(f, "%s: %d occurrences\n", element.second.c_str(), (int)(element.first / prob));
      else
        fprintf(f, "%s: %d occurrences\n", element.second.c_str(), element.first);
  }
}

void finish(Fase* fase)
{
  delete fase;
  delete G;
  fclose(outFile);
}

void myfinish(MyFase* fase)
{
  delete fase;
  delete G;
  fclose(outFile);
}

void genPermGraph(int *perm, Graph *orig, Graph *dest) {
  dest->createGraph(K, UNDIRECTED);
  for (int i=0; i!=K; ++i) {
    for (int j=0; j!=K; ++j) {
      if (orig->hasEdge(perm[i],perm[j])) {
        dest->addEdge(i,j);
      }
    }
  }
}

void printPermGraph(Graph *g) {
  for (int i=0; i!=K; ++i) {
    for (int j=0; j!=K; ++j) {
      printf("%d", g->hasEdge(i,j));
    }
  }
  putchar('\n');
}

void insertSubgraphPerm(int depth, int *perm, int *used, Graph *g, IGtrie *igtrie) {
  if (depth == K) {
    // generate cur permutation's graph
    Graph *_g = new GraphMatrix();
    genPermGraph(perm, g, _g);
    printPermGraph(_g);

    // insert permutation graph in the igtrie
    int _depth = 1;
    int nodeLabel = 0;
    long long int label = 0;
    int vsub[K]; for (int i=0; i!=K; ++i) vsub[i] = i;

    Label::init(_g, dir);

    while (_depth != K) {
      label = Label::updateLabel(vsub, vsub[_depth], _depth);
      nodeLabel = igtrie->insertLabel(nodeLabel, label, Label::repDigits(_depth));
      _depth++;
    }

    igtrie->incrementLabel(nodeLabel, 1);

    delete _g;
    return;
  }

  // recursion
  for (int i=0; i!=K; ++i) {
    if (used[i]) {
      continue;
    }

    used[i] = 1;
    perm[depth] = i;
    insertSubgraphPerm(depth+1, perm, used, g, igtrie);
    used[i] = 0;
  }
}

void insertSubgraph(Graph *g, IGtrie* igtrie) {
  int perm[K];
  int used[K] = {0};

  insertSubgraphPerm(0, perm, used, g, igtrie);
}

void readSubgraph(Graph *g) {
  char s[K*K+1];
  for (int i=0; i!=K; ++i) {
    scanf("%s", &s[i*K]);
  }
  GraphUtils::strToGraph(g, s, K, dir);
}

void simulateInsert(vector<int>& vsub, Graph* g, IGtrie *igtrie) {
  int *_vsub = &vsub[0];
  int depth = 1;
  int nodeLabel = 0;
  long long int label = 0;

  while (depth != K) {
    label = Label::updateLabel(_vsub, _vsub[depth], depth);
    nodeLabel = igtrie->insertLabel(nodeLabel, label, Label::repDigits(depth));
    depth++;
  }

  // igtrie->incrementLabel(nodeLabel, 1);
}

void expand(int depth, vector<int>& vsub, int *used, int **vext, int *vextSz,
            Graph *g, vector< vector<int> >& acc) {
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

    expand(depth+1, vsub, used, vext, vextSz, g, acc);
    used[curNode] = 0;
  }
}

void insertSubgraph_1(Graph *g, IGtrie *igtrie) {
  vector<int> vsub(K);
  vector< vector<int> > perm;
  int used[K] = {0};
  int **vext = NULL;
  int vextSz[K];

  Label::init(g, dir);

  vext = new int*[K];
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

    expand(1, vsub, used, vext, vextSz, g, perm);
    used[i] = 0;
  }

  for (auto p: perm) {
    /*
    for (int v: p) {
      cout << v << " ";
    }
    cout << "\n";
    */
    simulateInsert(p, g, igtrie);
  }

  for (int i=1; i!=K; ++i) {
    delete vext[i];
  }
  delete vext;
}

void iter1() {
  MyFase *fase = new MyFase(G, dir, K);

  int ni;
  scanf("%d", &ni);
  while (ni--) {
    Graph *g = new GraphMatrix();
    readSubgraph(g);
    fase->setQuery(g);
    // insertSubgraph(g, igtrie);
    //insertSubgraph_1(g, igtrie);
    delete g;
  }

  // fase->igtrieEnum();
  fase->runCensus();

  if (detailed) {
    cout << "\n\tDetailed Output:\n";
    for (auto element : fase->subgraphCount())
      cout << element.second << ": " << element.first << "\n";
  }

  myfinish(fase);
}

int main(int argc, char **argv)
{
  init();
  read(argc, argv);

  if (K <= 2 || K >= MAXMOTIF)
  {
    fprintf(stderr, "Subgraph size needs to be between 3 and %d...\n", MAXMOTIF - 1);
    return 1;
  }

  iter1();
  /*
  Random::init(time(NULL));
  Fase* fase = new Fase(G, dir);
  initSamplingProbabilities(fase);

  Timer::start();
  fase->runCensus(K);
  Timer::stop();

  output(fase);
  finish(fase);
  */

  return 0;
}
