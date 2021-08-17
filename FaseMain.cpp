#include "Common.h"
#include "Fase.h"
#include "DynamicGraph.h"
#include "GraphMatrix.h"
#include "GraphUtils.h"
#include "Timer.h"
#include "Graph.h"
#include "Random.h"

using namespace std;

Graph *G;
int K = 0, zeroBased = 1;
double sampProb[MAXMOTIF], prob;
bool dir = false, detailed = false, draw = false, samp = false, largeScale = false, monitor = false, monitor2 = false;
char ifilename [200];
char ufilename [200];
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
  printf("------------ FaSE Usage --------------\nMain Settings: ./FASE -s <Subgraph Size> -i <input file> -u <stream file> [arguments...]\n\n\tAll commands:\n-h : Displays this help information\n-s <Integer> : Subgraph Size\n-i <Filename> : Name of input file (Format in Readme.txt)\n-u <Filename> : Name of stream file (Format in Readme.txt)\n-d : Directed Subgraph (Default undirected)\n-o : Name of output file (Default is stdout)\n-dt : Detailed Result (Displays all subgraph types and occurrences)\n-ls : Use a large scale representation (default is adjacency matrix)\n-z : Use 0-based input (Suitable for input files starting at node 0)\n-m : Monitor mode (reports new occurrences for each update)\n-p <P1> <P2> ... <Ps> : Sets the sampling probabilities by depth (note that -s must have been selected first)\n-q : Ignore arguments and prompt input\n--------------------------------------\n");
}

void read(int argc, char **argv)
{
  int E, V, i, check = 0, itera = 0;
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

    if (argv[i][1] == 'm')
    {
      if (argv[i][2] == '2')
        monitor2 = true;
      else
        monitor = true;
    }

    if (argv[i][1] == 'i')
    {
      i++;
      strcpy(ifilename, argv[i]);
      check |= (1 << 0);
      continue;
    }

    if (argv[i][1] == 'u')
    {
      i++;
      strcpy(ufilename, argv[i]);
      check |= (1 << 1);
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
      check |= (1 << 2);
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
    if (check != (1 << 3) - 1)
    {
      K = 0;
      if (check != 0)
        printf("Warning: Incorrect number of necessary arguments provided\n");
      displayHelp();
      return;
    }
    GraphUtils::readFileTxt(G, ifilename, dir, false, zeroBased);
    G->sortNeighbours();
    // G->makeArrayNeighbours();
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
  // G->makeArrayNeighbours();

  // Input filename
  printf("Insert stream file name: ");
  scanf(" %s", ufilename);

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

  // Monitor mode
  printf("Monitor mode? (y/N) ");
  char chmonitor;
  scanf(" %c", &chmonitor);
  if (chmonitor == 'y' || chmonitor == 'Y')
    monitor = true;

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
  fprintf(f, "Stream: %s\n", ufilename);
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
  fprintf(f, "Computation Time (ms): %0.4lf\n", Timer::elapsed() * 1000);

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
    for (auto element : fase->subgraphCount(monitor || monitor2))
      if (samp && fabs(prob) > 10e-7)
        fprintf(f, "%s: %d occurrences\n", element.first.c_str(), (int)(element.second / prob));
      else
        fprintf(f, "%s: %d occurrences\n", element.first.c_str(), element.second);
  }
}

void outputOccur(Fase *fase, int u = -1, int v = -1, bool increment = true)
{
  FILE *f = outFile;
  bool isMonitor = monitor || monitor2;
  auto counters = fase->subgraphCount(isMonitor);
  int motifCount = fase->getMotifCount();

  if (u == -1)
    fprintf(f, "census: ");
  else
    fprintf(f, "%c(%d,%d): ", "-+"[increment], u, v);

  fprintf(f, "%d occurrences\n", motifCount);

  for (auto elem : counters)
    fprintf(f, "%s: %d\n", elem.first.c_str(), elem.second);

  fprintf(f, "\n");
}

void finish(Fase* fase)
{
  delete fase;
  delete G;
  fclose(outFile);
}

bool readSubgraph(Graph *g)
{
  string s;
  cin >> s;

  int size = (int)sqrt(s.size() / 1.0);

  if (size <= 2 || size > K) {
    return false;
  }

  GraphUtils::strToGraph(g, s.c_str(), size, dir);
  return true;
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

  Fase* fase = new Fase(G, dir, K);
  initSamplingProbabilities(fase);

  // Subgraph input
  // TODO: input from file inside read() ?
  int ni;
  Graph *g = NULL;
  scanf("%d", &ni);
  while (ni--) {
    g = new GraphMatrix();
    if (readSubgraph(g)) {
      if (monitor2)
        fase->setQuery2(g);
      else
        fase->setQuery(g);

      delete g;
      g = NULL;
    }
  }

  delete g;
  g = NULL;

  /*
  Timer::start();
  fase->runCensus();
  Timer::stop();
  output(fase);
  */

  if (monitor || monitor2) {
    fase->setupMonitor();
  } else {
    fase->runCensus();
    outputOccur(fase);
  }

  FILE *f = fopen(ufilename, "r");
  if (!f)
    exit(EXIT_FAILURE);

  char op;
  int u, v, V = G->numNodes();
  bool inc;

  while (fscanf(f, "%c %d %d\n", &op, &u, &v) == 3) {
    u -= zeroBased;
    v -= zeroBased;

    if (op != 'A' && op != 'R')
    {
      cout << "Unknown stream command '" << op << "'\n";
      continue;
    }
    else
      inc = (op == 'A');

    if (u >= V || v >= V)
    {
      cout << "Edge (" << u << "," << v << ") exceeds graph size\n";
      continue;
    }

    if ( u == v                    ||
         (inc && G->hasEdge(u, v)) ||
         (!inc && !G->hasEdge(u, v)) )
    {
      cout << "Irrelevant (" << u << "," << v << ") update\n";
      continue;
    }

    if (monitor)
      fase->monitor(u, v, inc);
    else if (monitor2)
      fase->monitor2(u, v, inc);
    else
      fase->updateCensus(u, v, inc);

    outputOccur(fase, u, v, inc);
  }

  fclose(f);
  finish(fase);

  return 0;
}
