/* -------------------------------------------------

//                                                 
//  88888888888           ad88888ba   88888888888  
//  88                   d8"     "8b  88           
//  88                   Y8,          88           
//  88aaaaa  ,adPPYYba,  `Y8aaaaa,    88aaaaa      
//  88"""""  ""     `Y8    `"""""8b,  88"""""      
//  88       ,adPPPPP88          `8b  88           
//  88       88,    ,88  Y8a     a8P  88           
//  88       `"8bbdP"Y8   "Y88888P"   88888888888  
//                                                 
//

Pedro {Paredes, Ribeiro} - DCC/FCUP

----------------------------------------------------
Base FaSE implementation

---------------------------------------------------- */

#include "Fase.h"

int Fase::K;
long long int Fase::MotifCount = 0;
int Fase::typeLabel;
bool Fase::directed;
Graph* Fase::G;
int *Fase::sub;
int Fase::subNum;
int Fase::graphSize;
short **Fase::extCpy;
char Fase::globStr[MAXS];
char Fase::s[20 * 20 + 1];

long long int Fase::getTypes()
{
  return GTrie::getCanonicalNumber();
}

long long int Fase::getLeafs()
{
  return GTrie::getClassNumber();
}

long long int Fase::getNodes()
{
  return GTrie::getNodeNumber();
}

void Fase::listTree(FILE* f)
{
  return GTrie::listGtrie(f);
}

void Fase::listClasses(FILE* f)
{
  return GTrie::listClasses(f);
}

/*! Generate the next piece of the LSLabel
    \param w the newly added vertex*/
char* Fase::LSLabel(int w, int subSize)
{
  return LSLabeling::Label(sub, subSize, w, typeLabel, globStr, directed);
}

void Fase::destroy()
{
  int i;
  GTrie::destroy();
  delete[] sub;
  for (i = 0; i < K; i++)
    delete[] extCpy[i];
  delete [] extCpy;
}

/*! List all k-subgraphs in a larger graph
    \param G the graph to be explored
    \param K the size of the subgraphs*/
void Fase::EnumerateSubgraphs(Graph *_G, int _K)
{
  K = _K;
  G = _G;
  sub = new int[K];
  graphSize = G->numNodes();
  GTrie::init();
  LSLabeling::init(G);
  int i, j, extNum = 0;
  extCpy = new short*[K];

  for (i = 0; i < K; i++)
    extCpy[i] = new short[graphSize];

  for (i = 0; i < graphSize; i++)
  {
    sub[0] = i;
    int *nei = G->arrayNeighbours(i);
    int neiNum = G->numNeighbours(i);
    extNum = 0;
    for (j = 0; j < neiNum; j++)
      if (nei[j] > i)
        extCpy[0][extNum++] = nei[j];
    ExtendSubgraph(extNum, 1);
  }
  //  printf("\t\tWith G-Tries:\n");
  //   printf("Found %d Motifs\n", MotifCount);
  //  printf("LS-Classes: %lld\n", GTrie::getClassNumber());
  //  GTrie::listClasses();
  //  GTrie::listGtrie();
}

/*! Generates the next node on the Fase-tree
    \param ext the list of possible next chosen vertices
    \param extNum the number of elements in ext
    \param subSize the size of the enumerating subgraph
    \param v the label of the base vertex*/
void Fase::ExtendSubgraph(int extNum, int subSize)
{
  int i;

  if (subSize == K - 1)
  {
    for (i = 0; i < extNum; i++)
    {
      int exti = extCpy[subSize - 1][i];
      GTrie::insert(LSLabel(exti, subSize));
      sub[subSize] = exti;
      MotifCount++;

      if (!GTrie::getCurrentLeaf())
      {
        s[0] = '\0';
        Isomorphism::canonicalStrNauty(G, sub, s);
        GTrie::setCanonicalLabel(s);
      }

      GTrie::jump();    
    }
    return;
  }

  int j, o;
  int extCpyNum;
  int exti, *eExcl, eExclNum;
  
  memcpy(extCpy[subSize], extCpy[subSize - 1], extNum * sizeof(short));

  for (i = extNum - 1; i >= 0; i--)
  {
    extCpyNum = i;
    exti = extCpy[subSize - 1][i];
    eExcl = G->arrayNeighbours(exti);
    eExclNum = G->numNeighbours(exti);

    for (j = 0; j < eExclNum; j++)
    {
      int eEj = eExcl[j];
      if (eEj <= sub[0])
        continue;
      for (o = 0; o < subSize; o++)
        if (eEj == sub[o] || G->isConnected(eEj, sub[o]))
          break;
      if (o == subSize)
        extCpy[subSize][extCpyNum++] = eEj;
    }

    GTrie::insert(LSLabel(exti, subSize));
    sub[subSize] = exti;

    ExtendSubgraph(extCpyNum, subSize + 1);
 
    GTrie::jump();
  }
}
