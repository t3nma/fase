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
#include "Graphlets5.h"

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
long long int Fase::count[50];
long long int **Fase::type;
int Fase::clique[6] = {0, 0, 1, 3, 31, 511};
int Fase::cliqueCount;

//int types[MAXGRAPHS];
int Fase::graphlets[30];
long long int **Fase::orbits;


int mapLabelToGraph[32] = {-1, 12, 5, 6,-1,-1,-1,-1,-1, 1,-1,
                            2, 7, 4,-1, 3,-1,15,19,16,21,
                           18,20,17,-1, 8,-1, 9,14,11,-1,
                           10};

int mapGraphToGraphlet[22] = {-1, 4, 6, 7, 6, 1, 2, 3, 6, 7, 8,
                           7, 1,-1, 6, 3, 6, 7, 5, 4, 6,
                           3};

int mapGraphToOrbits[22][4] = {{-1,-1,-1,-1},{7,6,6,6}, {11,10,9,10}, {13,12,12,13}, {11,9,10,10}, {1, 2, 1, -1}, {3, 3, 3, -1}, {5,4,5,4},
                                  {11,10,10,9}, {13,13,12,12}, {14,14,14,14}, {13,12,13,12}, {2, 1, 1, -1}, {-1,-1,-1,-1}, {10,10,11,9},
                                  {5,5,4,4}, {10,11,9,10}, {12,13,12,13}, {8,8,8,8}, {6,7,6,6},
                                  {11,9,10,10}, {4,5,5,4}};

FILE* typeFile;

/****** ORCA ******/
#include <unordered_map>

long long **orbit;

int *tri;

typedef long long int64;
typedef pair<int,int> PII;

struct PAIR {
    int a, b;
    PAIR(int a0, int b0) { a=min(a0,b0); b=max(a0,b0); }
};

bool operator<(const PAIR &x, const PAIR &y) {
    if (x.a==y.a) return x.b<y.b;
    else return x.a<y.a;
}

bool operator==(const PAIR &x, const PAIR &y) {
    return x.a==y.a && x.b==y.b;
}

struct hash_PAIR {
    size_t operator()(const PAIR &x) const {
        return (x.a<<8) ^ (x.b<<0);
    }
};

struct TRIPLE {
    int a, b, c;
    TRIPLE(int a0, int b0, int c0) {
        a=a0; b=b0; c=c0;
        if (a>b) swap(a,b);
        if (b>c) swap(b,c);
        if (a>b) swap(a,b);
    }
};

bool operator<(const TRIPLE &x, const TRIPLE &y) {
    if (x.a==y.a) {
        if (x.b==y.b) return x.c<y.c;
        else return x.b<y.b;
    } else return x.a<y.a;
}
bool operator==(const TRIPLE &x, const TRIPLE &y) {
    return x.a==y.a && x.b==y.b && x.c==y.c;
}

struct hash_TRIPLE {
    size_t operator()(const TRIPLE &x) const {
        return (x.a<<16) ^ (x.b<<8) ^ (x.c<<0);
    }
};

//unordered_map<PAIR, int, hash_PAIR> common2;
unordered_map<TRIPLE, int, hash_TRIPLE> common3;
//unordered_map<PAIR, int, hash_PAIR>::iterator common2_it;
unordered_map<TRIPLE, int, hash_TRIPLE>::iterator common3_it;

#define common3_get(x) (((common3_it=common3.find(x))!=common3.end())?(common3_it->second):0)
//#define common2_get(x) (((common2_it=common2.find(x))!=common2.end())?(common2_it->second):0)

int **common_x ;
PAIR *edges;
int **inc;

int **common2_adjm;

/****** ORCA ******/

void Fase::destroy()
{
  int i;
  delete[] sub;
  for (i = 0; i < K; i++)
    delete[] extCpy[i];
  delete [] extCpy;
}

void Fase::GraphletsCount(Graph *_G, int _K)
{
  K = _K;
  G = _G;
  sub = new int[K];
  graphSize = G->numNodes();
  cliqueCount = 0;
  memset(count, 0, sizeof count);
  //memset(type, 0, sizeof type);
  int i, j, extNum = 0;
  extCpy = new short*[K];

  for (i = 0; i < K; i++)
    extCpy[i] = new short[graphSize];


  type = (long long int**)malloc(G->numNodes()*sizeof(long long int*));
  for (int i=0; i < G->numNodes(); i++) type[i] = (long long int*)malloc(MAXGRAPHS*sizeof(long long int));

  orbits = (long long int**)malloc(G->numNodes()*sizeof(long long int*));
  for (int i=0; i < G->numNodes(); i++) orbits[i] = (long long int*)malloc(73*sizeof(long long int));

  /*
  memset(types, 0, sizeof types);
  typeFile = fopen("types.txt", "w");
  */

  buildCommonNodes();

  for (i = 0; i < graphSize; i++)
  {
    sub[0]     = i;
    int *nei   = G->arrayNeighbours(i);
    int neiNum = G->numNeighbours(i);

    orbits[i][0] = neiNum;

    extNum = 0;
    for (j = 0; j < neiNum; j++)
      if (nei[j] > i)
        extCpy[0][extNum++] = nei[j];
    GraphletsExtendSubgraph(extNum, 0);

    solveEquations(i);
    calcOrbitFrequency(i);
  }

  //solveEquations();
  //calcOrbitFrequency();
  //calcGraphletFrequency();

  //fclose(typeFile);
}

void Fase::GraphletsExtendSubgraph(int extNum, int node)
{
  int graph, graphlet;
  int x, a, b, c, xa, xb, xc, ab, ac, bc;

  int ncases;
  int* mycase;

  bool** adjM    = G->adjacencyMatrix();
  int * deg      = G->arrayNumNeighbours();

  int node1   = node;
  int extNum1 = extNum;

  int i1, j1, o1, i2, j2, o2, i3;
  int extCpyNum1, extCpyNum2;
  int exti1, *eExcl1, eExclNum1, exti2, *eExcl2, eExclNum2, exti3, *eExcl3, eExclNum3;

  memcpy(extCpy[1], extCpy[0], extNum1 * sizeof(short));

  // K = 2
  for (i1 = extNum1 - 1; i1 >= 0; i1--)
  {
    extCpyNum1 = i1;
    exti1      = extCpy[0][i1];
    eExcl1     = G->arrayNeighbours(exti1);
    eExclNum1  = G->numNeighbours(exti1);

    for (j1 = 0; j1 < eExclNum1; j1++)
    {
      int eEj1 = eExcl1[j1];
      if (eEj1 <= sub[0]) continue;
      for (o1 = 0; o1 < 1; o1++)
        if (eEj1 == sub[o1] || G->isConnected(eEj1, sub[o1])) break;
      if (o1 == 1)
        extCpy[1][extCpyNum1++] = eEj1;
    }
    sub[1] = exti1;

    // K = 3
    int node2   = (node1 << 1);
    int extNum2 = extCpyNum1;

    memcpy(extCpy[2], extCpy[1], extNum2 * sizeof(short));

    for (i2 = extNum2 - 1; i2 >= 0; i2--)
    {
      extCpyNum2 = i2;
      exti2      = extCpy[1][i2];
      eExcl2     = G->arrayNeighbours(exti2);
      eExclNum2  = G->numNeighbours(exti2);

      for (j2 = 0; j2 < eExclNum2; j2++)
      {
        int eEj2 = eExcl2[j2];
        if (eEj2 <= sub[0])
          continue;
        for (o2 = 0; o2 < 2; o2++)
          if (eEj2 == sub[o2] || G->isConnected(eEj2, sub[o2]))
            break;
        if (o2 == 2)
          extCpy[2][extCpyNum2++] = eEj2;
      }

      int nm = 0;

      bool *p = adjM[exti2];
      for (int j2 = 0; j2 < 2; j2++)
        nm |= ((int)(*(p + sub[j2])) << j2);

      int myType = (node2 << 2) | nm;
      graph    = mapLabelToGraph[myType];
      graphlet = mapGraphToGraphlet[graph];

      int* orbitsandstuff = mapGraphToOrbits[graph];

      sub[2] = exti2;

      orbits[sub[0]][orbitsandstuff[0]]++;
      orbits[sub[1]][orbitsandstuff[1]]++;
      orbits[sub[2]][orbitsandstuff[2]]++;


      // K = 4
      int node3   = myType;
      int extNum3 = extCpyNum2;

      memcpy(extCpy[3], extCpy[2], extNum3 * sizeof(short));

      for (i3 = 0; i3 < extNum3; i3++)
      {
        exti3  = extCpy[2][i3];
        sub[3] = exti3;


        int nm = 0;
        bool *p = adjM[exti3];
        //for (int j3 = 0; j3 < 3; j3++)
          nm |= ((int)(*(p + sub[0]))) | ((int)(*(p + sub[1])) << 1) | ((int)(*(p + sub[2])) << 2) ;

        int myType = (node3 << 3) | nm;

        //MotifCount++;

        graph    = mapLabelToGraph[myType];
        graphlet = mapGraphToGraphlet[graph];

        int* orbitsandstuff = mapGraphToOrbits[graph];

        orbits[sub[0]][orbitsandstuff[0]]++;
        orbits[sub[1]][orbitsandstuff[1]]++;
        orbits[sub[2]][orbitsandstuff[2]]++;
        orbits[sub[3]][orbitsandstuff[3]]++;

        /** Graphlet 3 **/
        if(graphlet == 3){
            // Orbit 4
            ncases = 2; //Graphlets5::getNCases(graph, 4);

            for(int i = 0; i < ncases; i++){
              mycase      = Graphlets5::getCase(graph, 4, i);

              x  = sub[mycase[0]]; a  = sub[mycase[1]]; b  = sub[mycase[2]]; c  = sub[mycase[3]];

              bc = inc[b][c];

              orbits[x][35] += common_x[a][c] - 1;
              orbits[x][34] += common_x[x][c];
              orbits[x][27] += tri[bc];
              orbits[x][18] += deg[b] - 2;
              orbits[x][15] += deg[c] - 1;
            }
        }

        /** Graphlet 4 **/
        else if(graphlet == 4){
            // Orbit 6

            // Orbit 7
            mycase      = Graphlets5::getCase(graph, 7, 0);

            x  = sub[mycase[0]];

            orbits[x][23] += deg[x] - 3;
        }

        /** Graphlet 5 **/
        else if(graphlet == 5){
            // Orbit 8
            ncases = 4; //Graphlets5::getNCases(graph, 8);

            for(int i = 0; i < ncases; i++){
              mycase      = Graphlets5::getCase(graph, 8, i);

              x  = sub[mycase[0]]; a  = sub[mycase[1]]; b  = sub[mycase[2]]; c  = sub[mycase[3]];

              xa = inc[x][a]; xb = inc[x][b];

              orbits[x][53] += tri[xa] + tri[xb];
              orbits[x][50] += common_x[x][c] - 2;
            }
        }

        /** Graphlet 6 **/
        else if(graphlet == 6){
            // Orbit 11
            mycase      = Graphlets5::getCase(graph, 11, 0);

            x  = sub[mycase[0]]; c  = sub[mycase[3]];
            xc = inc[x][c];

            orbits[x][44] += tri[xc];

            // Orbit 10

            // Orbit 9
            mycase      = Graphlets5::getCase(graph, 9, 0);

            a  = sub[mycase[1]]; b  = sub[mycase[2]]; c  = sub[mycase[3]];
            ab = inc[a][b]; ac = inc[a][c];

            orbits[x][45] += common2_adjm[b][c] - 1;
            orbits[x][39] += tri[ab] + tri[ac] - 2;
            orbits[x][31] += deg[a] - 3;
            orbits[x][24] += deg[b] + deg[c] - 4;
        }

        /** Graphlet 7 **/
        else if(graphlet == 7){
            // Orbit 13
            ncases = Graphlets5::getNCases(graph, 13);

            for(int i = 0; i < ncases; i++){
              mycase      = Graphlets5::getCase(graph, 13, i);

              x  = sub[mycase[0]]; a  = sub[mycase[1]]; b  = sub[mycase[2]]; c  = sub[mycase[3]];
              xa = inc[x][a]; xb = inc[x][b]; xc = inc[x][c];

              orbits[x][68] += common3_get(TRIPLE(a,b,c)) - 1;
              orbits[x][64] += common2_adjm[b][c] - 2;
              orbits[x][61] += tri[xb] + tri[xc] - 2;
              orbits[x][55] += tri[xa] - 2;
            }

            // Orbit 12
        }

        /** Graphlet 8 **/
        else if (graphlet == 8)
        {
          //Calculo Clique-K
          //if(GTrieNode::adjM[a][neigh] && GTrieNode::adjM[b][neigh] && GTrieNode::adjM[c][neigh]) GTrie::orbit_freq[72]++;
          eExcl3    = G->arrayNeighbours(exti3);
          eExclNum3 = G->numNeighbours(exti3);
          for (int j3 = 0; j3 < eExclNum3; j3++)
          {
            int eEj3 = eExcl3[j3];

            bool *p = adjM[eEj3];
            int fl = 1;
            for (int l = 0; fl && l < 3; l++)
              fl &= ((sub[l] >= eEj3) & (int)*(p + sub[l]));

            orbits[x][72] += fl * 5;
           }
           // Orbit 14

           ncases = 4; //Graphlets5::getNCases(graph, 14);

           for(int i = 0; i < ncases; i++){
             mycase      = Graphlets5::getCase(graph, 14, i);

             x  = sub[mycase[0]]; a  = sub[mycase[1]]; b  = sub[mycase[2]]; c  = sub[mycase[3]];

             xa = inc[x][a]; xb = inc[x][b]; xc = inc[x][c];

             orbits[x][70] += common3_get(TRIPLE(a, b, c)) - 1;
             orbits[x][67] += tri[xa] + tri[xb] + tri[xc] - 6;
             orbits[x][58] += deg[x] - 3;
           }
        }
       }
     }
   }
}

void Fase::buildCommonNodes(){
    bool** adjM    = G->adjacencyMatrix();
    int ** fastnei = G->matrixNeighbours();

    common2_adjm = (int**)malloc(G->numNodes()*sizeof(int*));
    for (int i=0; i < G->numNodes(); i++)
        common2_adjm[i] = (int*)malloc(G->numNodes()*sizeof(int));

    printf("computing common nodes\n");

    int m = G->numEdges();
    edges = (PAIR*)malloc(m/2*sizeof(PAIR));

    inc = (int**)malloc(G->numNodes()*sizeof(int*));
    for (int i=0; i < G->numNodes(); i++) inc[i] = (int*)malloc(G->numNodes()*sizeof(int)); //no need, should be g->numNeighbours(x)

    common_x = (int**)malloc(G->numNodes()*sizeof(int*));
    for (int i=0; i < G->numNodes(); i++) common_x[i] = (int*)malloc(G->numNodes()*sizeof(int));

    int edge = 0;

    for (int x = 0; x < graphSize; x++) {
        for (int n1 = 0; n1 < G->numNeighbours(x); n1++) {
            int a = fastnei[x][n1];
            if(x < a) {
                edges[edge] = PAIR(x, a);
                inc[x][a]=edge; inc[a][x]=edge;
                edge++;
            }
            for (int n2 = n1 + 1; n2 < G->numNeighbours(x); n2++) {
                int b   = fastnei[x][n2];
                //PAIR ab = PAIR(a,b);
                //common2[ab]++;
                common2_adjm[a][b]++;
                common2_adjm[b][a]++;
                for (int n3 = n2 + 1; n3 < G->numNeighbours(x); n3++) {
                    int c  = fastnei[x][n3];
                    int st = adjM[a][b] + adjM[a][c] + adjM[b][c];
                    if (st < 2) continue;
                        TRIPLE abc = TRIPLE(a,b,c);
                        /*TRIPLE acb = TRIPLE(a,c,b);
                        TRIPLE bac = TRIPLE(b,a,c);
                        TRIPLE bca = TRIPLE(b,c,a);
                        TRIPLE cab = TRIPLE(c,a,b);
                        TRIPLE cba = TRIPLE(c,b,a);*/
                        common3[abc]++;
                        //common3[cba]=common3[cab]=common3[bca]=common3[bac]=common3[acb] = common3[abc];
                }
            }
            for (int na = 0; na < G->numNeighbours(a);na++) {
                int b = fastnei[a][na];
                if (b!=x  && !(adjM[x][b])) {
                    common_x[x][b]++;
                }
            }
        }
   }

    m   = m/2;
    tri = (int*)calloc(m,sizeof(int));
    for (int i = 0; i < m; i++) {
        int x=edges[i].a, y=edges[i].b;
        for (int xi = 0, yi = 0; xi < G->numNeighbours(x) && yi < G->numNeighbours(y); ) {
            if (fastnei[x][xi] == fastnei[y][yi]) { tri[i]++; xi++; yi++; }
            else if (fastnei[x][xi] < fastnei[y][yi]) { xi++; }
            else { yi++; }
        }
     }
}

void Fase::solveEquations(int x){
    orbits[x][70] = (orbits[x][70]-4*orbits[x][72]);
        orbits[x][68] = (orbits[x][68]-3*orbits[x][70]);
        orbits[x][67] = (orbits[x][67]-12*orbits[x][72]-6*orbits[x][70]);
        orbits[x][64] = (orbits[x][64]-3*orbits[x][70]-1*orbits[x][68]-1*orbits[x][68]);
        orbits[x][61] = (orbits[x][61]-6*orbits[x][70]-2*orbits[x][68]-2*orbits[x][67])/2;
        orbits[x][58] = (orbits[x][58]-4*orbits[x][72]-3*orbits[x][70]-1*orbits[x][67]);
        orbits[x][55] = (orbits[x][55]-3*orbits[x][70]-2*orbits[x][67])/3;
        orbits[x][53] = (orbits[x][53]-2*orbits[x][68]-2*orbits[x][64]-2*orbits[x][64]);
        orbits[x][50] = (orbits[x][50]-1*orbits[x][68]-2*orbits[x][64])/3;
        orbits[x][45] = (orbits[x][45]-1*orbits[x][67]-1*orbits[x][64]-3*orbits[x][58]);
        orbits[x][44] = (orbits[x][44]-1*orbits[x][67]-2*orbits[x][61])/4;
        orbits[x][39] = (orbits[x][39]-2*orbits[x][67]-2*orbits[x][61]-6*orbits[x][58])/2;
        orbits[x][35] = (orbits[x][35]-2*orbits[x][61]-1*orbits[x][53]-2*orbits[x][45])/2;
        orbits[x][34] = (orbits[x][34]-2*orbits[x][61]-2*orbits[x][53])/2;
        orbits[x][31] = (orbits[x][31]-1*orbits[x][67]-2*orbits[x][61]-3*orbits[x][58]-4*orbits[x][44]-2*orbits[x][39]);
        orbits[x][27] = (orbits[x][27]-2*orbits[x][61]-1*orbits[x][53]-2*orbits[x][45])/2;
        orbits[x][24] = (orbits[x][24]-2*orbits[x][67]-2*orbits[x][64]-2*orbits[x][61]-6*orbits[x][58]-1*orbits[x][53]-2*orbits[x][45]-2*orbits[x][39]);
        orbits[x][23] = (orbits[x][23]-1*orbits[x][55]-1*orbits[x][39]-1*orbits[x][31])/4;
        orbits[x][18] = (orbits[x][18]-2*orbits[x][61]-1*orbits[x][53]-4*orbits[x][45]-2*orbits[x][35]-2*orbits[x][27]-1*orbits[x][24])/2;
        orbits[x][15] = (orbits[x][15]-2*orbits[x][61]-2*orbits[x][53]-2*orbits[x][45]-2*orbits[x][35]-2*orbits[x][34]-2*orbits[x][27]);
}



void Fase::calcOrbitFrequency(int i){
    int myType, myGraph, freq, k;

    for(myType = 0; myType < 32; myType++){
        myGraph = mapLabelToGraph[myType];
        if(myGraph == -1) continue; //should remove -1s

        freq = type[i][myType];

        graphlets[mapGraphToGraphlet[myGraph]] += freq;
        //for(k = 0; k < K; k++)
          //orbits[i][mapGraphToOrbits[myGraph][k]] += freq;
    }

    // G9
    orbits[i][16] = orbits[i][15];
        orbits[i][17] = orbits[i][15]/2;

        // G10
        orbits[i][19] = orbits[i][18] * 2;
        orbits[i][20] = orbits[i][18];
        orbits[i][21] = orbits[i][18];

        // G11
        orbits[i][22] = orbits[i][23] * 4;

        // G12
        orbits[i][25] = orbits[i][24]/2;
        orbits[i][26] = orbits[i][24];

        // G13
        orbits[i][28] = orbits[i][27];
        orbits[i][29] = orbits[i][27] * 2;
        orbits[i][30] = orbits[i][27];

        // G14
        orbits[i][32] = orbits[i][31];
        orbits[i][33] = orbits[i][31]/2;

        // G15

        // G16
        orbits[i][36] = orbits[i][35];
        orbits[i][37] = orbits[i][35] * 2;
        orbits[i][38] = orbits[i][35];

        // G17
        orbits[i][40] = orbits[i][39] * 2;
        orbits[i][41] = orbits[i][39];
        orbits[i][42] = orbits[i][39];

        // G18
        orbits[i][43] = orbits[i][44] * 4;

        // G19
        orbits[i][46] = orbits[i][45];
        orbits[i][47] = orbits[i][45];
        orbits[i][48] = orbits[i][45] * 2;

        // G20
        orbits[i][49] = 3*orbits[i][50]/2;

        // G21
        orbits[i][51] = orbits[i][53];
        orbits[i][52] = orbits[i][53]/2;

        // G22
        orbits[i][54] = 3*orbits[i][55]/2;

        // G23
        orbits[i][56] = orbits[i][58];
        orbits[i][57] = orbits[i][58] * 3;

        // G24
        orbits[i][59] = orbits[i][61] * 2;
        orbits[i][60] = orbits[i][61] * 2;

        // G25
        orbits[i][62] = orbits[i][64]/2;
        orbits[i][63] = orbits[i][64];

        // G26
        orbits[i][65] = orbits[i][67]/2;
        orbits[i][66] = orbits[i][67];

        // G27
        orbits[i][69] = orbits[i][68]/4;

        // G28
        orbits[i][71] = 3*orbits[i][70]/2;
}

void Fase::calcGraphletFrequency(){
    //graphlets[0]+= orbits[0]/2;
}
