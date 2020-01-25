/*
  Rural Postman Probem
  Dado um Grafo simples G = (V, E) e R, um subconjunto de E.
  Encontrar um caminho fechado de custo mínimo que passe por todas as arestas de R.
  S(X) = conjunto de arestas que possui uma ponta em um vértice do conjunto X (subconjunto de V) e outra em V\X
  
  Minimizar: somatorio(xe * cost(e)), "e" pertence ao conjunto da arestas E. cost(e) = peso da aresta "e".

  Restrições
  (1): somatorio(x[e]) - 2*z[i] = 0 (garante que o caminho seja fechado, grau par)
  (2): Seja W = conjunto formado por todos os subconjuntos de Vr = {V1, V2, ..., Vp};
        Então, para todo w e W.
        somatorio(xe) >= 2, para todo xe que pertence a S(w)
        Garante que o grafo seja conexo

  Resolvendo modelo utilizando Gurobi.
*/
#include "gurobi_c++.h"
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <sstream>
#include <stack>
#include <cstring>
#include <set>
#include <vector>
#include <algorithm>
#include "dsu.h"
#define pb push_back
#define f first
#define s second
using namespace std;
#define N 10005
#define INF 200000000
#define eps 0.01
typedef pair<int, int> pii;
string itos(int i) {stringstream s; s << i; return s.str(); }

int n, m, req, nreq, cost[N], C[N], especial[N], V_especiais;

struct edge
{
  int v, cost, req, id, x;
  edge(int a = 0, int x = 0, int b = 0, int c = 0, int fr = 0)
  {
    v = a;
    cost = x;
    req = b;
    id = c;
    x = fr;
  }
};

vector<edge> grafo[N], G[N];
vector<edge> edges;

/*
  Dada uma possível solução válida (de acordo com a restrição (1)), encontrar a componente
  ue possui tamanho mínimo. A componente é inserida no vector subset.

  Complexidade O(N * alpha(N))
*/
void procurar_componentes(int n, int m, double *x, vector<int> &subset)
{
  DSU *UF = create_DSU(n);

   build(UF, n);
  for(int i = 0; i < m; i++)
  {
    if((x[i] >= 1.0 - eps) or edges[i].req == 1)
     join(UF, edges[i].x, edges[i].v);
  }

  int menor = INF, opt = -1;

  for(int i = 0; i < n; i++)
    if(especial[i])
      menor = min(menor, get_sz(UF, i));

  for(int i = 0; i < n; i++)
    if(especial[i] and get_sz(UF, i) == menor)
    {
      opt = Find(UF, i);
      break;
    }

  for(int i = 0; i < n; i++)
    if(especial[i] and Find(UF, i) == opt) subset.push_back(i);

  clear(UF);
}

/*
  Callback. Sempre que uma possível solução é encontrada, pego sua componente de tamanho mínimo e 
  adiciona a restrição (2) para essa componente em relação aos outros vértices.

  Complexidade O(N*log(N))
*/
class subtourelim: public GRBCallback
{
  public:
    GRBVar* vars;
    int n, m;
    subtourelim(GRBVar* xvars, int xn, int xm)
    {
      vars = xvars;
      n  = xn;
      m = xm;
    }
  protected:
    void callback() {
      try {
        if (where == GRB_CB_MIPSOL) {
          //Encontrei uma possível solução, preciso verificar se o grafo está conexo.
          double *x = new double[m];

          x = getSolution(vars, m);

          vector<int> subset;

          procurar_componentes(n, m, x, subset);

          int len = (int)(subset.size());

          //Minha componente não possui todos os vértices
          if (len < (int)V_especiais) {
            GRBLinExpr expr = 0;
            for(int i = 0; i < (int)subset.size(); i++){
              int x = subset[i];
              for(auto w: grafo[x]){
                //pego todas as arestas que possuem uma ponta em subset e outra fora.
                if(binary_search(subset.begin(), subset.end(), w.v) or !especial[w.v]) continue;
                expr += vars[w.id];
              }
            }
            addLazy(expr >= 2);
          }

          subset.clear();
          delete[] x;
        }
      } catch (GRBException e) {
        cout << "Error number: " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
      } catch (...) {
        cout << "Error during callback" << endl;
      }
    }
    
};
vector<edge> Er;
int dist[N][N];

/*
  Esta função é responsável por simplificar o meu grafo.
  Primeiro, para todas as arestas pertencentes a R, adiciono elas no meu grafo auxiliar Gr = (Vr, R)
  então, seja dist[u][v] = menor caminho entre os vértices (u, v) no grafo original, então cria-se uma aresta entre 
  esses dois vértices com custo igual ao menor caminho.
  Depois removemos arestas (u, v) se existe um k, tal que dist[u][v] = dist[u][k] + dist[k][v], u != k, u != v.

  Complexidade O(N^3)
*/
void Simplificar_Grafo()
{
  const int inf = 0x3f3f3f3f;
  vector<int> vr;
  vector<edge> arestas;
  set<pair<pii,int> > ed;

  //Adiciona arestas requeridas = R
  for(auto w: Er)
  {
    int v = w.v, cost = w.cost, req = w.req, id = w.id, x = w.x;
    grafo[x].push_back({v, cost, req, m, x});
    grafo[v].push_back({x, cost, req, m, v});
    m++;
    ed.insert({{min(x, v),max(x,v)}, cost});
    vr.pb(x), vr.pb(v);
    especial[x]=especial[v]=1;
    arestas.pb(w);
  }

  sort(vr.begin(), vr.end());
  vr.erase(unique(vr.begin(), vr.end()), vr.end());
  memset(dist, 0x3f3f3f3f, sizeof dist);

  for(int i = 0; i < n; i++)
  {
    dist[i][i] = 0;
    for(auto w: G[i])
      dist[i][w.v] = min(w.cost, dist[i][w.v]),
      dist[w.v][i] = min(w.cost, dist[w.v][i]);
  }

  //Calcula as distâncias mínimas.
  for(int k = 0; k < n; k++)
    for(int i = 0; i < n; i++)
      for(int j = 0; j < n; j++)
        dist[i][j] = min(dist[i][j], dist[i][k] + dist[k][j]);

  // remove arestas paralelas e arestas não necessárias
  for(int i = 0; i < (int)vr.size(); i++){
    for(int j = i + 1; j < (int)vr.size(); j++){
      int a = vr[i], b = vr[j];
      if(dist[a][b] == inf) continue;
      edge A; A.v=b,A.cost=dist[a][b],A.req=0,A.id=m,A.x=a;
      edge B; B.v=a,B.cost=dist[b][a],B.req=0,B.id=m,B.x=b;
      int rem =0;
      if(ed.count({{min(a,b),max(a,b)}, dist[a][b]})) continue;

      for(int k = 0; k < (int)vr.size(); k++){
        if(k == i or k == j) continue;
        int vk = vr[k];
        if(dist[a][b] == dist[a][vk] + dist[vk][b]) rem=1;
      }
      if(rem) continue;
      m++;
      grafo[a].pb(A);
      grafo[b].pb(B);
      arestas.pb(A);
    }
  }
  edges = arestas;

  for(int j = 0; j < (int)edges.size(); j++){
    C[j] = edges[j].cost;
    edges[j].id = j;
  }

  V_especiais = (int)vr.size();
}

/*
  N = número de vértices
  req = Número de arestas que pertencem a R
  nreq = Número de arestas que não pertencem a R
*/
int main()
{
  cin>>n>>req>>nreq;

  //Lendo arestas que pertencem a R. liga os vértices(a, b) e tem custo = c
  for(int i = 0, a, b, c; i < req; i++)
  {
    cin>>a>>b>>c;
    edge A; A.v=b,A.cost=c,A.req=1,A.x=a, A.id=i;
    edge B; B.v=a,B.cost=c,B.req=1,B.x=b, B.id = i;
    G[a].push_back(A);
    G[b].push_back(B);
    cost[a] ++;
    cost[b] ++;
    Er.push_back(A);
  }

  //Lendo arestas que não pertencem a R. Liga os vértices (a, b) e tem custo = c
  for(int i = req, a, b, c; i < req+nreq; i++)
  {
    cin>>a>>b>>c;
    edge A; A.v=b,A.cost=c,A.req=0,A.id=i,A.x=a;
    edge B; B.v=a,B.cost=c,B.req=0,B.id=i,B.x=b;
    G[a].push_back(A);
    G[b].push_back(B);
  }

 Simplificar_Grafo();

  GRBEnv *env = NULL;
  GRBVar *x = new GRBVar[m];
  GRBVar *z = new GRBVar[n];

  try
  {
    env = new GRBEnv();
    GRBModel model = GRBModel(*env);

    model.set(GRB_IntParam_LazyConstraints, 1);

    /* x[i] = quantidade de vezes que a aresta "i" será utilizada.
      Se é uma aresta requerida, deve ser usada pelo menos uma vez.
      Õbserve que deve ser usada <= 2 vezes, pois uma aresta (a, b) pode ser usada nos sentidos:
      a->b  e  b->a.
    */

    for(int i = 0; i < m; i++){
      if(edges[i].req) x[i] = model.addVar(1.0, 2.0, edges[i].cost, GRB_INTEGER, "x(" + itos(i) + ")");
      else x[i] = model.addVar(0.0, 2.0, edges[i].cost, GRB_INTEGER, "x(" + itos(i) + ")");
    }

    //z[i] apenas garante a paridade para encontrarmos um caminho fechado
    for(int i = 0; i < n; i++)
      z[i] = model.addVar(0.0, n, 0, GRB_INTEGER, "z(" +itos(i)+")");

    /*
      Restrição (1):
      somatorio( x[e] ) - 2*z[i] = 0
    */
    for(int i = 0; i < n; i++)
    {
      if(!especial[i]) continue;
      GRBLinExpr expr = 0;
      for(auto aresta: grafo[i]){
        expr += x[aresta.id];
      }

      expr -= 2*z[i];

      model.addConstr(expr == 0);
    }

    subtourelim cb = subtourelim(x, n, m);

    model.setCallback(&cb);

    // Optimize model

    model.optimize();

    // imprimir resposta
    if (model.get(GRB_IntAttr_SolCount) > 0) {
      double *sol = new double[m];
      sol = model.get(GRB_DoubleAttr_X, x, m);

      vector<int> subset;

      procurar_componentes(n, m, sol, subset);

      int len = subset.size();

      assert(len == V_especiais);

      double ans = 0;
      for(int i = 0; i < m; i++){
        if ((sol[i] >= 1.0 - eps) or (edges[i].req)){
          ans += sol[i] * edges[i].cost;
        }
      }
      cout<<"ans = "<<ans<<"\n";
    }
  } catch (GRBException e) {
    cout << "Error number: " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch (...) {
    cout << "Error during optimization" << endl;
  }
}