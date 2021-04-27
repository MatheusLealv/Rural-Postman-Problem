/*
  Stochastic Rural Postman Problem
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
#define eps (1e-6)
#define rep(i, a, b) for(int i = a; i < (b); ++i)
#define all(x) begin(x), end(x)
#define sz(x) (int)(x).size()
typedef long long ll;
typedef pair<int, int> pii;
typedef vector<int> vi;

string itos(int i) {stringstream s; s << i; return s.str(); }
int qtdCenarios = 0;

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

vi eulerWalk(vector<vector<pii>>& gr, int nedges, int src=0) {
  int n = sz(gr);
  vi D(n), its(n), eu(nedges), ret, s = {src};
  D[src]++;
  while (!s.empty()) {
    int x = s.back(), y, e, &it = its[x], end = sz(gr[x]);
    if (it == end){ ret.push_back(x); s.pop_back(); continue; }
    ++it;
    y=gr[x][it].f;
    e=gr[x][it].s;
    if (!eu[e]) {
      D[x]--, D[y]++;
      eu[e] = 1; s.push_back(y);
    }}
  for (int x : D) if (x < 0 || sz(ret) != nedges+1) return {};
  return {ret.rbegin(), ret.rend()};
}

struct Cenario{
  double probability;
  int n, req, nreq, m, V_especiais;
  vector<int> cost, C, especial;
  vector< vector<edge> > G, grafo;
  vector<edge> edges, Er;
  vector<int> isGood;

  vector<int> sol_x, sol_y;

  GRBVar *x;// = new GRBVar[S[cenarioID].m];
  GRBVar *y;// = new GRBVar[S[cenarioID].m];
  GRBVar *z;// = new GRBVar[S[cenarioID].n];
  double *sol;
  double *solY;
  Cenario(double prob_ = 1.0, int n_ = 0, int req_ = 0, int nreq_ = 0){
    n = n_;
    req = req_;
    nreq = nreq_;
    probability = prob_;
    m = 0;
    cost.resize(req + nreq + 1), C.resize(req + nreq + 1), especial.resize(n+1);
    G.resize(n + 1), grafo.resize(n + 1);
    V_especiais = 0;
    
    //*******************
    isGood.resize(n + 1);
    z = new GRBVar[n_];
    x = new GRBVar[req_ + nreq_];
    y = new GRBVar[req_ + nreq_];

    // cout<<"CREATE "<<req_ + nreq_<<"\n";
    sol = new double[req_ + nreq_];
    solY = new double[req_ + nreq_];
  }

  void addReqEdge(int a, int b, int c, int i){
    edge A; A.v=b,A.cost=c,A.req=1,A.x=a, A.id=i;
    edge B; B.v=a,B.cost=c,B.req=1,B.x=b, B.id = i;
    G[a].push_back(A);
    G[b].push_back(B);
    cost[a] ++;
    cost[b] ++;
    Er.push_back(A);

    edges.pb(A); // **********
    m++;
  }

  void addNoReqEdge(int a, int b, int c, int i){
    edge A; A.v=b,A.cost=c,A.req=0,A.id=i,A.x=a;
    edge B; B.v=a,B.cost=c,B.req=0,B.id=i,B.x=b;
    G[a].push_back(A);
    G[b].push_back(B);

    edges.pb(A); //**************
    m++;
  }

  /*
    Dada uma possível solução válida (de acordo com a restrição (1)), encontrar a componente
    ue possui tamanho mínimo. A componente é inserida no vector subset.

    Complexidade O(N * alpha(N))
  */
  void procurar_componentes(double *x, double *y, vector<int> &subset, DSU *&UF){
    for(int i = 0; i < m; i++){
      if((x[i] + y[i]>= 1.0 - eps) or edges[i].req == 1)
       join(UF, edges[i].x, edges[i].v);
    }

    int menor = INF, opt = -1;

    for(int i = 0; i < n; i++)
      if(especial[i]){
        menor = min(menor, get_sz(UF, i));
      }

    for(int i = 0; i < n; i++)
      if(especial[i] and get_sz(UF, i) == menor){
        opt = Find(UF, i);
        break;
      }

    for(int i = 0; i < n; i++)
      if(especial[i] and Find(UF, i) == opt) subset.push_back(i);
  }

  void get_base_solution(){
    sol_x.resize(sz(edges));
    sol_y.resize(sz(edges));
    for(int i=0;i<sz(sol_x);i++)sol_x[i]=0;
    for(int i = 0; i < sz(sol_y);i++) sol_y[i]=0;

    DSU *UF = create_DSU(n);
    vector<int> grau(n + 1);
    int tot=0;
    vector<vector<pii>> graph(n+1);
    auto E = edges;
    sort(all(E), [&](edge A, edge B){
      return A.cost<B.cost;
    });
    for(auto e: E){
      if(e.req or (Find(UF, e.x) != Find(UF, e.v))){
        grau[e.x] ++;
        grau[e.v] ++;
        join(UF, e.x, e.v);
        sol_x[e.id] ++;
        sol_x[e.id] ++;
        graph[e.x].pb({e.v, 2*e.id});
        graph[e.v].pb({e.x, 2*e.id});
        tot=max(tot, 2*e.id);
      }
    }

  }
};
vector<Cenario> S;

/*
  Callback. Sempre que uma possível solução é encontrada, pego sua componente de tamanho mínimo e 
  adiciona a restrição (2) para essa componente em relação aos outros vértices.

  Complexidade O(N*log(N))
*/

class subtourelim: public GRBCallback
{
  public:
    subtourelim()
    {

    }
  protected:
    void callback() {
      try {
        if (where == GRB_CB_MIPSOL) {
          //Encontrei uma possível solução, preciso verificar se o grafo está conexo.
          for(int cenarioID = 1; cenarioID <= qtdCenarios; cenarioID++){
            int n = S[cenarioID].n, m = S[cenarioID].m;

            double *x = new double[m];
            double *y = new double[m];
            x = getSolution(S[cenarioID].x, m);
            y = getSolution(S[cenarioID].y, m);

            vector<int> subset;
            DSU *UF = create_DSU(n);
            build(UF, n);
            S[cenarioID].procurar_componentes(x, y, subset, UF);

            int len = (int)(subset.size());

            int cont = 0;
            for(auto vertice: subset){
              if(S[cenarioID].isGood[vertice]){
                cont ++;
              }
            }
            //Minha componente não possui todos os vértices
            if (cont < S[cenarioID].V_especiais) {
              GRBLinExpr expr = 0;
              for(int i = 0; i < (int)subset.size(); i++){
                int x = subset[i];
                for(auto w: S[cenarioID].grafo[x]){
                  //pego todas as arestas que possuem uma ponta em subset e outra fora.
                  assert(w.v != x);
                  if(binary_search(subset.begin(), subset.end(), w.v)) continue;
                  expr += S[cenarioID].x[w.id];
                  expr += S[cenarioID].y[w.id];
                }
              }
              addLazy(expr >= 2);
            }
            clear(UF);
            subset.clear();
            delete[] x;
          }
        }
      } catch (GRBException e) {
        cout << "#1 Error number: " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
      } catch (...) {
        cout << "Error during callback" << endl;
      }
    }
    
};

void solveModel(int qtdCenarios){
 // S[cenarioID].simplificarGrafo(); *******************


  // ****************
  for(int cenarioID = 0; cenarioID <= qtdCenarios; cenarioID ++){
    S[cenarioID].grafo = S[cenarioID].G;
    for(auto edge: S[cenarioID].Er){
      int a = edge.x, b = edge.v;
       S[cenarioID].isGood[a] = 1;
       S[cenarioID].isGood[b] = 1;
    } 
    for(int i = 0; i < S[cenarioID].n; i++){
      if(S[cenarioID].isGood[i]){
        S[cenarioID].V_especiais ++;
        S[cenarioID].especial[i] = 1;
      }
    }
  }

  // ******************8

  GRBEnv *env = NULL;

  try
  {
    env = new GRBEnv();
    GRBModel model = GRBModel(*env);
    model.set(GRB_IntParam_LazyConstraints, 1);

    for(int cenarioID = 0; cenarioID <= qtdCenarios; cenarioID++){
      cout<<"CRIANDO VARIAVEIS DO CENÁRIO "<<cenarioID<<"\n";
      if(cenarioID > 0){
        S[cenarioID].get_base_solution();
      }
      for(int i = 0; i < S[cenarioID].m; i++){
        string varName = "x(" + itos(i) + "," + itos(cenarioID) + ")";
        string varName2 = "y(" + itos(i) + "," + itos(cenarioID) + ")";
        double P = S[cenarioID].probability;

         S[cenarioID].x[i] = model.addVar(0.0, 2.0, P*S[cenarioID].edges[i].cost, GRB_INTEGER, varName);

        //SETTTTT
        // if(cenarioID > 0) S[cenarioID].x[i].set(GRB_DoubleAttr_Start, sol_x[i]);

       if(cenarioID > 0) S[cenarioID].y[i] = model.addVar(0.0, 2, 0, GRB_INTEGER, varName2);

        //model.update(); //***
        
        if(cenarioID > 0){
          GRBLinExpr expr = S[cenarioID].y[i] - S[0].x[i];
          model.addConstr(expr <= 0, "C(" + itos(i) + "," + itos(cenarioID) + ")");
        }

        if(S[cenarioID].edges[i].req){
          GRBLinExpr expr = S[cenarioID].x[i] + S[cenarioID].y[i];
          model.addConstr(expr >= 1.0 - eps, "C2(" + itos(i) + "," + itos(cenarioID) + ")");
        }
        // cout<<"CRIA VAR "<<cenarioID<<" "<<varName<<"\n";
      }
      //z[i] apenas garante a paridade para encontrarmos um caminho fechado


      //model.update(); //***
      /*
        Restrição (1):
        somatorio( x[e] ) - 2*z[i] = 0
      */
      if(cenarioID > 0){
        for(int i = 0; i < S[cenarioID].n; i++){
         string varName = "z(" + itos(i) + "," + itos(cenarioID) + ")";
         S[cenarioID].z[i] = model.addVar(0.0, sz(S[cenarioID].grafo[i]), 0, GRB_INTEGER, varName);
        }
        for(int i = 0; i < S[cenarioID].n; i++) {
         // if(!S[cenarioID].especial[i]) continue; **************
          GRBLinExpr expr = 0;
          for(auto aresta: S[cenarioID].grafo[i]){
            expr += S[cenarioID].x[aresta.id];
            expr += S[cenarioID].y[aresta.id];
          }

          expr -= 2*S[cenarioID].z[i];

          model.addConstr(expr == 0, "C3(" + itos(i) + "," + itos(cenarioID) + ")");
        }
       //model.update(); //***
      }
      S[cenarioID].sol_y.resize(S[cenarioID].m);
    }

    model.update();
    cout<<"\n---- TODAS AS VARIÁVEIS FORAM CRIADAS ----\n";
    cout<<"\n---- COMPRANDO ARESTAS NO PRIMEIRO ESTÁGIO ----\n";

    S[0].sol_x.resize(S[0].m);
      vector<int> visitada(S[0].m);
    while(true){
      
      bool take = false;
      for(int i = 0; i < S[0].m; i++){
        if(visitada[i]) continue;

        double buy0 = 0.0, buy2 = 2*S[0].edges[i].cost, buy1=S[0].edges[i].cost;
        for(int cenarioID = 1; cenarioID <= qtdCenarios; cenarioID ++){
          buy1 -= S[cenarioID].probability*min(1.0, 1.0*S[cenarioID].sol_x[i])*S[cenarioID].edges[i].cost;
          buy2 -= S[cenarioID].probability*min(2.0, 1.0*S[cenarioID].sol_x[i])*S[cenarioID].edges[i].cost;
        }
        double mn = min({buy0, buy1, buy2});
        if(abs(mn) <= eps) continue;
        if(abs(buy2-mn) <= eps or abs(buy1-mn) <= eps){
          visitada[i]=1;
          take = true;
          int q = (abs(buy2-mn) <= eps?2:1);
          S[0].sol_x[i] = q;

          for(int cenarioID = 1; cenarioID <= qtdCenarios; cenarioID ++){
              S[cenarioID].sol_y[i] = min(S[cenarioID].sol_x[i], q);
              S[cenarioID].sol_x[i] = S[cenarioID].sol_x[i] - S[cenarioID].sol_y[i]; 
          }
        }
      }
      if(!take) break;
    }
    // for(int cenarioID = 1; cenarioID <= qtdCenarios; cenarioID ++) S[cenarioID].get_base_solution();
  cout<<"---- ARESTAS DO PRIMEIRO ESTÁGIO DECIDIDAS ----\n";
  for(int cenarioID = 0; cenarioID <= qtdCenarios; cenarioID ++){
    cout<<"SET CENÁRIO #"<<cenarioID<<"\n";
    for(int i = 0; i < S[cenarioID].m; i++){
     S[cenarioID].x[i].set(GRB_DoubleAttr_Start, S[cenarioID].sol_x[i]);
     if(cenarioID > 0)S[cenarioID].y[i].set(GRB_DoubleAttr_Start, S[cenarioID].sol_y[i]);
    }
  }  

  cout<<"--- RESOLVENDO O MODELO ---\n";

    model.update();
    subtourelim cb = subtourelim();
    model.setCallback(&cb);

    // Optimize model

    model.write("out333.lp");
    model.optimize();
    double expectedValue = 0.0;

    // imprimir resposta
    if (model.get(GRB_IntAttr_SolCount) > 0) {
      
      for(int cenarioID = 0; cenarioID <= qtdCenarios; cenarioID ++){
        S[cenarioID].sol = model.get(GRB_DoubleAttr_X, S[cenarioID].x, S[cenarioID].m);
        S[cenarioID].solY = model.get(GRB_DoubleAttr_X, S[cenarioID].y, S[cenarioID].m);

        vector<int> subset;

        DSU *UF = create_DSU(S[0].n);
        S[cenarioID].procurar_componentes(S[cenarioID].sol, S[cenarioID].solY, subset, UF);

        int len = subset.size();
        int conta = 0;
        for(auto vertice: subset){
          if(S[cenarioID].isGood[vertice]){
            conta ++;
          }
        }

        assert(conta == S[cenarioID].V_especiais);

        for(int i = 0; i < S[cenarioID].m; i++){
              expectedValue += S[cenarioID].probability*S[cenarioID].sol[i] * S[cenarioID].edges[i].cost;
        }
      }
      cout<<"ans = "<<expectedValue<<"\n";
    }
  }
  catch (GRBException e) {
    cout << "#2 Error number: " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch (...) {
    cout << "Error during optimization" << endl;
  }
}

void input(){
  int n, m, inflacao;
  string inutil;
  cin>>inutil>>n;
  cin>>inutil>>m;
  cin>>inutil>>qtdCenarios;
  cin>>inutil>>inflacao>>inutil;
  cin>>inutil>>inutil;
  Cenario S0 = Cenario(1.00, n, 0, m);
  S.push_back(S0);

  for(int j = 0, a, b, c; j < m; j++){
    cin>>a>>b>>c;
    a--,b--;
    S[0].addNoReqEdge(a, b, c, j);
  }

  /*
  exemplo:
  CENÁRIO 1
  PROBABILIDADE: 599/3905 = 0.153393
  QTD_ARESTAS_REQUERIDAS: 282
  QTD_ARESTAS_NAO_REQUERIDAS: 315
  */
  for(int cenario = 1; cenario <= qtdCenarios; cenario ++){
    double prob;
    int req, noreq;
    cin>>inutil>>inutil;
    cin>>inutil>>inutil>>inutil>>prob;
    cin>>inutil>>req;
    cin>>inutil>>noreq;
    Cenario Si = Cenario(prob, n, req, noreq);
    S.push_back(Si);

    for(int j = 0; j < m; j++){
      int a, b, c, tipo;
      cin>>a>>b>>c>>tipo;
      a--,b--;
      if(tipo) S[cenario].addReqEdge(a, b, c, j);
      else S[cenario].addNoReqEdge(a, b, c, j);
    }
  }
}


int main(){
  input();
   // S[0].simplificarGrafo();
  solveModel(qtdCenarios);

}