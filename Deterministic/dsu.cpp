#include "dsu.h"
#include <stdio.h>
#include <stdlib.h>
#include <vector>
using namespace std;

struct DSU
{
  vector<int> pai, peso, sz;
};

void build(DSU *UF, int n){
	UF->sz.resize(n);
	UF->pai.resize(n);
	UF->peso.resize(n);
	for(int i = 0; i < n; i++)
	  UF->pai[i] = i, UF->peso[i] = 0, UF->sz[i] = 1;
}

int get_sz(DSU *UF, int x){
	x = Find(UF, x);
	return UF->sz[x];
}
DSU *create_DSU(int n){
	DSU *UF = new DSU;
	build(UF, n);
	return UF;
}

void clear(DSU *UF){
	UF->pai.clear();
	UF->peso.clear();
	UF->sz.clear();
	delete UF;
}

int Find(DSU *UF, int x)
{
	if(UF->pai[x] == x) return x;
	return UF->pai[x] = Find(UF,UF->pai[x]);
}

void join(DSU *UF, int a, int b)
{
	a = Find(UF,a), b = Find(UF,b);
	if(a == b) return;
	if(UF->peso[a] < UF->peso[b]) swap(a, b);
	UF->pai[a] = b;
	UF->sz[b] += UF->sz[a];
	if(UF->peso[a]==UF->peso[b])UF->peso[b]++;
}
