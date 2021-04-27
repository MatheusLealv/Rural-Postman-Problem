#include "dsu.h"
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>

using namespace std;

struct DSU
{
  vector<int> pai, peso, sz;
  vector<bool> has_good;
};

void build(DSU *UF, int n){
	UF->sz.resize(n+1);
	UF->pai.resize(n+1);
	UF->peso.resize(n+1);
	UF->has_good.resize(n + 1);
	for(int i = 0; i < n+1; i++)
	  UF->pai[i] = i, UF->peso[i] = 0, UF->sz[i] = 1, UF->has_good[i] = 0;
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
	UF->has_good.clear();
	delete UF;
}

int Find(DSU *UF, int x){
	if(UF->pai[x] == x) return x;
	return UF->pai[x] = Find(UF,UF->pai[x]);
}

void join(DSU *UF, int a, int b){
	a = Find(UF,a), b = Find(UF,b);
	if(a == b) return;
	if(UF->peso[a] < UF->peso[b]) swap(a, b);
	UF->pai[a] = b;
	if(UF->has_good[a])UF->has_good[b] = UF->has_good[a];
	UF->sz[b] += UF->sz[a];
	if(UF->peso[a]==UF->peso[b])UF->peso[b]++;
}

bool hasGood(DSU *UF, int x){
	x = Find(UF, x);
	return UF->has_good[x];
}

void setGood(DSU *UF, int x){
	UF->has_good[x] = 1;
}