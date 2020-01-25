#ifndef DSU_H
#define DSU_U

typedef struct DSU DSU;

void build(DSU *UF, int n);

DSU *create_DSU(int n);

void clear(DSU *UF);

int get_sz(DSU *UF, int x);

int Find(DSU *UF, int x);

void join(DSU *UF, int a, int b);

#endif
