#include<iostream>

// A union-find element keeps its rank and the id of its parent.

struct UF_Elmt {
  int parent;                // -1 if root
  char rank;                  // rank
};

class UnionFind {
 public:
  UnionFind();
  void init(int n);       // allocates n singleton sets
  ~UnionFind();          // deallocates
  int unionsets(int i1, int i2);
  int find(int i);
  void print();
  
  //returns 1 (true) if i and j are in the same set and 0 (false)
  //otherwise, in which case a union is performed on the sets of i and j.

  int equiv(int i, int j);

 private:
  UF_Elmt *A;
  int n;
};
