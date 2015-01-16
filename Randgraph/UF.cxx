
#include "UF.h"

void UnionFind::init(int num_elmts) {      // n is number of UF elements
  int i;

  A = new UF_Elmt[num_elmts];
  n=num_elmts;

  for(i=0;i<n;i++) {
    A[i].parent = -1;     // i is a root, a singleton set
    A[i].rank   = 0;
  }
}

UnionFind::UnionFind() {
  A = 0;
  n = 0;
}

UnionFind::~UnionFind() {
  if (n > 0)
    delete[] A;
}

void UnionFind::print() {
  int i;

  for(i=0;i<n;i++) {
    cout << i << " --> " << A[i].parent << "   (rank " << (int) A[i].rank << ")\n";
  }
}

// no error detection.  assumes i1 and i2 are roots
// returns the new root
int UnionFind::unionsets(int i1, int i2) {
  if(A[i1].rank <= A[i2].rank) {
    A[i1].parent = i2;
    if(A[i1].rank == A[i2].rank) A[i2].rank++;
    return(i2);
  } else {
    A[i2].parent = i1;
    return(i1);
  }
}

// Halves the path.
int UnionFind::find(int i) {
  int p,gp;

  while(A[i].parent != -1) {
    p = A[i].parent;
    if(A[p].parent != -1) gp = A[p].parent; else return(p); 
                                                   // an odd length path will
    A[i].parent = gp;                              // terminate here
    i = gp;
  }
  return i;     // an even length one here
}

// true if root of i is the same as root of j.
int UnionFind::equiv(int i, int j) {
  int ri = find(i);
  int rj = find(j);

  if (ri == rj) return(1);   // return true if i and j are equivalent
  unionsets(ri,rj);         // otherwise, union them
  return(0);                 // and return false
}

/*
int main() {
  UnionFind UF;

  UF.allocate(5);
  UF.unionsets(0,1);
  UF.unionsets(2,3);
  UF.find(1);
  UF.print();
  cout << UF.equiv(0,4) << "\n";
  cout << UF.equiv(0,4) << "\n";
  UF.print();

}
*/



