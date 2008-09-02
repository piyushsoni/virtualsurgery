// GLOBAL VARIABLES
int showSkeletOnly=1;
int rings=1;                               // number of rings for colorcoding
int r=1;                                // radius of spheres for displaying vertices
boolean showNormals=false, showVertices=false, showEdges=false, showTriangles=true,  showSelectedTriangle=true, showLabels=false, showPath=false;  // flags for rendering
boolean showSkeleton=true, showSelectedLake=true, showOtherLakes=true, showDistance=false, showEB=true, showEBrec=true, showClusters=true;
//Mesh M = new Mesh();     // creates a default triangle mesh
int MaxMeshes =10;
Mesh[] Meshes = new Mesh[MaxMeshes];
Mesh Current;
int currentMeshId=0, MeshesInstantiaded = 1;



//========================== class MESH ===============================
class Mesh {

//  ==================================== INIT, CREATE, COPY ====================================
  Mesh(int MeshID) {meshId = MeshID;}
 int meshId; 
 void declare() {
   for (int i=0; i<maxnv; i++) {G[i]=new pt(0,0,0); Nv[i]=new vec(0,0,0);};   // init vertices and normals
   for (int i=0; i<maxnt; i++) {Nt[i]=new vec(0,0,0);  };       // init triangle normals and skeleton lab els
   }
 void makeGrid (int w) { // make a 2D grid of vertices
  for (int i=0; i<w; i++) {for (int j=0; j<w; j++) { G[w*i+j].setTo(height*.8*j/(w-1)+height/10,height*.8*i/(w-1)+height/10,0);}}    
  for (int i=0; i<w-1; i++) {for (int j=0; j<w-1; j++) {                  // define the triangles for the grid
    V[(i*(w-1)+j)*6]=i*w+j;       V[(i*(w-1)+j)*6+2]=(i+1)*w+j;       V[(i*(w-1)+j)*6+1]=(i+1)*w+j+1;
    V[(i*(w-1)+j)*6+3]=i*w+j;     V[(i*(w-1)+j)*6+5]=(i+1)*w+j+1;     V[(i*(w-1)+j)*6+4]=i*w+j+1;}; };
  nv = w*w;
  nt = 2*(w-1)*(w-1); 
  nc=3*nt;
  }
 void init() {init(true);}
 void init(boolean doOTable) { 
   for (int i=0; i<maxnt; i++) {Nt[i]=new vec(0,0,0); visible[i]=true;}; 
   for (int v=0; v<maxnv; v++) {col[v] = grey; } 
   if (doOTable) computeO(); 
   computeValenceAndResetNormals(); 
   computeTriNormals(); 
   computeVertexNormals(); 
   c=0; sc=0;  
  }
  // ============================================= DISPLAY =======================================
pt Cbox = new pt(width/2,height/2,0);                   // mini-max box center
float Rbox=1000;                                        // Radius of enclosing ball
boolean showEdges=false;
boolean showDistance=false;
boolean showNormals=false;
boolean showPath=false;
boolean showVertices=false;
boolean showEBrec=false;
boolean showTriangles=true;
boolean showColor = false;
boolean showSelectedVertex=false;
void computeBox() {
  pt Lbox =  G[0].make();  pt Hbox =  G[0].make();
  for (int i=1; i<nv; i++) { 
    Lbox.x=min(Lbox.x,G[i].x); Lbox.y=min(Lbox.y,G[i].y); Lbox.z=min(Lbox.z,G[i].z);
    Hbox.x=max(Hbox.x,G[i].x); Hbox.y=max(Hbox.y,G[i].y); Hbox.z=max(Hbox.z,G[i].z); 
    };
  Cbox.setToPoint(midPt(Lbox,Hbox));  Rbox=Cbox.disTo(Hbox);
//  println("Box r="+Rbox); print("C="); Cbox.write(); print("L="); Lbox.write(); print("H="); Hbox.write();
  };
pt cg(int c) {pt cPt = midPt(g(c),midPt(g(c),triCenter(t(c))));  return(cPt); };   // computes point at corner
void showCorner(int c, int r) {pt cPt = midPt(g(c),midPt(g(c),triCenter(t(c))));  cPt.show(r); };   // renders corner c as small ball
void showCornerAndNormal(int c, int r) {pt cPt = midPt(g(c),midPt(g(c),triCenter(t(c))));  noStroke(); cPt.show(r); 
                     stroke(magenta); vec N = Nt[t(c)].make(); N.makeUnit();  N.mul(10*r);  N.show(cPt);};   // renders corner c as small ball
void shade(int i) {if(visible[i]) { beginShape(TRIANGLES); G[V[3*i]].vert(); G[V[3*i+1]].vert(); G[V[3*i+2]].vert(); endShape(); };}; // shade tris

void showTriNormals() {vec N= new vec(0,0,0); for (int i=0; i<nt; i++) { N.setToVec(Nt[i]); N.makeUnit();  N.mul(5*r);  N.show(triCenter(i)); };  };
void showVertexNormals() {vec N= new vec(0,0,0); for (int i=0; i<nv; i++) {N.setToVec(Nv[i]); N.makeUnit(); N.mul(5*r);  N.show(G[i]); };  };
void show() {
  //int col=60;
  noSmooth(); noStroke();   smooth();
   if(showDistance) for(int t=0; t<nt; t++) {fill(60,120,(rings-Mt[t])*120/rings); shade(t);};  
  // if(showTriangles)  if(!showDistance) {fill(cyan); for(int t=0; t<nt; t++)  shade(t); noFill();};  

     if(showTriangles && !showDistance && !showColor) 
     {
       beginShape(TRIANGLES);  
       if(Current.meshId==meshId)
       {
         for(int t=0; t<nt; t++) 
         {
           if (hitTriangle[t]) fill(red); 
           else fill(metal); 
           G[V[3*t]].vert(); 
           G[V[3*t+1]].vert(); 
           G[V[3*t+2]].vert();
         }
       }
       else
       {
         for(int t=0; t<nt; t++) 
         {
           if (hitTriangle[t]) fill(red); 
           else fill(grey); 
           G[V[3*t]].vert(); 
           G[V[3*t+1]].vert(); 
           G[V[3*t+2]].vert();
         }         
       }
       endShape();
     }
   
   if(showTriangles && !showDistance &&  showColor) {beginShape(TRIANGLES);  for(int c=0; c<nc; c++) {fill(col[V[c]]); G[V[c]].vert();} endShape();}

  if (showEdges) {stroke(dblue); for(int i=0; i<nc; i++) if(visible[t(i)]) drawEdge(i); };  
  stroke(red); showBorder();
  if (showVertices) {noStroke(); noSmooth();fill(blue); for (int v=0; v<nv; v++)  G[v].show(r); noFill();};
  if (showLabels) { fill(black); for (int i=0; i<nv; i++) {G[i].label(str(i),labelD); }; };
  if (showNormals) {stroke(blue); showTriNormals(); stroke(magenta); showVertexNormals(); };                // show triangle normals
  if (showSelectedTriangle) {
    noStroke(); fill(green); shade(t(c)); 
   // fill(blue); stroke(dblue); showCornerAndNormal(c,r);  fill(cyan); noStroke(); showCornerAndNormal(prevc,r);  
    }; 
  if (showSelectedVertex) {
     noStroke(); fill(black); g(c).show(r);  
  }
 // fill(dred); mark.show(r); // stroke(red); line(eye.x,eye.y,eye.z,mark.x,mark.y,mark.z);
  }

//  ==========================================================  VERTICES ===========================================
 int maxnv = 45000;                         //  max number of vertices
 int nv = 0;                              // current  number of vertices
 pt[] G = new pt [maxnv];                   // geometry table (vertices)
 vec[] Nv = new vec [maxnv];                 // vertex normals or laplace vectors
 int[] Mv = new int[maxnv];                  // vertex markers
 int [] Valence = new int [maxnv];          // vertex valence (count of incident triangles)
 boolean [] Border = new boolean [maxnv];   // vertex is border
 boolean [] VisitedV = new boolean [maxnv];  // vertex visited
 
 color[] col = new color[maxnv];
 
void computeVertexNormals() {  // computes the vertex normals as sums of the normal vectors of incident tirangles scaled by area/2
  for (int i=0; i<nv; i++) {Nv[i].setTo(0,0,0);};  // resets the valences to 0
  for (int i=0; i<3*nt; i++) {Nv[v(i)].add(Nt[t(i)]);};
  for (int i=0; i<nv; i++) {Nv[i].makeUnit();}; 
  };
int addVertex(pt P) { G[nv].setTo(P); nv++; return nv-1;};
int addVertex(float x, float y, float z) { G[nv].x=x; G[nv].y=y; G[nv].z=z; nv++; return nv-1;};
void move() {g().addScaledVec(pmouseY-mouseY,N());}
//  ==========================================================  EDGES ===========================================
void findShortestEdge() {  // assumes manifold
  float md=d(g(p(c)),g(n(c))); int ma=c;
  for (int a=0; a<nc; a++) if (vis(a)&&(d(g(p(a)),g(n(a)))<md)) {ma=a; md=d(g(p(a)),g(n(a)));}; 
  c=ma;
  } 
void drawEdge(int c) {line(g(p(c)).x,g(p(c)).y,g(p(c)).z,g(n(c)).x,g(n(c)).y,g(n(c)).z); };  // draws edge of t(c) opposite to corner c
void showBorder() {for (int i=0; i<nc; i++) {if (visible[t(i)]&&b(i)) {drawEdge(i);}; }; };         // draws all border edges

//  ==========================================================  TRIANGLES ===========================================
 int maxnt = maxnv*2;                       // max number of triangles
 int nt = 0;                   // current and max number of triangles
 vec[] Nt = new vec [maxnt];                // triangles normals
 boolean[] visible = new boolean[maxnt];    // set if triangle visible
 int[] Mt = new int[maxnt];                 // triangle markers for distance and other things   
 boolean [] VisitedT = new boolean [maxnt];  // triangle visited
pt triCenter(int i) {return(triCenterFromPts( G[V[3*i]], G[V[3*i+1]], G[V[3*i+2]] )); };  pt triCenter() {return triCenter(t());}  // computes center of triangle t(i) 
vec triNormal(int i) { return(triNormalFromPts(G[V[3*i]], G[V[3*i+1]], G[V[3*i+2]])); };  vec triNormal() {return triNormal(t());} // computes triangle t(i) normal * area / 2
void computeTriNormals() {for (int i=0; i<nt; i++) {Nt[i].setToVec(triNormal(i)); }; };             // caches normals of all tirangles
void writeTri (int i) {println("T"+i+": V = ("+V[3*i]+":"+v(o(3*i))+","+V[3*i+1]+":"+v(o(3*i+1))+","+V[3*i+2]+":"+v(o(3*i+2))+")"); };
void hitTriangle() {
   prevc=c;       // save for geodesic 
   float smallestDepth=10000000;
  boolean hit=false;
  for (int t=0; t<nt; t++) {
    if (rayHitTri(eye,mark,g(3*t),g(3*t+1),g(3*t+2))) {
      hit=true;
      float depth = rayDistTriPlane(eye,mark,g(3*t),g(3*t+1),g(3*t+2));
      if ((depth>0)&&(depth<smallestDepth)) {smallestDepth=depth;  c=3*t;};
      }; 
    };
  if (hit) {
    pt X = eye.make(); X.addScaledVec(smallestDepth,eye.vecTo(mark));
    mark.setToPoint(X);
    float distance=X.disTo(g(c));
    int b=c;
    if (X.disTo(g(n(c)))<distance) {b=n(c); distance=X.disTo(g(b)); };
    if (X.disTo(g(p(c)))<distance) {b=p(c);};
    c=b;
    println("c="+c+", pc="+prevc+", t(pc)="+t(prevc));
    };
  }
void addTriangle(int i, int j, int k) {V[nc++]=i; V[nc++]=j; V[nc++]=k; nt++;}

boolean[] hitTriangle = new boolean[maxnt];
void hitTriangle2() {
   pt P = eye, Q = S(P, S(2, V(P, M(g(p(c)), g(n(c))))));
   for (int t = 0; t < nt; t++) hitTriangle[t] = edgeHitTriangle(P,Q,t); //rayHitTri(P,Q,g(3*t),g(3*t+1),g(3*t+2)); 
}

boolean edgeHitTriangle(pt P, pt Q, int t) {
   int vA = v(3*t), vB = v(3*t+1), vC = v(3*t+2);
   boolean sP, sQ, sA, sB, sC;
   sQ=(tetVolume(P, G[vA], G[vB], G[vC])<0); sP=(tetVolume(Q, G[vA], G[vB], G[vC])>0); 
   sA=(vC>vB != tetVolume(P, Q, G[min(vB, vC)], G[max(vB, vC)])>0);  
   sB=(vC>vA != tetVolume(P, G[min(vA, vC)], Q, G[max(vA, vC)])>0);
   sC=(vB>vA != tetVolume(P, G[min(vA, vB)], G[max(vA, vB)], Q)>0);
   return (sQ == sP && sQ == sA && sQ == sB && sQ == sC);
  }
  
  // =========================================== TRIANGLE WALK =======================================
  
  int startT = 0;
  pt rayP = null, rayQ = null;
  
  void setRay() {
    C.setMark();
    rayP = P(eye); rayQ = P(mark);
  }
  
  void walkToRay() {
    if (startT > nt) startT = 0;
    if (rayP == null) return;
    for (int t = 0; t < nt; t++) hitTriangle[t] = false;
    int t = walkToRay(startT, rayP, rayQ, true);
    
    int c = 3*startT; pt P = rayP, Q = rayQ;
    pt M = M(g(c), g(n(c)), g(p(c))); vec Np = C(V(P, M), V(P, Q)); vec U = U(C(Np, V(P, Q)));
    pt A = T(P, dot(U, V(P, M)), U); pt B = S(A, V(P, Q));
    fill(blue, 20); noStroke();
    beginShape(TRIANGLES);
       A.vert(); P.vert(); Q.vert(); 
       A.vert(); Q.vert(); B.vert();
    endShape();
  }
  
  
  // Starting at triangle <start> walk to a triangle that intersects with ray P+sPQ
  int walkToRay(int start, pt P, pt Q, boolean debug) {
     // declare variables. M = midpoint of start triangle. Np = normal vector of plane spanned by PM and PQ
     int c = start*3; vec PQ = V(P, Q); pt M = M(g(c), g(n(c)), g(p(c))); vec Np = C(V(P, M), PQ);
     if (debug) hitTriangle[start] = true; 

     boolean sc = dot(V(P, g(c)), Np) < 0, sn = dot(V(P, g(n(c))), Np) < 0, sp = dot(V(P, g(p(c))), Np) < 0;
     
     if (!sn && sp) c = o(c);  if (!sp && sc) c = l(c);   if (!sc && sn) c = r(c);  // chose direction for first triangle
     
     boolean fail = false;
     while(mixed(PQ, V(P, g(n(c))), V(P, g(p(c)))) < 0) {             // test if we reached the ray
        if (debug) hitTriangle[t(c)] = true; 
        if (t(c) == start) {fail = true; break;}                      // test for endless loop (no straight walk exists)
        if (dot(V(P, g(c)), Np) < 0) c = l(c); else c = r(c);         // advance to next triangle 
     }
     if (fail) return -1; else return t(o(c));
  } 
  
  // If there's no straight line walk from start to an intersecting triangle
  // try starting at a number of randomly selected triangles
  int walkToRaySafe(int start, pt P, pt Q, boolean debug) {
    int r, t = start, n = 0; 
    while(-1 == (r = walkToRay(t, P, Q, debug)) && n < 20) {t = int(random(nt)); n++;} 
    return r;
  }
  
  int lastT = 0;
  void selectTriangle() {
    if (lastT > nt) lastT = 0;
    int t, st = lastT, n = 0; float depth;
    do {
      t = walkToRay(st, eye, mark, false);
      if (t != -1 ) depth = rayDistTriPlane(eye,mark,g(3*t),g(3*t+1),g(3*t+2))*d(eye,mark); else depth = 0;
      //println(depth + ": " + d(mouse3D, eye) + "+" + (Rbox*0.03));
      st = int(random(nt)); n++;
    } while ((t == -1 || depth > d(mouse3D, eye) + Rbox*0.03) && n < 20); 
    if (t != -1 ) {c = t*3; lastT = t; }
  }
  
  
 color lerpColorRGB(color c1, color c2, float amt) {
     int r = int(red(c1) * (1-amt) + red(c2) * amt); 
     int g = int(green(c1) * (1-amt) + green(c2) * amt); 
     int b = int(blue(c1) * (1-amt) + blue(c2) * amt);

     return 255 + 256 * (0 + int(b) + 256*(0 + int(g) + 256*(0 + int(r) + 0)));
  }
  
  void paint(color clr) {
      IntQueue q = new IntQueue(); IntSet closed = new IntSet();
      q.add(c); q.add(n(c)); q.add(p(c));
      while (!q.isEmpty()) {
         int c0 = q.removeFirst(), c = c0; float d = decay(d(g(c), mouse3D)); 
         col[v(c)] = lerpColor(col[v(c)], clr, d*0.02); //println(d);
         do {
            if (!closed.contains(n(c)) && decay(d(g(n(c)), mouse3D)) > 0) {q.add(n(c)); closed.add(n(c));}
            c = s(c); 
         } while (c != c0);
      }
   }
  
   float decay(float d) {
       return constrain(1-10*d/Rbox, 0, 1);
   }  
  
  
  void floodDelete() {
      IntQueue q = new IntQueue(); q.add(t(c));
      while (!q.isEmpty()) {
         int t = q.removeFirst();
         for (int c = 3*t; c < 3*t+3; c++) {int t1 = t(o(c)); if (!hitTriangle[t1] && visible[t1]) {visible[t1] = false; q.add(t1);}} 
      }
      excludeInvisibleTriangles();
      compactVO();
      compactV();
      fanHoles();
  }
  
// ============================================= CORNER OPERATORS =======================================
 int nc = nt*3;                             // current number of corners (3 per triangle)
 int c = 0;                                 // current corner shown in image and manipulated with keys: n, p, o, l, r
 int sc=0;                                  // saved value of c
 int[] V = new int [3*maxnt];               // V table (triangle/vertex indices)
 int[] O = new int [3*maxnt];               // O table (opposite corner indices)
 int[] W = new int [3*maxnt];               // mid-edge vertex indices for subdivision (associated with corner opposite to edge)
 int[] Tc = new int[3*maxnt];               // corner type

int t (int c) {int r=int(c/3); return(r);};            int t() {return t(c);}
int n (int c) {int r=3*int(c/3)+(c+1)%3; return(r);};  int n() {return n(c);}
int p (int c) {int r=3*int(c/3)+(c+2)%3; return(r);};  int p() {return p(c);}
int v (int c) {return(V[c]);};                         int v() {return v(c);}
int o (int c) {return(O[c]);};                         int o() {return o(c);}
int l (int c) {return(o(n(c)));};                      int l() {return l(c);}
int r (int c) {return(o(p(c)));};                      int r() {return r(c);}
int s (int c) {return n(l(c));};
pt g (int c) {return(G[V[c]]);}; pt g() {return g(c);}            // shortcut to get the point of the vertex v(c) of corner c
vec N (int c) {return(Nv[V[c]]);}; vec N() {return N(c);}            // shortcut to get the point of the vertex v(c) of corner c
int w (int c) {return(W[c]);};               // temporary indices to mid-edge vertices associated with corners during subdivision
boolean nb(int c) {return O[c]!=-1 ;};  boolean nb() {return nb(c);}     // not border
boolean b(int c) {return O[c]==-1 ;};              // border: returns true if corner has no opposite
boolean vis(int c) {return visible[t(c)]; };   // true if tiangle of c is visible
  
void previous() {c=p();};
void next() {c=n();};
void opposite() {if(nb()) {c=o(c);};};
void left() {next(); opposite();};
void right() {previous(); opposite();};
void back() {opposite();};
void turn() {left(); next(); };

void writeCorner (int c) {println("c="+c+", n="+n(c)+", p="+p(c)+", o="+o(c)+", v="+v(c)+", t="+t(c)+"."+", nt="+nt+", nv="+nv ); }; 
void writeCorner () {writeCorner (c);}
void writeCorners () {for (int c=0; c<nc; c++) {println("T["+c+"]="+t(c)+", visible="+visible[t(c)]+", v="+v(c)+",  o="+o(c));};}

// ============================================= O TABLE CONSTRUCTION =========================================
void computeOnaive() {                         // sets the O table from the V table, assumes consistent orientation of triangles
  for (int i=0; i<3*nt; i++) {O[i]=-1;};  // init O table to -1: has no opposite (i.e. is a border corner)
  for (int i=0; i<3*nt; i++) {  for (int j=i+1; j<3*nt; j++) {       // for each corner i, for each other corner j
      if( (v(n(i))==v(p(j))) && (v(p(i))==v(n(j))) ) {O[i]=j; O[j]=i;};};}; // make i and j opposite if they match         
  };

void computeO() { 
  int nIC [] = new int [maxnv];                   // number of incident corners
  int maxValence=0;
  for (int c=0; c<nc; c++) {O[c]=-1;};  // init O table to -1: has no opposite (i.e. is a border corner)
  for (int v=0; v<nv; v++) {nIC[v]=0; };
  for (int c=0; c<nc; c++) {nIC[v(c)]++;}
  for (int v=0; v<nv; v++) {if(nIC[v]>maxValence) {maxValence=nIC[v]; };};
  println(" Max valence = "+maxValence+". ");
  int IC [][] = new int [maxnv][maxValence];                   // incident corners
  for (int v=0; v<nv; v++) {nIC[v]=0; };
  for (int c=0; c<nc; c++) {IC[v(c)][nIC[v(c)]++]=c;}
  for (int c=0; c<nc; c++) {
    for (int i=0; i<nIC[v(p(c))]; i++) {
      int a = IC[v(p(c))][i];
      for (int j=0; j<nIC[v(n(c))]; j++) {
         int b = IC[v(n(c))][j];
         if ((b==n(a))&&(c!=n(b))) {O[c]=n(b); O[n(b)]=c; };
         };
      };
    };
  }

// ============================================================= ARCHIVAL ============================================================
 boolean flipOrientation=false;            // if set, save will flip all triangles

void saveMesh() {
  String [] inppts = new String [nv+1+nt+1];
  int s=0;
  inppts[s++]=str(nv);
  for (int i=0; i<nv; i++) {inppts[s++]=str(G[i].x)+","+str(G[i].y)+","+str(G[i].z);};
  inppts[s++]=str(nt);
  if (flipOrientation) {for (int i=0; i<nt; i++) {inppts[s++]=str(V[3*i])+","+str(V[3*i+2])+","+str(V[3*i+1]);};}
    else {for (int i=0; i<nt; i++) {inppts[s++]=str(V[3*i])+","+str(V[3*i+1])+","+str(V[3*i+2]);};};
  saveStrings("mesh.vts",inppts);  println("saved on file");
  };

void loadMesh() {
  println("loading fn["+fni+"]: "+fn[fni]); 
  String [] ss = loadStrings(fn[fni]);
  String subpts;
  int s=0;   int comma1, comma2;   float x, y, z;   int a, b, c;
  nv = int(ss[s++]);
    print("nv="+nv);
    for(int k=0; k<nv; k++) {int i=k+s; 
      comma1=ss[i].indexOf(',');   
      x=float(ss[i].substring(0, comma1));
      String rest = ss[i].substring(comma1+1, ss[i].length());
      comma2=rest.indexOf(',');    y=float(rest.substring(0, comma2)); z=float(rest.substring(comma2+1, rest.length()));
      G[k].setTo(x,y,z);
    };
  s=nv+1;
  nt = int(ss[s]); nc=3*nt;
  println(", nt="+nt);
  s++;
  for(int k=0; k<nt; k++) {int i=k+s;
      comma1=ss[i].indexOf(',');   a=int(ss[i].substring(0, comma1));  
      String rest = ss[i].substring(comma1+1, ss[i].length()); comma2=rest.indexOf(',');  
      b=int(rest.substring(0, comma2)); c=int(rest.substring(comma2+1, rest.length()));
      V[3*k]=a;  V[3*k+1]=b;  V[3*k+2]=c;
    }
  };
  
//Read the standard obj file for vertices and triangles  
void readObj(String fileName)
{
  try
  {
    //TMesh ret = new TMesh();
    //Property<Vector3f> normals = new Property<Vector3f>("normal");
    println("Filename is " +fileName);
    
    BufferedReader bfr = new BufferedReader(new FileReader(fileName));
    boolean readVT = false;
    boolean readVN = false;
   
    String line;
    int nIndex = 0;
    
    while((line = bfr.readLine()) != null)
    {
    if(line.startsWith("#"))
    continue;
    /*else if(line.startsWith("vt"))
    readVT = true;
    else if(line.startsWith("vn"))
    {
    if(!readVN)
    ret.addVertexProperty(normals);
    readVN = true;
    String rest = line.substring(3);
    StringTokenizer st = new StringTokenizer(rest, " ");
    float x = Float.parseFloat(st.nextToken());
    float y = Float.parseFloat(st.nextToken());
    float z = Float.parseFloat(st.nextToken());

    normals.set(nIndex, new Vector3f(x, y, z));
    nIndex++;
    }*/
    else if(line.startsWith("v "))
    {
      String[] vertices = line.split(" ");
      addVertex(Float.parseFloat(vertices[1]), Float.parseFloat(vertices[2]), Float.parseFloat(vertices[3]));
    }
    else if(line.startsWith("f "))
    {
      String[] triangles = line.replace('/', ' ').split(" ");
      addTriangle(Integer.parseInt(triangles[1])-1, Integer.parseInt(triangles[4])-1, Integer.parseInt(triangles[7])-1);

    }
    
    }
      //if(nIndex != ret.vertices().size())
      //ret.computeVertexNormals();
      //System.out.println("Num Normals: " + nIndex);
      //System.out.println("Num Vertices: " + ret.vertices().size());
      //return(ret);
      nc = nt*3;
      saveMesh();
  }
  catch(IOException e)
  {
    e.printStackTrace();
  }
  //return(null);
  println(" nv = "+ nv +"nt = "+nt);
}

//  ==========================================================  SIMPLIFICATION ===========================================

void flip() {flip(c);}
void flip(int c) {      // flip edge opposite to corner c
    V[n(o(c))]=v(c); V[n(c)]=v(o(c));
    int co=o(c); O[co]=r(c); O[r(c)]=co; O[c]=r(co); O[r(co)]=c; O[p(c)]=p(co); O[p(co)]=p(c);
  }
  
void flipWhenLonger() {  for (int c=0; c<3*nt; c++) {
  if (nb(c)) {if (g(n(c)).disTo(g(p(c)))>g(c).disTo(g(o(c)))) {flip(c);}; }; };
  } 


void collapse() {collapse(c);}
void collapse(int c) {                                   // collapse edge opposite to corner c, does not check anything !!! assumes manifold
   int b=n(c), oc=o(c), vpc=v(p(c));
   visible[t(c)]=false; visible[t(oc)]=false;
   for (int a=b; a!=p(oc); a=n(l(a))) V[a]=vpc;
   O[l(c)]=r(c); O[r(c)]=l(c);     O[l(oc)]=r(oc); O[r(oc)]=l(oc); 
  }

//  ==========================================================  HOLES ===========================================
pt holeCenter = new pt (0,0,0);
vec holeNormal = new vec(0,0,1);
pt  centerOfHole() {pt C=new pt(0,0,0); int nb=0; for (int i=0; i<nc; i++) {if (visible[t(i)]&&b(i)) {nb++; C.addPt(g(p(i)));}; }; C.mul(1./nb); return C;};         // draws all border edges
vec  normalOfHole(pt C) {vec N=new vec(0,0,0); for (int i=0; i<nc; i++) {if (visible[t(i)]&&b(i)) N.add(cross(C.vecTo(g(p(i))),C.vecTo(g(n(i))))); }; N.makeUnit(); return N;};         // draws all border edges
void excludeInvisibleTriangles () {for (int b=0; b<nc; b++) {if (!visible[t(o(b))]) {O[b]=-1;};};}
void hole() {holeCenter.setTo(centerOfHole()); holeNormal.setTo(normalOfHole(holeCenter)); };
void fanHoles() {
  println("FANHOLES: nv="+nv +", nt="+nt +", nc="+nc );
  for (int t=0; t<nt; t++) {VisitedT[t]=false;};
  int lnt=nt;
  int L=0;
  for (int cc=0; cc<nc; cc++) {
   if (visible[t(cc)]&&(!VisitedT[t(cc)]) && (!nb(cc) )) {L++;
      print("<"); G[nv].setTo(0,0,0); int hl=fanHole(cc,L); G[nv].mul(1.0/float(hl)); nv++; println("> hl="+hl);
     };
   };
   for (int t=lnt; t<nt; t++) {visible[t]=true;};
   nc=3*nt;
   println("Filled "+L+" holes");
  }
  
int fanHole(int cc, int L) {
 int hl=0; int o=0;  int f=cc;
 print("."); VisitedT[t(f)]=true; hl++; o=3*nt; V[o]=nv; V[n(o)]=v(p(f)); V[p(o)]=v(n(f)); O[o]=f; O[f]=o; nt++; G[nv].addPt(g(p(f)));
 int lc=p(o); 
 f=n(f);  while(nb(f)) {f=n(o(f)); }; 
 while (f!=cc) {    
       print("."); VisitedT[t(f)]=true; hl++;  o=3*nt; V[o]=nv; V[n(o)]=v(p(f)); V[p(o)]=v(n(f)); O[o]=f; O[f]=o; nt++; G[nv].addPt(g(p(f)));
      O[n(o)]=lc; O[lc]=n(o); lc=p(o);
      f=n(f);  while(nb(f)&&(f!=cc)) {f=n(o(f)); }; 
      }; 
  O[lc]=n(o(cc)); O[n(o(cc))]=lc; c=lc;
  return(hl);
  }
  
void compactVO() {  
  println("COMPACT TRIANGLES: nv="+nv +", nt="+nt +", nc="+nc );
  int[] U = new int [nc];
  int lc=-1; for (int c=0; c<nc; c++) {if (visible[t(c)]) {U[c]=++lc; }; };
  for (int c=0; c<nc; c++) {if (nb(c)) {O[c]=U[o(c)];} else {O[c]=-1;}; };
  int lt=0;
  for (int t=0; t<nt; t++) {
    if (visible[t]) {
      V[3*lt]=V[3*t]; V[3*lt+1]=V[3*t+1]; V[3*lt+2]=V[3*t+2]; 
      O[3*lt]=O[3*t]; O[3*lt+1]=O[3*t+1]; O[3*lt+2]=O[3*t+2]; 
       visible[lt]=true; 
      lt++;
      };
    };
nt=lt; nc=3*nt;    
  println("      ...  NOW: nv="+nv +", nt="+nt +", nc="+nc );
  }

void compactV() {  
  println("COMPACT VERTICES: nv="+nv +", nt="+nt +", nc="+nc );
  int[] U = new int [nv];
  boolean[] deleted = new boolean [nv];
  for (int v=0; v<nv; v++) {deleted[v]=true;};
  for (int c=0; c<nc; c++) {deleted[v(c)]=false;};
  int lv=-1; for (int v=0; v<nv; v++) {if (!deleted[v]) {U[v]=++lv; }; };
  for (int c=0; c<nc; c++) {V[c]=U[v(c)]; };
  lv=0;
  for (int v=0; v<nv; v++) {
    if (!deleted[v]) {G[lv].setToPoint(G[v]);  deleted[lv]=false; 
      lv++;
      };
    };
 nv=lv;
 println("      ...  NOW: nv="+nv +", nt="+nt +", nc="+nc );
  }

// =========================================== GEODESIC MEASURES, DISTANCES =============================
 boolean[] P = new boolean [3*maxnt];       // marker of corners in a path to parent triangle
 int[] Distance = new int[maxnt];           // triangle markers for distance fields 
 int[] SMt = new int[maxnt];                // sum of triangle markers for isolation
 int prevc = 0;                             // previously selected corner

void computeDistance(int maxr) {
  int tc=0;
  int r=1;
  for(int i=0; i<nt; i++) {Mt[i]=0;};  Mt[t(c)]=1; tc++;
  for(int i=0; i<nv; i++) {Mv[i]=0;};
  while ((tc<nt)&&(r<=maxr)) {
      for(int i=0; i<nc; i++) {if ((Mv[v(i)]==0)&&(Mt[t(i)]==r)) {Mv[v(i)]=r;};};
     for(int i=0; i<nc; i++) {if ((Mt[t(i)]==0)&&(Mv[v(i)]==r)) {Mt[t(i)]=r+1; tc++;};};
     r++;
     };
  rings=r;
  }
  
void computeIsolation() {
  println("Starting isolation computation for "+nt+" triangles");
  for(int i=0; i<nt; i++) {SMt[i]=0;}; 
  for(c=0; c<nc; c+=3) {println("  triangle "+t(c)+"/"+nt); computeDistance(1000); for(int j=0; j<nt; j++) {SMt[j]+=Mt[j];}; };
  int L=SMt[0], H=SMt[0];  for(int i=0; i<nt; i++) { H=max(H,SMt[i]); L=min(L,SMt[i]);}; if (H==L) {H++;};
  c=0; for(int i=0; i<nt; i++) {Mt[i]=(SMt[i]-L)*255/(H-L); if(Mt[i]>Mt[t(c)]) {c=3*i;};}; rings=255;
  for(int i=0; i<nv; i++) {Mv[i]=0;};  for(int i=0; i<nc; i++) {Mv[v(i)]=max(Mv[v(i)],Mt[t(i)]);};
  println("finished isolation");
  }
  
void computePath() {
  for(int i=0; i<nt; i++) {Mt[i]=0;}; Mt[t(prevc)]=1; // Mt[0]=1;
  for(int i=0; i<nc; i++) {P[i]=false;};
  int r=1;
  boolean searching=true;
  while (searching) {
     for(int i=0; i<nc; i++) {
       if (searching&&(Mt[t(i)]==0)&&(o(i)!=-1)) {
         if(Mt[t(o(i))]==r) {
           Mt[t(i)]=r+1; 
           P[i]=true; 
           if(t(i)==t(c)){searching=false;};
           };
         };
       };
     r++;
     };
  for(int i=0; i<nt; i++) {Mt[i]=0;};
  rings=1;
  int b=c;
  int k=0;
   while (t(b)!=t(prevc)) {rings++;  
   if (P[b]) {b=o(b); print(".o");} else {if (P[p(b)]) {b=r(b);print(".r");} else {b=l(b);print(".l");};}; Mt[t(b)]=rings; };
  }

// ============================================================= SMOOTHING ============================================================
void computeValenceAndResetNormals() {      // caches valence of each vertex
  for (int i=0; i<nv; i++) {Nv[i].setTo(0,0,0); Valence[i]=0;};  // resets the valences to 0
  for (int i=0; i<nc; i++) {Valence[v(i)]++; };
  }

void computeLaplaceVectors() {  // computes the vertex normals as sums of the normal vectors of incident tirangles scaled by area/2
  computeValenceAndResetNormals();
  for (int i=0; i<3*nt; i++) {Nv[v(p(i))].add(g(p(i)).vecTo(g(n(i))));};
  for (int i=0; i<nv; i++) {Nv[i].div(Valence[i]);}; 
  };
  
void tuck(float s) {for (int i=0; i<nv; i++) {G[i].addScaledVec(s,Nv[i]);}; };  // displaces each vertex by a fraction s of its normal

// ============================================================= SUBDIVISION ============================================================
void splitEdges() {            // creates a new vertex for each edge and stores its ID in the W of the corner (and of its opposite if any)
  for (int i=0; i<3*nt; i++) {  // for each corner i
    if(b(i)) {G[nv]=midPt(g(n(i)),g(p(i))); W[i]=nv++;}
    else {if(i<o(i)) {G[nv]=midPt(g(n(i)),g(p(i))); W[o(i)]=nv; W[i]=nv++; }; };  // if this corner is the first to see the edge
    };
  };
  
void bulge() {              // tweaks the new mid-edge vertices according to the Butterfly mask
  for (int i=0; i<3*nt; i++) {
    if((nb(i))&&(i<o(i))) {    // no tweak for mid-vertices of border edges
     if (nb(p(i))&&nb(n(i))&&nb(p(o(i)))&&nb(n(o(i))))
      {G[W[i]].addScaledVec(0.25,midPt(midPt(g(l(i)),g(r(i))),midPt(g(l(o(i))),g(r(o(i))))).vecTo(midPt(g(i),g(o(i))))); };
      }; 
    };
  };
  
void splitTriangles() {    // splits each tirangle into 4
  for (int i=0; i<3*nt; i=i+3) {
    V[3*nt+i]=v(i); V[n(3*nt+i)]=w(p(i)); V[p(3*nt+i)]=w(n(i));
    V[6*nt+i]=v(n(i)); V[n(6*nt+i)]=w(i); V[p(6*nt+i)]=w(p(i));
    V[9*nt+i]=v(p(i)); V[n(9*nt+i)]=w(n(i)); V[p(9*nt+i)]=w(i);
    V[i]=w(i); V[n(i)]=w(n(i)); V[p(i)]=w(p(i));
    };
  nt=4*nt; nc=3*nt;
  };

// ============================================================= COMPRESSION ============================================================
// deleted

 } // ==== END OF MESH CLASS
  
float log2(float x) {float r=0; if (x>0.00001) { r=log(x) / log(2);} ; return(r);}
vec labelD=new vec(-10,-10, 2);           // offset vector for drawing labels
int maxr=1;


class IntQueue {
   private LinkedList q = new LinkedList();
  
   void add(int i) {q.add(new Integer(i)); }
   int removeFirst() { return ((Integer) q.removeFirst()).intValue(); }
   boolean isEmpty() {return q.isEmpty(); }
}

class IntSet {
   private HashSet s = new HashSet();
  
  void add(int i) {s.add(new Integer(i)); }
  boolean contains(int i) {return s.contains(new Integer(i)); } 
  
}
