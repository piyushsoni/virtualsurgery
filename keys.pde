boolean showHelpText=false;
  void showHelp() {
    fill(yellow,50); rect(0,0,height,height); pushMatrix(); translate(20,20); fill(0);
        text("                    MESH VIEWER written by Jarek Rossignac in June 2006, updated in February 2008",0,0); translate(0,20);
        translate(0,20);
        text("Click in this window to start. Press SPACE at any time to show/hide this help text",0,0); translate(0,20);
        text("'k:action' means press key 'k' to perform the 'action' ",0,0); translate(0,20);
        text("'*k:action' means click(&drag) the mouse while holding the key 'k' to perform 'action' ",0,0); translate(0,20);
        text("VIEW: *:rotate, *z:zoom, *v:pan, Z:home, j:jump, J:jumping(on/off)",0,0); translate(0,20);
        text("CORNER: *s:select, p:previous, n:next, o:opposite, l:left, r:right, t:turn",0,0); translate(0,20);
        text("VERTEX: V:vertices(show/hide), N:normals(show/hide), *m:move ",0,0); translate(0,20);
        text("EDGE: E:edges(show/hide), f:flip, F:flip(allLonger), C:findShortest c:collapse, L:collapse(shortest)",0,0); translate(0,20);
        text("TRIANGLE: T:triangles(show/hide), *d_elete(hide/show)",0,0); translate(0,20);
        text("MESH: R:refine, S:smooth, H:heal, M:mend(fillHoles)",0,0); translate(0,20);
        text("COMPRESS: b:begin, a:advance(oneTriangle), K:compress(mesh), B:show(triangleColors)",0,0); translate(0,20);
        text("DISTANCE: D:distance(how/hide) ,:smaller, .:larger, 0:zero, I:isolation, P:path(between selections) ",0,0); translate(0,20);
        text("FILES: 'g' next file. 'G' read model from file",0,0); translate(0,20);
        text("   ",0,0); translate(0,20);
        text("A:compressModel   (creates files mesh.con and mesh.geo in default directory)",0,0); translate(0,20);
        text("O:uncompressModel (loads files mesh.con and mesh.geo from default directory)",0,0); translate(0,20);
  
     popMatrix(); noFill();
        }
  void keys() {
  if (key==' ') {showHelpText=!showHelpText;};
  if (key==',') {maxr--; if (maxr<1) maxr=1; Current.computeDistance(maxr);};
  if (key=='.') {maxr++; Current.computeDistance(maxr);};
  if (key=='0') {maxr=0; Current.computeDistance(maxr);};
  if (key=='1') {Current.startT = Current.t(Current.c);}
  if (key=='2') {Current.setRay();}
  if (key=='a') {};
  if (key=='b') {};                 
  if (key=='c') {Current.collapse(); Current.left();};  // decimate one vertex by collapsing edge opposite to C
  if (key=='d') {C.setMark(); Current.hitTriangle();  Current.visible[Current.t(Current.c)]=!Current.visible[Current.t(Current.c)];};  // delete/undelete (show/hide) the triangle picked with mouse   
  if (key=='e') {Current.floodDelete();};
  if (key=='f') {Current.flip();};                 // flips edge opposite to c
  if (key=='g') {fni=(fni+1)%fniMax; println("Will read model "+fn[fni]);};
  if (key=='i') {};   
  if (key=='h') { }; 
  if (key=='j') {}
  if (key=='k') {};  
  if (key=='l') {Current.left(); if (jumps) C.jump(Current);};
 if (key=='m')  {Current.showColor = !Current.showColor;};   // also used in updateView to move vertex
  if (key=='n') {Current.next();};  
  if (key=='o') {Current.opposite(); if (jumps) C.jump(Current); };  
  if (key=='p') {Current.previous();};  
  if (key=='q') {colIndex = (colIndex + 1) % colorlist.length; }  
  if (key=='r') {Current.right(); if (jumps) C.jump(Current);}; 
  if (key=='s') {C.setMark(); Current.hitTriangle();/* C.F.setToPoint(mark); C.pullE(); C.pose(); */ Current.writeCorner(); };   // select current corner
  if (key=='t') {Current.turn();};  
  if (key=='u') {};   
 // if (key=='v') {};  // used in updateView 
  if (key=='w') {Current.writeCorner();};  
  if (key=='x') {};   
 // if (key=='y') {try {Current.readCLERS(Current.conFile);} catch (Exception e) {}};   
//  if (key=='z') {};   // used in updateView

  if (key=='A') {};  
  if (key=='B') {}; 
  if (key=='C') {Current.findShortestEdge(); if (jumps) C.jump(Current);};  // collpases shortest edge
  if (key=='D') {Current.showDistance=!Current.showDistance; if(Current.showDistance) {Current.computeDistance(maxr); println("showing distance");} else println("not showing distance");}; 
  if (key=='E') {Current.showEdges=!Current.showEdges; if(Current.showEdges) println("showing wires (edges) "); else println("not showing wires (edges)"); };  
  if (key=='F') {Current.flipWhenLonger(); println("flipped edges when this would shorten them");}; 
  if (key=='G') {println("loading fn["+fni+"]: "+fn[fni]); Current.loadMesh();  Current.init(); initView(Current); Current.fanHoles(); println("Loaded a model"); };  
  if (key=='H') {Current.excludeInvisibleTriangles();  Current.compactVO(); Current.compactV();  println("healed mesh by removing deleted triangles and compacting"); }; 
  if (key=='I') {Current.computeIsolation(); println("Computed isolation"); }; 
  if (key=='J') {jumps=!jumps;  if(jumps) println("will track corner c automatically"); else println("will not track corner c automatically"); }; 
  if (key=='K') {};  
  if (key=='L') {Current.findShortestEdge(); Current.collapse(); Current.left();}; 
//  if (key=='M') {Current.fanHoles();  println("mended holes by filling them with fans");   }; 
//  if (key=='N') {Current.showNormals=!Current.showNormals; if(Current.showNormals) println("showing normals"); else println("not showing normals"); };  
  if(key=='M') {--currentMeshId; if(currentMeshId<0) currentMeshId= MeshesInstantiaded -1; Current = Meshes[currentMeshId]; println("Current mesh id: " + currentMeshId); }
  if (key=='N') {++currentMeshId; if(currentMeshId>=MaxMeshes) currentMeshId=0; if(Meshes[currentMeshId] == null) {Meshes[currentMeshId] = new Mesh(currentMeshId); Meshes[currentMeshId].declare(); Meshes[currentMeshId].loadMesh(); Meshes[currentMeshId].init(); initView(Meshes[currentMeshId]); ++MeshesInstantiaded; } Current = Meshes[currentMeshId]; println("Current mesh id: " + currentMeshId); }
  if (key=='O') {}; 
  if (key=='P') {Current.computePath(); Current.showDistance=true; println("showing path"); println("not showing path"); };
  if (key=='Q') { }; 
  if (key=='R') {Current.splitEdges(); Current.bulge(); Current.splitTriangles(); Current.init(); print("Refined the mesh");};   // refine mesh
  if (key=='S') {Current.computeLaplaceVectors(); Current.tuck(0.6); Current.computeLaplaceVectors(); Current.tuck(-0.6); println("Smoothed the mesh");};  
  if (key=='T') {Current.showTriangles=!Current.showTriangles; if(Current.showTriangles) println("showing triangles"); else println("not showing triangles"); }; 
  if (key=='U') {}; 
  if (key=='V') {Current.showVertices=!Current.showVertices; if(Current.showVertices) println("showing vertices"); else println("not showing vertices"); };  
  if (key=='W') {Current.saveMesh();};  
  if (key=='X') {String S="mesh"+"-####.tif"; saveFrame(S); println("saved a picture");};   ;
  if (key=='Y') {}; 
  if (key=='Z') {C.F.setToPoint(Cbox); C.D=Rbox*2; C.U.setTo(0,1,0); C.E.setToPoint(C.F); C.E.addVec(new vec(0,0,1)); C.pullE(); C.pose(); println("resert view"); }; 

  if (keyCode==LEFT) {Current.left(); Current.right(); Current.left(); if (jumps) {C.jump(Current);};};
  if (keyCode==RIGHT) {Current.right(); Current.left(); Current.right(); if (jumps) {C.jump(Current);};};
  if (keyCode==DOWN) {Current.back(); Current.left(); Current.right(); Current.right(); Current.left(); Current.back(); if (jumps) {C.jump(Current);}; };
  if (keyCode==UP) {Current.left(); Current.right(); Current.right(); Current.left(); if (jumps) {C.jump(Current);};};
  }
  
void updateView() {
  if (keyCode==SHIFT) {C.Pan(); C.pullE(); };
  if (keyCode==CONTROL) {C.turn(); C.pullE(); };  
  if (key=='z') {C.pose(); C.zoom(); C.pullE(); };
   if (key=='1') {C.pose(); C.fly(1.0); C.pullE(); };  
   if (key=='2') {C.pose(); C.Pan(); C.pullE(); };
   if (key=='3') {C.pose(); C.fly(-1.0); C.pullE();  };
   if (key=='v') {C.pose(); C.Turn(); C.pullE(); };
   if (key=='5') {C.pan(); C.pullE(); }; 
   //if (key=='m') {Current.move();}; 
  }

pt Mouse = new pt(0,0,0);                 // current mouse position
float xr, yr = 0;                         // mouse coordinates relative to center of window
int px=0, py=0;                           // coordinats of mouse when it was last pressed
