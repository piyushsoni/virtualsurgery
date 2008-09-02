// Triangle mesh viewer + corner table + subdivision + smoothing + compression + simplificatioin + geodesics + isolation
// Written by Jarek Rossignac June 2006. Modified February 2008

import processing.opengl.*;                // load OpenGL
import javax.media.opengl.*; 
import javax.media.opengl.glu.*; 
import java.nio.*;
String [] fn=  {"horse.vts","femur.vts","ballJoint2m.vts", "bunny.vts","torus.vts","G3.vts.txt", "G3b.vts.txt", "tet.vts","fandisk.vts","squirrel.vts","venus.vts"};
int fni=0; int fniMax=fn.length;  

// ** SETUP **
void setup() { size(800, 800, OPENGL); setColors(); sphereDetail(6); //smooth();
  PFont font = loadFont("Courier-14.vlw"); textFont(font, 12);  // font for writing labels on screen
  Meshes[0] = new Mesh(0);
  Current = Meshes[0];
  Current.declare(); 
  //Current.loadMesh();
  Current.readObj("C:\\Documents and Settings\\psoni6\\Root\\Development\\Workspace\\Processing\\Cut Arteries\\VirtualSurgery\\data\\chopb.obj");
  Current.init(); 

  initView(Current);
  //while(!keyPressed){}
//  println("Exiting");
  //System.exit(0);  
  } 
 
// ** DRAW **
void draw() {
  mouse3D = fromMouse();
  background(white); 
  perspective(PI/2.0,width/height,1.0,6.0*Rbox); 
  if (showHelpText) {camera(); translate(-290,-290,0); scale(1.7,1.7,1.0); showHelp(); showColors();  return; };
  lights(); directionalLight(0,0,128,0,1,0); directionalLight(0,0,128,0,0,1);
  translate(float(height)/2, float(height)/2, 0.0);     // center view wrt window  
  if ((!keyPressed)&&(mousePressed)) {C.pan(); C.pullE(); };
  if ((keyPressed)&&(mousePressed)) {if (key=='p') Current.paint(colorlist[colIndex]); else updateView();}; 
  C1.track(C); C2.track(C1); C2.apply(); C2.setMark(); 

  Current.selectTriangle();
  for(int i=0; i<MeshesInstantiaded ;i++)
  {
    Meshes[i].show();
  }  
  //Current.show();  
  Current.walkToRay();

   
  }

//***      KEY ACTIONS (details in keys tab)
void keyPressed() { keys(); };
void mousePressed() {C.anchor(); C.pose();   if (keyPressed&&(key=='m')) {C.setMark(); Current.hitTriangle();};   };   // record where the cursor was when the mouse was pressed
void mouseReleased() {C.anchor(); C.pose(); };  // reset the view if any key was pressed when mouse was released 









