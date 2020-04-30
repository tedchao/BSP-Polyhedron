void keyPressed() 
  {
//  if(key=='`') picking=true; 
  if(key=='n') P.nextCorner();
  if(key=='z') P.swing();


  if(key=='r') P.reorderPmcTree();
  if(key=='@') ;
  if(key=='?') scribeText=!scribeText;
  if(key=='!') snapPicture();
  if(key=='~') filming=!filming;
  if(key==']') showBalls=!showBalls;
  if(key=='f') {P.setPicekdLabel(key);}
  if(key=='t') {P.translateArrow(Of);}
  if(key=='s') semiTransparent=!semiTransparent;
  //if(key=='b') {P.setPicekdLabel(key);}
  if(key=='c') {P.setPicekdLabel(key);}
  if(key=='F') {P.addPt(Of,'f');}
  if(key=='S') {P.addPt(Of,'s');}
  //if(key=='B') {P.addPt(Of,'b');}
  if(key=='C') {P.addPt(Of,'c');}
  if(key=='m') {method=(method+1)%5;}
  if(key=='{') {showCurve=!showCurve;}
  if(key=='\\') {showKeys=!showKeys;}
  if(key=='}') {showPath=!showPath;}
  if(key=='|') {showCorrectedKeys=!showCorrectedKeys;}
  if(key=='=') {showTube=!showTube;}

  if(key=='3') {P.resetOnCircle(3,300); Q.copyFrom(P);}
  if(key=='4') {P.resetOnCircle(4,400); Q.copyFrom(P);}
  if(key=='5') {P.resetOnCircle(5,500); Q.copyFrom(P);}
  if(key=='^') track=!track;
  if(key=='q') Q.copyFrom(P);
  if(key=='p') P.copyFrom(Q);
  if(key==',') {level=max(level-1,0); f=0;}
  if(key=='.') {level++;f=0;}

  if(key=='e') {R.copyFrom(P); P.copyFrom(Q); Q.copyFrom(R);}
  if(key=='d') {P.set_pv_to_pp(); P.deletePicked();}
  if(key=='i') P.insertClosestProjection(Of); // Inserts new vertex in P that is the closeset projection of O
  //if(key=='W') {P.saveVertices("data/vertices"); P.saveEdges("data/edges"); P.savefaces("data/faces");}  // save data
  if(key=='L') {P.loadPts("data/pts"); Q.loadPts("data/pts2");}   // loads saved model
  //if(key=='w') P.savePts("data/pts");   // save vertices to pts
  if(key=='l') P.loadPts("data/pts"); 
  //if(key=='a') {animating=!animating; P.setFifo();}// toggle animation
  if(key=='^') showVecs=!showVecs;
  if(key=='#') exit();
  if(key=='=') {}
  
  if(key=='P') showPlane=!showPlane;
  if(key=='[') planeSize-=20;
  if(key==']') planeSize+=20;
  
  if(key==';') intersectSize++;
  if(key==':') {
    intersectSize--;
    if (intersectSize < 2) intersectSize = 2;
  }
  change=true;   // to save a frame for the movie when user pressed a key 
  }

void mouseWheel(MouseEvent event) 
  {
  dz -= event.getAmount(); 
  change=true;
  }
  
pt Pick = P();
void mousePressed() 
  {
  if (!keyPressed) {
    P.set_pv_to_pp(); println("picked vertex "+P.pp);
    P.reorderPmcTree();
  }
  if(keyPressed && key=='a') {P.addPt(Of);}  
  if(keyPressed && key=='b') {P.setPickedTo2(P(Pick),sphere);}
  }
  
void mouseMoved() 
  {
    if (keyPressed && key==' ') {rx-=PI*(mouseY-pmouseY)/height; ry+=PI*(mouseX-pmouseX)/width;};
    if (keyPressed && key=='`') dz+=(float)(mouseY-pmouseY); // approach view (same as wheel)
  }
  
void mouseDragged() 
  {
    if (!keyPressed) {
      P.translateArrow(Of);
      P.rotateArrow(Of);
    }
  }

// **** Header, footer, help text on canvas
void displayHeader()  // Displays title and authors face on screen
    {
    scribeHeader(title,0);
    scribeHeader("Number of planes: " + P.numberOfPlanes(),1);
    scribeHeader("Number of vertices: " + P.numberOfVertices(),2);
    scribeHeader("Number of edges: " + P.numberOfEdges(),3);
    scribeHeader("Number of faces: " + P.numberOfFaces(),4);
    scribeHeaderRight(name); 
    fill(white); image(myFace, width-myFace.width,25,myFace.width/2,myFace.height/2);
    fill(white); image(myFace2, width-myFace2.width/2,25,myFace2.width/2,myFace2.height/2);
    }
void displayFooter()  // Displays help text at the bottom
    {
    //scribeFooter(guide,1); 
    //scribeFooter(menu,0);
    scribeFooter(P.printBSPtree(),0);
    }

String title ="3D Polygon Generator", name ="Jacob Hanbeen Kim, Cheng-Kang Chao",
       menu="?:help, !:picture, ~:(start/stop)capture, space:rotate, `/wheel:closer, t/T:target, a:anim, #:quit",
       guide="click&drag:pick&slide on floor, xz/XZ:move/ALL, e:exchange, q/p:copy, l/L:load, w/W:write, m:subdivide method"; // user's guide
