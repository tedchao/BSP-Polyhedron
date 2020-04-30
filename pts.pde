
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

class pts // class for manipulaitng and displaying pointclouds or polyloops in 3D 
{ 
  int maxnv = 16000;                 //  max number of vertices
  pt[] G = new pt [maxnv];           // geometry table (vertices)
  
  char[] L = new char [maxnv];       // labels of points
  vec [] LL = new vec[ maxnv];       // displacement vectors
  Boolean loop=true;                 // used to indicate closed loop 3D control polygons
  int pv =0,                         // picked vertex index,
    iv=0,                            //  insertion vertex index
    dv = 0,                          // dancer support foot index
    nv = 0,                          // number of vertices currently used in P
    pp=1;                            // index of picked vertex

  vec Up = V(0, 0, 1);               // Up vector
  
  
  // New Variables ==============================================================
  halfPlane[] BSP_tree = new halfPlane[5]; 
  ArrayList<vertices> verticesList = new ArrayList<vertices>();
  ArrayList<edge> edgeList = new ArrayList<edge>();
  ArrayList<face> faceList = new ArrayList<face>();
  
  // Corners
  int cc = 0;
  int cf = 0;
  
  HashMap<Integer, Integer[]> S = new HashMap<Integer, Integer[]>();
  
  // ============================================================================
  
  int[] PMC_tree = {0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
  
  int numberOfPlanes() {
    return floor(nv / 2);
  }
  int numberOfVertices() {
    return verticesList.size();
  }
  int numberOfEdges() {
    return edgeList.size();
  }
  int numberOfFaces() {
    return faceList.size();
  }
  
  int[] checkArraySize(int[] array, int index) {
    if (index >= array.length) {
      int[] temp = new int[array.length * 2];
      for(int i = 0; i < temp.length; i++) {
        if (i >= array.length) {
          temp[i] = -1;
        } else {
          temp[i] = array[i];
        }
      }
      array = temp;
    }
    
    return array;
  }
  
  pts() {
  }
  pts declare() 
  {
    for (int i=0; i<maxnv; i++) G[i]=P(); 
    for (int i=0; i<maxnv; i++) LL[i]=V(); 
    return this;
  }     // init all point objects
  pts empty() {
    nv=0; 
    pv=0; 
    return this;
  }                                 // resets P so that we can start adding points
  pts addPt(pt P, char c) { 
    G[nv].setTo(P); 
    pv=nv; 
    L[nv]=c; 
    nv++;  
    return this;
  }          // appends a new point at the end
  
  // appends a new point at the end
  // bound the length between A and B
  pts addPt(pt P) {
    P.z = 0;
    if (nv%2 == 0){
      G[nv].setTo(P);
    }
    else{
      vec AB = U(G[nv-1],P);
      G[nv].setTo(P(G[nv-1], V(r, AB)));
      addToPMC(floor(nv/2));
    }
    pv=nv;
    L[nv]='f';
    nv++;
    return this;
  }          
  pts addPt(float x, float y) { 
    G[nv].x=x; 
    G[nv].y=y; 
    pv=nv; 
    nv++; 
    return this;
  } // same byt from coordinates
  pts copyFrom(pts Q) {
    empty(); 
    nv=Q.nv; 
    for (int v=0; v<nv; v++) G[v]=P(Q.G[v]); 
    return this;
  } // set THIS as a clone of Q

  pts resetOnCircle(int k, float r)  // sets THIS to a polyloop with k points on a circle of radius r around origin
  {
    empty(); // resert P
    pt C = P(); // center of circle
    for (int i=0; i<k; i++) addPt(R(P(C, V(0, -r, 0)), 2.*PI*i/k, C)); // points on z=0 plane
    pv=0; // picked vertex ID is set to 0
    return this;
  } 
  // ********* PICK AND PROJECTIONS *******  
  int SETppToIDofVertexWithClosestScreenProjectionTo(pt M)  // sets pp to the index of the vertex that projects closest to the mouse 
  {
    pp=0; 
    for (int i=1; i<nv; i++) if (d(M, ToScreen(G[i]))<=d(M, ToScreen(G[pp]))) pp=i; 
    return pp;
  }
  pts showPicked(int r) {
    fill(brown,30);
    show(G[pv], r); 
    return this;
  }
  pt closestProjectionOf(pt M)    // Returns 3D point that is the closest to the projection but also CHANGES iv !!!!
  {
    pt C = P(G[0]); 
    float d=d(M, C);       
    for (int i=1; i<nv; i++) if (d(M, G[i])<=d) {
      iv=i; 
      C=P(G[i]); 
      d=d(M, C);
    }  
    for (int i=nv-1, j=0; j<nv; i=j++) { 
      pt A = G[i], B = G[j];
      if (projectsBetween(M, A, B) && disToLine(M, A, B)<d) {
        d=disToLine(M, A, B); 
        iv=i; 
        C=projectionOnLine(M, A, B);
      }
    } 
    return C;
  }

  // ********* MOVE, INSERT, DELETE *******  
  pts insertPt(pt P) { // inserts new vertex after vertex with ID iv
    for (int v=nv-1; v>iv; v--) {
      G[v+1].setTo(G[v]);  
      L[v+1]=L[v];
    }
    iv++; 
    G[iv].setTo(P);
    L[iv]='f';
    nv++; // increments vertex count
    return this;
  }
  pts insertClosestProjection(pt M) {  
    pt P = closestProjectionOf(M); // also sets iv
    insertPt(P);
    return this;
  }
  pts deletePicked() 
  {
    for (int i=pv; i<nv; i++) 
    {
      G[i].setTo(G[i+1]); 
      L[i]=L[i+1];
    }
    pv=max(0, pv-1); 
    nv--;  
    return this;
  }
  pts setPt(pt P, int i) { 
    G[i].setTo(P); 
    return this;
  }

  pts drawBalls(float r) {
    for (int v=0; v<nv; v++) show(G[v], r); 
    return this;
  }
  pts showPicked(float r) {
    show(G[pv], r); 
    return this;
  }
  
  pts showPts(float r){
    int halfnv = floor(nv/2);
    fill(dgreen);
    for (int v = 0; v < halfnv; v++){ 
      show(G[2*v], r/10);
      show(G[n(2*v)], r/20);
      if(semiTransparent){   // show transparent semisphere
        fill(yellow,100); show(G[2*v], sphere);
        fill(red,100); show(G[n(2*v)], sphere);
      }
    }
    
    return this;
  }
  
  pts drawClosedCurve(float r) 
  {
    fill(dgreen);
    for (int v=0; v<nv; v++) show(G[v], r*3);    
    fill(magenta);
    for (int v=0; v<nv-1; v++) stub(G[v], V(G[v], G[v+1]), r, r);  
    stub(G[nv-1], V(G[nv-1], G[0]), r, r);
    pushMatrix(); //translate(0,0,1); 
    scale(1, 1, 0.03);  
    fill(grey);
    for (int v=0; v<nv; v++) show(G[v], r*3);    
    for (int v=0; v<nv-1; v++) stub(G[v], V(G[v], G[v+1]), r, r);  
    stub(G[nv-1], V(G[nv-1], G[0]), r, r);
    popMatrix();
    return this;
  }
  pts set_pv_to_pp() {
    pv=pp; 
    return this;
  }
  pts movePicked(vec V) { 
    println("Updating pos");
    G[pv].add(V); 
    return this;
  }      // moves selected point (index p) by amount mouse moved recently
  pts setPickedTo(pt Q) { 
    G[pv].setTo(Q); 
    return this;
  }      // moves selected point (index p) by amount mouse moved recently
  
  pts setPickedTo2(pt Q, float r) { 
    pt B = P(G[pv]);
    float dist = d(Q, B);
    if (dist <= r) {
      G[pv].setTo(Q);
    }
    return this;
  }      // moves selected point (index p) by amount mouse moved recently
  
  pts moveAll(vec V) {
    for (int i=0; i<nv; i++) G[i].add(V); 
    return this;
  };   
  pt Picked() {
    return G[pv];
  } 
  pt Pt(int i) {
    if (0<=i && i<nv) return G[i]; 
    else return G[0];
  } 

  // ********* I/O FILE *******  
  void savefaces(String fn) 
  {
    String [] inppts = new String [floor(nv/2)+1];
    int s=0;
    inppts[s++]="Number of planes:" + str(floor(nv/2)) + ", " + "points" + ", " + "outward normal";
    for (int i=0; i<floor(nv/2); i++) {
      vec N = V(G[2*i], G[n(2*i)]);
      inppts[s++]="(" + str(G[i].x)+", "+str(G[i].y)+", "+str(G[i].z)+ ")" + ", "+ "(" + str(N.x)+", "+str(N.y)+", "+str(N.z)+ ")";
    }
    saveStrings(fn, inppts);
  };
  
  /*
  void saveVertices(String fn)
  {
    String [] v = new String [intersectCounter+1];
    int s=0;
    v[s++]= "ID" +"  "+"nv:"+str(intersectCounter) +"  "+ "0=convex, 1=concave, 2=on solid(mixed), 3=not on solid" ;
    for (int i=0; i<intersectCounter; i++) {
      v[s++]=str(i)+": "+str(intersectionPoints[i].x)+","+str(intersectionPoints[i].y)+","+str(intersectionPoints[i].z)+","+convexVertex[i];
    }
    saveStrings(fn, v);
  };
  
  void saveEdges(String fn)
  {
    String [] v = new String [edgeList.size()+1];
    int s=0;
    v[s++]= "ID  " + "ID_of_vertices" +"  "+"ne:"+str(edgeList.size()) +"  "+ "ID_of_planes   " + "0=convex, 1=concave, 2=on solid(mixed), 3=not on solid" ;
    for (int i=0; i<edgeList.size(); i++) {
      v[s++]=str(i)+": "+str(edges[i][0])+" "+str(edges[i][1]) + ";  " + str(edgesPlanes[i][0]) + " "+ str(edgesPlanes[i][1]) + ";  " +str(edgeConvex[i]);
    }
    saveStrings(fn, v);
  };*/

  void loadPts(String fn) 
  {
    println("loading: "+fn); 
    String [] ss = loadStrings(fn);
    String subpts;
    int s=0;   
    int comma, comma1, comma2;   
    float x, y;   
    int a, b, c;
    nv = int(ss[s++]); 
    print("nv="+nv+"\n");
    for (int k=0; k<nv; k++) 
    {
      int i=k+s; 
      //float [] xy = float(split(ss[i],",")); 
      String [] SS = split(ss[i], ","); 
      G[k].setTo(float(SS[0]), float(SS[1]), float(SS[2]));
      L[k]=SS[3].charAt(0);
    }
    pv=0;
  };

  // Dancer
  void setPicekdLabel(char c) {
    L[pp]=c;
  }



  void setFifo() 
  {
    _LookAtPt.reset(G[dv], 60);
  }              


  void next() {
    dv=n(dv);
  }
  int n(int v) {
    return (v+1)%nv;
  }
  int p(int v) {
    if (v==0) return nv-1; 
    else return v-1;
  }

  pts subdivideDemoInto(pts Q) 
  {
    Q.empty();
    for (int i=0; i<nv; i++)
    {
      Q.addPt(P(G[i])); 
      Q.addPt(P(G[i], G[n(i)])); 
      //...
    }
    return this;
  }  
  
  
  void displaySkater() 
  {
    if (showCurve) {
      fill(yellow); 
      for (int j=0; j<nv; j++) caplet(G[j], 6, G[n(j)], 6);
    }
    pt[] B = new pt [nv];           // geometry table (vertices)
    for (int j=0; j<nv; j++) B[j]=P(G[j], V(0, 0, 100));
    if (showPath) {
      fill(lime); 
      for (int j=0; j<nv; j++) caplet(B[j], 6, B[n(j)], 6);
    } 
    if (showKeys) {
      fill(cyan); 
      for (int j=0; j<nv; j+=4) arrow(B[j], G[j], 3);
    }

    if (animating) f=n(f);
    if (showSkater) 
    {
      // ....
    } else {
      fill(red); 
      arrow(B[f], G[f], 20);
    } //
  }
  
  //--------my implementations---------
  // show arrows and planes
  void showArrow(){
    int halfnv = floor(nv/2);
    fill(blue); 
    for (int i=0; i<halfnv; i++) arrow(G[2*i], G[n(2*i)], 6);
  }

  void showPlane(float d, color[] rainbow){
    int halfnv = floor(nv/2);
    
    for (int i=0; i<halfnv; i++){
      pt A = P(G[2*i]);
      pt B = P(G[n(2*i)]);
      vec Normal = U(A, B);
      
      // tool vector to find points on plane
      vec S = V(Normal.z, Normal.x, Normal.y);
      vec PlaneVector1 = U(cross(Normal, S));
      vec PlaneVector2 = U(cross(Normal, PlaneVector1));
      
      // find one point on a plane
      pt C = P(A, V(d, PlaneVector1));
      pt rC = P(A, V(-d, PlaneVector1));
      
      // shrink to unit vector and find a point near A
      pt D = P(A, V(d, PlaneVector2));
      pt rD = P(A, V(-d, PlaneVector2));
   
      fill(rainbow[i%7], 60 + 20*ceil(i/7)); show(C, D, rC, rD);
    }
  }
  
  void rotateArrow(pt Q) {
    int selected, center = -1;
    if (pv % 2 == 0) {
      selected = pv;
      center = pv + 1;
    } else {
      center = pv - 1;
      selected = pv;
    }    
    
    pt Selected = P(G[selected]);
    pt Center = P(G[center]);
    
    vec A = V(Center, Selected);
    vec B = V(Center, Q);
    
    //float angle = acos(dot(A,B) / (norm(A)*norm(B)));
    float angle = angle(A,B);
    
    float newX = Selected.x*cos(angle) - Selected.y*sin(angle);
    float newY = Selected.x*sin(angle) + Selected.y*cos(angle);
    
    pt C = P(Center, V(norm(A),U(Center, Q)));
    
    G[selected] = C;
  }
  
  void translateArrow(pt Q) {
    int aIdx, bIdx = -1;
    if (pv % 2 == 0) {
      aIdx = pv;
      bIdx = pv + 1;
    } else {
      aIdx = pv - 1;
      bIdx = pv;
    }    
    
    pt A = P(G[aIdx]);
    pt B = P(G[bIdx]);
    
    vec dir = U(A,Q);
    
    //G[aIdx].setTo(P(A, V(1,U(A,B))));
    //G[bIdx].setTo(P(B, V(1,U(A,B))));
    G[aIdx].setTo(P(A, V(1,dir)));
    G[bIdx].setTo(P(B, V(1,dir)));
  }
  
  boolean PMC(pt X, pt P, vec N) {
    float condition = (dot(V(P,X), N));
    if (condition < 0) return true;
    return false;
  }
  
  void addToPMC(int index) {
    if(index > 0) {
      PmcTree(0, index, G[index*2], G[0], V(G[0],G[1]));
      computeVertices();
    }
    println("Tree: ");
    for (int t : PMC_tree) { print(t+ ",");}
  }
  
  void PmcTree(int idx, int current, pt X, pt P, vec N) {
    PMC_tree = checkArraySize(PMC_tree, (2*idx)+2);
    if (PMC(X, P, N)) {
      int R = PMC_tree[(2*idx)+2];
      if (R != -1) {
        PmcTree((2*idx)+2, current, P(G[current*2]), P(G[R*2]), V(P(G[R*2]),P(G[R*2+1])));
      } else
      {
        PMC_tree[(2*idx)+2] = current;
      }
      
    }
    else {
      int L = PMC_tree[(2*idx)+1];
      if (L != -1) {
        PmcTree((2*idx)+1, current, P(G[current*2]), P(G[L*2]), V(P(G[L*2]),P(G[L*2+1])));
      } else
      {
        PMC_tree[(2*idx)+1] = current;
      }
    }
  }
  
  void reorderPmcTree() {
    for (int i=1; i < PMC_tree.length; i++) {
      PMC_tree[i] = -1;
    }
    for (int i=1; i < nv/2; i++) {
      PmcTree(0, i, G[i*2], G[0], V(G[0],G[1])); 
    }
    
    computeVertices();
    //println("Reordered Tree: ");
    //for (int t : PMC_tree) { print(t+ ",");}
  }
  
  String printBSPtree() {
    return "BSP Tree: " + printBSPtreeRecurse(0);
  }
  
  String printBSPtreeRecurse(int root) {
    if (root >= PMC_tree.length || PMC_tree[root] == -1) return " ";
    
    
    String tree = PMC_tree[root] + "";
    
    String left = printBSPtreeRecurse(2*root+1);
    String right = printBSPtreeRecurse(2*root+2);
    
    if (left == " " && right == " ")
    {
      return tree;
    }
    else
    {
      return "(" + left + ") " + tree + " (" + right + ")";
    }  
  }
  
  void computeVertices() {
    //println("Start Recursion...");
    verticesList = new ArrayList<vertices>();
    //computeVerticesRecurse(0, 0);
    computeV();
    testTheVertex();
    computeEdges();
    computeCorners();
  }
  
  void computeV() {
    int heightOfTree = floor(log(PMC_tree.length + 1) / log(2));
    int lastParentNodeIdx = floor(pow(2, heightOfTree) - 1);
    
    //println(lastParentNodeIdx);
    //println(PMC_tree.length);
    for (int i = 0; i < lastParentNodeIdx; i++) {
      int depth = floor(log(i + 1) / log(2)) + 1;
      for (int j = floor(pow(2, depth) - 1); j < PMC_tree.length - 1; j++) {
        for (int k = j + 1; k < PMC_tree.length; k++) {
          vertices vertex = findVertex(i, j, k);
          if (vertex != null) {
            verticesList.add(vertex);
          }
        }
      }
    }
  }
  
  vertices findVertex(int i, int j, int k) {
    if (PMC_tree[i] == -1 || PMC_tree[j] == -1 || PMC_tree[k] == -1) {
      return null;
    }
    pt D = P(G[PMC_tree[k]*2]);
    pt M = P(G[PMC_tree[j]*2]);
    pt K = P(G[PMC_tree[i]*2]);
    
    vec D_Normal = V(D, G[PMC_tree[k]*2 + 1]);
    vec M_Normal = V(M, G[PMC_tree[j]*2 + 1]);
    vec K_Normal = V(K, G[PMC_tree[i]*2 + 1]);
     
    // the intersecting line between plane 1 and 2
    vec T = cross(M_Normal, D_Normal);
    
    // Check if there is intersection point
    //println("Checking... " + dot(T, K_Normal));
    if (dot(T, K_Normal) != 0) {
      
      println("vertices:", "(", i, ", ", j, ", ", k, ")");
      
      // Solving linear system Nc = A
      float A1 = dot(V(K.x, K.y, K.z), K_Normal);
      float A2 = dot(V(M.x, M.y, M.z), M_Normal);
      float A3 = dot(V(D.x, D.y, D.z), D_Normal);
      
      // component of Normal matrix N
      float n11 = K_Normal.x, n12 = K_Normal.y, n13 = K_Normal.z;
      float n21 = M_Normal.x, n22 = M_Normal.y, n23 = M_Normal.z;
      float n31 = D_Normal.x, n32 = D_Normal.y, n33 = D_Normal.z;
      
      // determinant of Normal Matrix N
      float detN = n11*(n22*n33 - n23*n32) - n21*(n12*n33 - n13*n32) + n31*(n12*n23 - n13*n22);
      
      // Intersection point c
      float c_x = (1/detN) * (A1*(n22*n33 - n23*n32) - A2*(n33*n12 - n32*n13) + A3*(n23*n12 - n22*n13));
      float c_y = (1/detN) * (-A1*(n33*n21 - n31*n23) + A2*(n33*n11 - n31*n13) - A3*(n23*n11 - n21*n13));
      float c_z = (1/detN) * (A1*(n32*n21 - n31*n22) - A2*(n32*n11 - n31*n12) + A3*(n22*n11 - n21*n12));
      
      // get intersection point
      pt X = P(c_x, c_y, c_z);
      
      //halfPlane h0 = new halfPlane(PMC_tree[i]*2 + 1, K_Normal);
      halfPlane h0 = new halfPlane(PMC_tree[i], K_Normal);
      halfPlane h1 = new halfPlane(PMC_tree[j], M_Normal);
      halfPlane h2 = new halfPlane(PMC_tree[k], D_Normal);
      
      println("Indices: ", "(", PMC_tree[k], ", ", PMC_tree[j],", ", PMC_tree[i], ")");

      vertices V = new vertices(X, h0, h1, h2); 
      println("push x: ", "(", X.x, ", ", X.y, ", ", X.z, ")");
      
      return V;
    }
    
    return null;
  }
  
  /*
  Old code - delete later
  void computeVerticesRecurse(int index, int depth) {
    //println("Recursing..." + index + ", " + depth);
    if (index >= PMC_tree.length) return;
    //println("VAL: " + PMC_tree[index]);
    if(PMC_tree[index] == -1) return;
    
    //println("Depth: " + depth);
    //println("Index: " + index);
    if (depth > 1) {
      for (int i = 1; i < depth - 1; i++) {
        int parentIdx_1 = floor((index - pow(2, i) + 1) / pow(2,i));
        for (int j = 1; j < depth - i + 1; j++) {
          int parentIdx_2 = floor((parentIdx_1 - pow(2, j) + 1) / pow(2,j));
          
          //println("Computing..." + i + ", " + j + ".... (" + index + ", " + parentIdx_1 + ", " + parentIdx_2 + ")");
          
          pt D = P(G[PMC_tree[index]*2]);
          pt M = P(G[PMC_tree[parentIdx_1]*2]);
          pt K = P(G[PMC_tree[parentIdx_2]*2]);
          
          vec D_Normal = V(D, G[PMC_tree[index]*2 + 1]);
          vec M_Normal = V(M, G[PMC_tree[parentIdx_1]*2 + 1]);
          vec K_Normal = V(K, G[PMC_tree[parentIdx_2]*2 + 1]);
          
          // the intersecting line between plane 1 and 2
          vec T = cross(M_Normal, D_Normal);
          
          // Check if there is intersection point
          println("Checking... " + dot(T, K_Normal));
          if (dot(T, K_Normal) != 0) {
            
            // Solving linear system Nc = A
            float A1 = dot(V(K.x, K.y, K.z), K_Normal);
            float A2 = dot(V(M.x, M.y, M.z), M_Normal);
            float A3 = dot(V(D.x, D.y, D.z), D_Normal);
            
            // component of Normal matrix N
            float n11 = K_Normal.x, n12 = K_Normal.y, n13 = K_Normal.z;
            float n21 = M_Normal.x, n22 = M_Normal.y, n23 = M_Normal.z;
            float n31 = D_Normal.x, n32 = D_Normal.y, n33 = D_Normal.z;
            
            // determinant of Normal Matrix N
            float detN = n11*(n22*n33 - n23*n32) - n21*(n12*n33 - n13*n32) + n31*(n12*n23 - n13*n22);
            
            // Intersection point c
            float c_x = (1/detN) * (A1*(n22*n33 - n23*n32) - A2*(n33*n12 - n32*n13) + A3*(n23*n12 - n22*n13));
            float c_y = (1/detN) * (-A1*(n33*n21 - n31*n23) + A2*(n33*n11 - n31*n13) - A3*(n23*n11 - n21*n13));
            float c_z = (1/detN) * (A1*(n32*n21 - n31*n22) - A2*(n32*n11 - n31*n12) + A3*(n22*n11 - n21*n12));
            
            // get intersection point
            pt X = P(c_x, c_y, c_z);
            
            halfPlane h0 = new halfPlane(PMC_tree[parentIdx_2]*2 + 1, K_Normal);
            halfPlane h1 = new halfPlane(PMC_tree[parentIdx_1]*2 + 1, M_Normal);
            halfPlane h2 = new halfPlane(PMC_tree[index]*2 + 1, D_Normal);
            vertices V = new vertices(X, h0, h1, h2); 
            
            verticesList.add(V);
          }
        }
      }
    }
    
    depth = depth + 1;
    // left
    computeVerticesRecurse(2*index+1, depth);
    // right
    computeVerticesRecurse(2*index+2, depth);
  }
  */
  void testTheVertex()
  {
    float s = 2;
    
    for (vertices V : verticesList) {
      
      pt X = V.getVertex();
      
      halfPlane h0 = V.getH0();
      halfPlane h1 = V.getH1();
      halfPlane h2 = V.getH2();
      
      vec K = h0.getNormal();
      vec M = h1.getNormal();
      vec D = h2.getNormal();
      
      pt x1 = P(P(P(X, M), K), D); // X + M + K + D
      pt x2 = P(P(P(X, M), M(K)), D); // X + M - K + D
      pt x3 = P(P(P(X, M(M)), M(K)), D); // X - M - K + D
      pt x4 = P(P(P(X, M(M)), K), D); // X - M + K + D
      
      pt x5 = P(P(P(X, M), K), M(D)); // X + M + K - D
      pt x6 = P(P(P(X, M), M(K)), M(D)); // X + M - K - D
      pt x7 = P(P(P(X, M(M)), M(K)), M(D)); // X - M - K - D
      pt x8 = P(P(P(X, M(M)), K), M(D)); // X - M + K - D
      
      x1 = P(X, s, U(X, x1));
      x2 = P(X, s, U(X, x2));
      x3 = P(X, s, U(X, x3));
      x4 = P(X, s, U(X, x4));
      x5 = P(X, s, U(X, x5));
      x6 = P(X, s, U(X, x6));
      x7 = P(X, s, U(X, x7));
      x8 = P(X, s, U(X, x8));
      
      int K_idx = h0.getPlane();
      int M_idx = h1.getPlane();
      int D_idx = h2.getPlane();
      
      V.setBit(0, vertexPMCtest(x1, K_idx, M_idx, D_idx));
      V.setBit(1, vertexPMCtest(x2, K_idx, M_idx, D_idx));
      V.setBit(2, vertexPMCtest(x3, K_idx, M_idx, D_idx));
      V.setBit(3, vertexPMCtest(x4, K_idx, M_idx, D_idx));
      V.setBit(4, vertexPMCtest(x5, K_idx, M_idx, D_idx));
      V.setBit(5, vertexPMCtest(x6, K_idx, M_idx, D_idx));
      V.setBit(6, vertexPMCtest(x7, K_idx, M_idx, D_idx));
      V.setBit(7, vertexPMCtest(x8, K_idx, M_idx, D_idx));
      
      //println(V.getBits());
    }
  }
  
  int[] test_tree;
  boolean vertexPMCtest(pt X, int index1, int index2, int index3)
  {
    test_tree = new int[PMC_tree.length];
    for (int i = 0; i < test_tree.length; i++) { test_tree[i] = -1; }
    test_tree[0] = 0;
    
    for (int i = 0; i < floor(nv/2); i++) {
      //if (i != floor(index1/2) && i != floor(index2/2) && i != floor(index3/2))
      //{
      //  vertexPMCtestRecurse(0, i, G[i*2], G[0], V(G[0],G[1]));
      //}
      vertexPMCtestRecurse(0, i, G[i*2], G[0], V(G[0],G[1]));
    }
    
    boolean result = vertexPMCtestRecurse(0, floor(nv/2), X, G[0], V(G[0],G[1]));
    //println("RESULT: " + result);
    //println("Test Tree: ");
    //for (int t : test_tree) { print(t+ ",");}
    //println("");
    
    return result;
  }
  
  boolean vertexPMCtestRecurse(int idx, int current, pt X, pt P, vec N)
  {
    test_tree = checkArraySize(test_tree, (2*idx)+2);
    boolean result = false;
    /*if (test_tree[idx] == -1)
    {
      test_tree[idx] = current;
    }*/
    if (PMC(X, P, N)) {
      int R = test_tree[(2*idx)+2];
      if (R != -1) {
        result = vertexPMCtestRecurse((2*idx)+2, current, P(G[current*2]), P(G[R*2]), V(P(G[R*2]),P(G[R*2+1])));
      } else
      {
        test_tree[(2*idx)+2] = current;
        result = true;
      }
    }
    else {
      int L = test_tree[(2*idx)+1];
      if (L != -1) {
        result = vertexPMCtestRecurse((2*idx)+1, current, P(G[current*2]), P(G[L*2]), V(P(G[L*2]),P(G[L*2+1])));
      } else
      {
        test_tree[(2*idx)+1] = current;
        result = false;
      }
    }
    
    return result;
  }
  
  // draw intersection point with small square planes, edges, and witness points
  void drawIntersect(float intersectSize) 
  {
    int d = 20;
    
    for (vertices V : verticesList) {
      pt intersetionPoint = V.getVertex();
      if (intersetionPoint != null) {
        fill(red); show(intersetionPoint, intersectSize);
      }
      
      for (halfPlane h : V.getAllPlanes()) {
        vec Normal = h.getNormal();
        
        if (Normal != null)
        {
          // tool vector to find points on plane
          vec S = V(Normal.z, Normal.x, Normal.y);
          vec PlaneVector1 = U(cross(Normal, S));
          vec PlaneVector2 = U(cross(Normal, PlaneVector1));
          
          // find one point on a plane
          pt C = P(intersetionPoint, V(d, PlaneVector1));
          pt rC = P(intersetionPoint, V(-d, PlaneVector1));
          
          // shrink to unit vector and find a point near A
          pt D = P(intersetionPoint, V(d, PlaneVector2));
          pt rD = P(intersetionPoint, V(-d, PlaneVector2));
       
          fill(rainbow[h.getPlane()%7], 60 + 20*ceil(h.getPlane()/7)); show(C, D, rC, rD);
        }
      }
    }
    
    /*
    for (int i = 0; i < witnessPoints.length; i++) {
      for (int j = 0; j < witnessPoints[i].length; j++) {
        if (witnessPoints[i][j] != null) {
          fill(blue); show(witnessPoints[i][j], 1);
        }
      }
    }*/
    
    for (edge e : edgeList) {
      vertices v0 = e.getV0();
      vertices v1 = e.getV1();
      
      pt A = v0.getVertex();
      pt B = v1.getVertex();
      
      if (A != null && B != null)
      {
        if (e.IsConvex()) {
          fill(lime);
        } else {
          fill(pink);
        }
        caplet(A, 1, B, 1);
      }
    }
  }
  
  
  ArrayList<halfPlane> SameIntersectPlanes(halfPlane[] kmd1, halfPlane[] kmd2){
    
    ArrayList<halfPlane> SameList = new ArrayList<halfPlane>();
    for(int i = 0; i < kmd1.length; i++){
      for(int j = 0; j < kmd2.length; j++){
        int a = kmd1[i].getPlane();
        int b = kmd2[j].getPlane();
        println("a and b: ", a, b);
        if (a == b){
          SameList.add(kmd1[i]);
        }
      }
    }
    return SameList;
  }
  
  void computeEdges() {
    if (verticesList.size() > 1) {
      edgeList = new ArrayList<edge>();
      
      for (int i = 0; i < verticesList.size()-1; i++) {
        vertices v1 = verticesList.get(i);
        
        println("pull x: ", "(", v1.getVertex().x, ", ", v1.getVertex().y, ", ", v1.getVertex().z, ")");

        halfPlane h0_1 = v1.getH0();
        halfPlane h1_1 = v1.getH1();
        halfPlane h2_1 = v1.getH2();
        
        int k1 = h0_1.getPlane();
        int m1 = h1_1.getPlane();
        int d1 = h2_1.getPlane();
        
        println("plane 1: ", "(", k1, ", ", m1, ", ", d1, ")");
        
        halfPlane[] kmd1 = {h0_1, h1_1, h2_1}; 
                
        edge bestEdge = null;
        double bestDist = Double.POSITIVE_INFINITY;;
        
        for (int j = i+1; j < verticesList.size(); j++) {
          vertices v2 = verticesList.get(j);
        
          halfPlane h0_2 = v2.getH0();
          halfPlane h1_2 = v2.getH1();
          halfPlane h2_2 = v2.getH2();
          
          int k2 = h0_2.getPlane();
          int m2 = h1_2.getPlane();
          int d2 = h2_2.getPlane();
          
          println("plane 2: ", "(", k2, ", ", m2, ", ", d2, ")");
          
          halfPlane[] kmd2 = {h0_2, h1_2, h2_2};    
          
          pt A = v1.getVertex();
          pt B = v2.getVertex();
          double dist = d(A,B);
          ArrayList<halfPlane> Planes = SameIntersectPlanes(kmd1, kmd2); 
          
          
          if (k1 == k2 && m1 == m2 && d1 == d2) {   // same point, maybe no need
            continue;  
          }
          
          // need further check if we need to check the distance
          // it is easy to render all edges for visualization
          // but for data structure, maybe not?
          
          //else if (dist < bestDist) {
          //  bestDist = dist; 
            //if(k1 == k2 && m1 == m2) {
            //  boolean convexity = false;
            //  bestEdge = new edge(v1, v2, h0_1, h1_1, convexity);
            //}
            //else if (m1 == m2 && d1 == d2) {
            //  boolean convexity = false;
            //  bestEdge = new edge(v1, v2, h1_1, h2_1, convexity);
            //}
            
            //else if (k1 == k2 && d1 == d2) {
            //  boolean convexity = false;
            //  bestEdge = new edge(v1, v2, h0_1, h2_1, convexity);
            //}
           else if (Planes != null && Planes.size() == 2){
              boolean convexity = false;
              bestEdge = new edge(v1, v2, Planes.get(0), Planes.get(1), convexity);
              edgeList.add(bestEdge);
              v1.addNeighbor(v2);
              v2.addNeighbor(v1);
           // }
          }
        }
        
        //if (bestEdge != null) {
        //  edgeList.add(bestEdge);
        //}
      }
    }
  }
  
  void computeCorners() {
    faceList = new ArrayList<face>();
    for (int i = 0; i < verticesList.size()-1; i++) {
        vertices v = verticesList.get(i);
        findCycle(v, v, v, 0, new ArrayList<vertices>());
    }
    
    
    S = new HashMap<Integer, Integer[]>();
    for(int f=0; f < faceList.size(); f++) {
      face F = faceList.get(f);
      
      ArrayList<vertices> vert = F.getVertices();
      pt A = vert.get(0).getVertex();
      pt B = vert.get(1).getVertex();
      pt C = vert.get(2).getVertex();
      pt D = P(0,0,0);
      if (vert.size() == 4) {
        D = vert.get(3).getVertex();
      }
      if(!cw(A,B,C,D)) {
        println("reversed");
        Collections.reverse(vert);
        A = vert.get(0).getVertex();
        B = vert.get(1).getVertex();
        C = vert.get(2).getVertex();
        D = P(0,0,0);
        if (vert.size() == 4) {
          D = vert.get(3).getVertex();
        }
        println(cw(A,B,C,D));
      }
      for (int k = 0; k < vert.size(); k++) {
        vertices v1 = vert.get(k);
        int j = k+1;
        if (j == vert.size()) {
          j = 0;
        }
        vertices v2 = vert.get(j);
        
        for (int n=0; n < faceList.size(); n++) {
          if (f != n) {
            face N = faceList.get(n);
            if (N.getVertices().contains(v1) && N.getVertices().contains(v2)) {
              v1.setRightFace(n);
              println(f,k,n);
              if (S.containsKey(f)) {
                Integer[] rightFace = S.get(f);
                rightFace[k] = n;
              } else {
                Integer[] rightFace = new Integer[] {-1,-1,-1};
                if (vert.size() == 4) { rightFace = new Integer[] {-1,-1,-1,-1}; }
                
                rightFace[k] = n;
                S.put(f, rightFace);
              }
              break;
            }
          }
        }
      }
    }
    
    for (Integer key : S.keySet()) {
      println();
      println(key);
      Integer[] values = S.get(key);
      for(Integer val : values) {
        print(val);
        print(",");
      }
    }
  }
  
  ArrayList<vertices> findCycle(vertices start, vertices parent, vertices current, int depth, ArrayList<vertices> visited) {
    if (depth < 5) {
      for (vertices n : current.getNeighbors()) {
        if (!visited.contains(n) && !n.equalsTo(parent)) {
          ArrayList<vertices> vertices = (ArrayList<vertices>) visited.clone();
          vertices.add(n);
        
          if (start.equalsTo(n)) {
            if (vertices.size() == 4) {
              if (current.getNeighbors().contains(vertices.get(0))) {
                vertices.remove(parent);
              } else if (parent.getNeighbors().contains(start)) {
                vertices.remove(n);
              }
            }
            face newF = new face(vertices);
            
            boolean exist = false;
            for (face F : faceList) {
              if (F.equalsTo(newF)) {
                exist = true;
                break;
              }
            }
            
            if (!exist) {
              faceList.add(newF);
            }
          } else {
            findCycle(start, current, n, depth + 1, vertices);
          }
        }
      }
    }
    
    return visited;
  }
  
  void showFaces() {
    for (face F : faceList) {
      ArrayList<vertices> vertices = F.getVertices();
      fill(red, 70); 
      if (vertices.size() == 3) {
        show(vertices.get(0).getVertex(), vertices.get(1).getVertex(), vertices.get(2).getVertex());
      } else {
        show(vertices.get(0).getVertex(), vertices.get(1).getVertex(), vertices.get(2).getVertex(), vertices.get(3).getVertex());
      }
    }
  }
  
  void showCorners() {
    color[] colors = new color[] {red, green, blue, yellow};
    for (face F : faceList) {
      ArrayList<vertices> vertices = F.getVertices();
      for (int i=0; i < vertices.size(); i++) {
        pt v = vertices.get(i).getVertex();
        int j = i+1;
        if (j >= vertices.size()) {
          j = 0;
        }
        
        int k = i-1;
        if (k < 0) {
          k = vertices.size() - 1;
        }
        pt v_next = vertices.get(j).getVertex();
        pt v_prev = vertices.get(k).getVertex();
        fill(colors[i]); show(P(0.8,v,0.1,v_next,0.1,v_prev), 3);
      }
    }
  }
  
  void showCurrentCorners() {
    if (faceList.size() != 0) {
      ArrayList<vertices> vertices = faceList.get(cf).getVertices();
      pt v = vertices.get(cc).getVertex();
      
      int j = cc+1;
        if (j >= vertices.size()) {
        j = 0;
      }
      
      int k = cc-1;
      if (k < 0) {
        k = vertices.size() - 1;
      }
      
      pt v_next = vertices.get(j).getVertex();
      pt v_prev = vertices.get(k).getVertex();
      
      noFill(); pen(black, 1); show(P(0.8,v,0.1,v_next,0.1,v_prev), 5);
    }
  }
  
  void nextCorner() {
    cc++;
    if (cc >= faceList.get(cf).getVertices().size()) {
      cc = 0;
    }
  }
  
  void swing() {
    vertices v = faceList.get(cf).getVertices().get(cc);
    println("cur", cc,cf);
    cf = S.get(cf)[cc];
    if (cf != -1) {
      cc = faceList.get(cf).getVertices().indexOf(v);
      println("updated",cc, cf);
    }
  }
  
} // end of pts class

class halfPlane {
  
  int plane;
  vec normal;
  
  halfPlane() {}
  
  halfPlane(int p, vec n) {
    this.plane = p;
    this.normal = n;
  }
  
  void setPlane(int idx) {
    this.plane = idx;
  }
  
  void setNormal(vec n) {
    this.normal = n;
  }
  
  int getPlane() {
    return this.plane;
  }
  
  vec getNormal() {
    return this.normal;
  }
}
class vertices {
  pt vertex;
  halfPlane h0, h1, h2;
  boolean[] bits;
  ArrayList<vertices> neighbors;
  
  int rightFace;
  
  vertices() {}
  
  vertices(pt v, halfPlane p1, halfPlane p2, halfPlane p3) {
    this.vertex = v;
    this.h0 = p1;
    this.h1 = p2;
    this.h2 = p3;
    this.bits = new boolean[8];
    this.neighbors = new ArrayList<vertices>();
    this.rightFace = -1;
  }
  
  vertices(pt v, halfPlane p1, halfPlane p2, halfPlane p3, boolean[] b) {
    this.vertex = v;
    this.h0 = p1;
    this.h1 = p2;
    this.h2 = p3;
    this.bits = b;
    this.neighbors = new ArrayList<vertices>();
    this.rightFace = -1;
  }
  
  void setVertex(pt v) {
    this.vertex = v;
  }
  
  void setH0(halfPlane plane) {
    this.h0 = plane;
  }
  
  void setH1(halfPlane plane) {
    this.h1 = plane;
  }
  
  void setH2(halfPlane plane) {
    this.h2 = plane;
  }
  
  void setBits(boolean[] val) {
    this.bits = val;
  }
  
  void setBit(int idx, boolean val) {
    this.bits[idx] = val;
  }
  
  void setRightFace(int idx) {
    this.rightFace = idx;
  }
  
  void addNeighbor(vertices v) {
    if (!this.neighbors.contains(v)) {
      this.neighbors.add(v);
    }
  }
  
  pt getVertex() {
    return this.vertex;
  }
  
  halfPlane getH0() {
    return this.h0;
  }
  
  halfPlane getH1() {
    return this.h1;
  }
  
  halfPlane getH2() {
    return this.h2;
  }
  
  halfPlane[] getAllPlanes() {
    return new halfPlane[] {h0, h1, h2};
  }
  
  boolean[] getBits() {
    return this.bits;
  }
  
  boolean getBit(int idx) {
    return this.bits[idx];
  }
  
  int getRightFace() {
    return this.rightFace;
  }
  
  ArrayList<vertices> getNeighbors() {
    return this.neighbors;
  }
  
  boolean equalsTo(vertices v) {
    if (this.vertex == v.getVertex()) {
      return true;
    }
    
    return false;
  }
}

class edge {
  
  vertices v0;
  vertices v1;
  halfPlane h0;
  halfPlane h1;
  boolean isConvex;
  
  edge() {}
  
  edge(vertices v0, vertices v1, halfPlane h0, halfPlane h1) {
    this.v0 = v0;
    this.v1 = v1;
    this.h0 = h0;
    this.h1 = h1;
    this.isConvex = false;
  }
  
  edge(vertices v0, vertices v1, halfPlane h0, halfPlane h1, boolean isConvex) {
    this.v0 = v0;
    this.v1 = v1;
    this.h0 = h0;
    this.h1 = h1;
    this.isConvex = isConvex;
  }
  
  void setV0(vertices v) {
    this.v0 = v;
  }
  
  void setV1(vertices v) {
    this.v1 = v;
  }
  
  void setH0(halfPlane h) {
    this.h0 = h;
  }
  
  void setH1(halfPlane h) {
    this.h1 = h;
  }
  
  void testConvexity() {
    boolean convex = true;
    /*
    implement
    */
    this.isConvex = convex;
  }
  
  vertices getV0() {
    return this.v0;
  }
  
  vertices getV1() {
    return this.v1;
  }
  
  halfPlane getH0() {
    return this.h0;
  }
  
  halfPlane getH1() {
    return this.h1;
  }
  
  boolean IsConvex() {
    return this.isConvex;
  }
  
}

class face {
  
  ArrayList<vertices> faceVertices;
  
  face() {
    this.faceVertices = new ArrayList<vertices>();
  }
  
  face(ArrayList<vertices> faceVertices) {
    this.faceVertices = faceVertices;
  }
  
  void addVertex(vertices v, int idx) {
    this.faceVertices.add(idx, v);
  }
  
  ArrayList<vertices> getVertices() {
    return this.faceVertices;
  }
  
  boolean equalsTo(face F) {
    if (F.getVertices().size() != this.faceVertices.size()) {
      return false;
    }
    
    int counter = 0;
    for (vertices v : F.getVertices()) {
      for (vertices v2 : this.faceVertices) {
        if (v.equalsTo(v2)) {
          counter++;
        }
      }
    }
    
    if (counter == this.faceVertices.size()) {
      return true;
    }
    
    return false;
  }
}
