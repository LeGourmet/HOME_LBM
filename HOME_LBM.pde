LBM simulation;

int GRID_SIZE_X = 400;
int GRID_SIZE_Y = 100;
int SCREEN_ZOOM = 3;

boolean paused = true;
COLOR_TYPE colorType = COLOR_TYPE.TYPE;

//int nextFrame=0;

void settings(){
    size(GRID_SIZE_X*SCREEN_ZOOM,GRID_SIZE_Y*SCREEN_ZOOM,P2D);
}

void setup(){  
  simulation = new LBM(GRID_SIZE_X,GRID_SIZE_Y);
  
  //simulation.setGlobalForceX(-0.000016f);
  simulation.setGlobalForceX(0.00001f);
  //simulation.setGlobalForceX(0.0016f);
  
  for(int i=0; i<simulation.getNx() ;i++)
    for(int j=0; j<simulation.getNy() ;j++) {
           if (i==0)                                         { simulation.setCell(i, j, CELL_TYPE.SOLID, 1.f, 0.f, 0.f, 0.f); } 
      else if (i==simulation.getNx()-1)                      { simulation.setCell(i, j, CELL_TYPE.SOLID, 1.f, 0.f, 0.f, 0.f); }
      else                                                   { simulation.setCell(i, j, CELL_TYPE.FLUID, 1.f, 0.f, 0.f, float(i)>(float(2*GRID_SIZE_Y)+0.1f*float(GRID_SIZE_Y)*cos(2.f*PI*float(j)/float(GRID_SIZE_Y))) ? 1.f : 0.f); }
      //     if (i==0)                                         { simulation.setCell(i, j, CELL_TYPE.EQUILIBRIUM, 1.f, 0.f, 0.f, 0.f); } 
      //else if (i==simulation.getNx()-1)                      { simulation.setCell(i, j, CELL_TYPE.EQUILIBRIUM, 1.f, 0.f, 0.f, 0.f); } 
      //else if(inSphere(i,j,GRID_SIZE_X/30,50,GRID_SIZE_Y/2)) { simulation.setCell(i, j, CELL_TYPE.SOLID, 1.f, 0.f, 0.f, 0.f); } 
      //else                                                   { simulation.setCell(i, j, CELL_TYPE.FLUID, 1.f, 0.f, 0.f, 0.f); }
    }
    
  frameRate(120);
}

void draw(){
  loadPixels();
  for(int i=0; i<simulation.getNx() ;i++) {
    for(int j=0; j<simulation.getNy() ;j++) {
      color col = simulation.getColor(i,j,colorType);
      for(int a=0; a<SCREEN_ZOOM ;a++) {
        for(int b=0; b<SCREEN_ZOOM ;b++) {
          pixels[(j*SCREEN_ZOOM+a) * width + (i*SCREEN_ZOOM+b)] = col;
        }
      }
    }
  }
  updatePixels();
    
  /*if(!paused && simulation.getT() == nextFrame) {
    String name = "./video/image_";
    int id = nextFrame/10;
         if(id<10) name += "00000000";
    else if(id<100) name += "0000000";
    else if(id<1000) name += "000000";
    else if(id<10000) name += "00000";
    save(name+id+".png");
    nextFrame += 10;
  }*/
    
  if(!paused) simulation.doTimeStep();
  
  float minP = Float.MAX_VALUE;
  float maxP = Float.MIN_VALUE;
  float minU = Float.MAX_VALUE;
  float maxU = Float.MIN_VALUE;
  float minPhi = Float.MAX_VALUE;
  float maxPhi = Float.MIN_VALUE;
  float totalPhi = 0.f;
  for(int i=0; i<simulation.getNx() ;i++) {
    for(int j=0; j<simulation.getNy() ;j++) {
      if(simulation.getType(i,j)==CELL_TYPE.SOLID || simulation.getType(i,j)==CELL_TYPE.EQUILIBRIUM) continue;
      minP = min(minP, simulation.getPressure(i,j));
      maxP = max(maxP, simulation.getPressure(i,j));
      float tmpU = sqrt(sq(simulation.getVelocityX(i,j))+sq(simulation.getVelocityY(i,j)));
      minU = min(minU, tmpU);
      maxU = max(maxU, tmpU);
      minPhi = min(minPhi, simulation.getPhi(i,j));
      maxPhi = max(maxPhi, simulation.getPhi(i,j));
      totalPhi += simulation.getPhi(i,j);
    }
  }
  
  text("Display : "+(colorType==COLOR_TYPE.PRESSURE ? "Pressure" : (colorType==COLOR_TYPE.VELOCITY ? "Velocity - Magnitude" : "Type")),10,15);
  text("Frame : "+simulation.getT(),10,30);
  text("P : ["+nf(minP,0,3)+", "+nf(maxP,0,3)+"]",10,45); 
  text("U : ["+nf(minU,0,3)+", "+nf(maxU,0,3)+"]",10,60);
  text("Phi : ["+nf(minPhi,0,3)+", "+nf(maxPhi,0,3)+"]",10,75);
  text("Total Phi : "+nf(totalPhi,0,3),10,90);
}

void keyPressed(){
  if(keyCode==' ') paused = !paused;
  if(key==TAB) colorType = COLOR_TYPE.values()[(colorType.ordinal() + 1) % COLOR_TYPE.values().length];
}

void mousePressed(){
  int radius = (GRID_SIZE_X+GRID_SIZE_Y)/100;
  
  for(int i=(mouseX/SCREEN_ZOOM-radius); i<(mouseX/SCREEN_ZOOM+radius) ;i++)
    for(int j=(mouseY/SCREEN_ZOOM-radius); j<(mouseY/SCREEN_ZOOM+radius) ;j++) {
      if(i>=0 && i<simulation.getNx() && j>=0 && j<simulation.getNy() && inSphere(i,j,radius,mouseX/SCREEN_ZOOM,mouseY/SCREEN_ZOOM))
        if(simulation.getType(i,j)==CELL_TYPE.FLUID) 
          simulation.setPressure(i, j, 2.f);
  }
}
