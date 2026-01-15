LBM simulation;

int GRID_SIZE_X = 200;
int GRID_SIZE_Y = 200;
int SCREEN_ZOOM = 3;

boolean paused = true;
COLOR_TYPE colorType = COLOR_TYPE.TYPE;

void settings(){
    size(GRID_SIZE_X*SCREEN_ZOOM,GRID_SIZE_Y*SCREEN_ZOOM,P2D);
}

// warning dont use eq type
void setup(){  
  simulation = new LBM(GRID_SIZE_X,GRID_SIZE_Y);
  
  simulation.setGlobalForceY(0.0016f);
  
  for(int i=0; i<simulation.getNx() ;i++)
    for(int j=0; j<simulation.getNy() ;j++){
      if (j==simulation.getNy()-1 || (sqrt(sq(i-100)+sq(j-100))<30)) {
        simulation.setCell(i, j, new Cell(CELL_TYPE.SOLID, 1.f, 0.f, 0.f, 1.f));
      } else if(i>30 && i<170 && j<50 && j>10) {
        simulation.setCell(i, j, new Cell(CELL_TYPE.L, 1.f, 0.f, 0.f, 1.f));
      } else {
        simulation.setCell(i, j, new Cell(CELL_TYPE.G, 1.f, 0.f, 0.f, 0.f));
      }
    }
  
  simulation.init();
  
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
  
  if(!paused) simulation.doTimeStep();
  
  float minP = Float.MAX_VALUE;
  float maxP = Float.MIN_VALUE;
  float minU = Float.MAX_VALUE;
  float maxU = Float.MIN_VALUE;
  for(int i=0; i<simulation.getNx() ;i++) {
    for(int j=0; j<simulation.getNy() ;j++) {
      if(simulation.getCell(i,j).getType()==CELL_TYPE.SOLID || simulation.getCell(i,j).getType()==CELL_TYPE.EQUILIBRIUM) continue;
      minP = min(minP, simulation.getCell(i,j).getDensity());
      maxP = max(maxP, simulation.getCell(i,j).getDensity());
      float tmpU = sqrt(sq(simulation.getCell(i,j).getVelocityX())+sq(simulation.getCell(i,j).getVelocityY()));
      minU = min(minU, tmpU);
      maxU = max(maxU, tmpU);
    }
  }
  
  text("Display : "+(colorType==COLOR_TYPE.PRESSURE ? "Density" : (colorType==COLOR_TYPE.VELOCITY ? "Velocity - Magnitude" : "Type")),10,15);
  text("Frame : "+simulation.getT(),10,30);
  text("P : ["+nf(minP,0,3)+", "+nf(maxP,0,3)+"]",10,45); 
  text("U : ["+nf(minU,0,3)+", "+nf(maxU,0,3)+"]",10,60);
}

void keyPressed(){
  if(keyCode==' ') paused = !paused;
  if(key==TAB) colorType = COLOR_TYPE.values()[(colorType.ordinal() + 1) % COLOR_TYPE.values().length];
}

void mousePressed(){
}
