LBM simulation;

int GRID_SIZE_X = 300;
int GRID_SIZE_Y = 100;
int SCREEN_ZOOM = 3;

void settings(){
    size(GRID_SIZE_X*SCREEN_ZOOM,GRID_SIZE_Y*SCREEN_ZOOM,P2D);
}

void setup(){  
  simulation = new LBM(GRID_SIZE_X,GRID_SIZE_Y);
  
  simulation.setGlobalForceX(0.00075f);
  
  for(int i=0; i<simulation.getNx() ;i++)
    for(int j=0; j<simulation.getNy() ;j++) {
      if(i==0){
        simulation.setCell(i,j,new Cell(CELL_TYPE.EQUILIBRIUM, 1.f, 0.f, 0.f, 0.f));
      } else if (i==simulation.getNx()-1) {
        simulation.setCell(i,j,new Cell(CELL_TYPE.EQUILIBRIUM, 1.f, 0.f, 0.f, 0.f));  
      } else if(inSphere(i,j,GRID_SIZE_X/30,50,GRID_SIZE_Y/2)){
        simulation.setCell(i,j,new Cell(CELL_TYPE.SOLID, 1.f, 0.f, 0.f, 0.f));
      } else {
        simulation.setCell(i,j,new Cell(CELL_TYPE.FLUID, 1.f, 0.f, 0.f, 0.f));
      }
    }
    
    frameRate(120);
}

void draw(){
  loadPixels();
  for(int i=0; i<simulation.getNx() ;i++) {
    for(int j=0; j<simulation.getNy() ;j++) {
      color col = simulation.getCell(i,j).getColor(colorType);
      for(int a=0; a<SCREEN_ZOOM ;a++) {
        for(int b=0; b<SCREEN_ZOOM ;b++) {
          pixels[(j*SCREEN_ZOOM+a) * width + (i*SCREEN_ZOOM+b)] = col;
        }
      }
    
    }
  }
  updatePixels();
    
  if(!paused) simulation.doTimeStep();
  
  text("framerate : "+int(frameRate), 10, 15);
}

void keyPressed(){
  if(key=='p' || key=='P') paused = !paused;
  if(key==TAB) colorType = COLOR_TYPE.values()[(colorType.ordinal() + 1) % COLOR_TYPE.values().length];
}

void mousePressed(){
  int radius = (GRID_SIZE_X+GRID_SIZE_Y)/100;
  
  for(int i=(mouseX/SCREEN_ZOOM-radius); i<(mouseX/SCREEN_ZOOM+radius) ;i++)
    for(int j=(mouseY/SCREEN_ZOOM-radius); j<(mouseY/SCREEN_ZOOM+radius) ;j++) {
      if(i>=0 && i<simulation.getNx() && j>=0 && j<simulation.getNy() && inSphere(i,j,radius,mouseX/SCREEN_ZOOM,mouseY/SCREEN_ZOOM))
        if(simulation.getCell(i,j).getType()==CELL_TYPE.FLUID) 
          simulation.getCell(i,j).setPressure(2.f);
  }
}
