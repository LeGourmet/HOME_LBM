import java.util.Stack;

LBM simulation;

int GRID_SIZE_X = 300;
int GRID_SIZE_Y = 300;
int SCREEN_ZOOM = 3;

boolean paused = true;
COLOR_TYPE colorType = COLOR_TYPE.TYPE;

int nextFrame=0;

void settings(){
    size(GRID_SIZE_X*SCREEN_ZOOM,GRID_SIZE_Y*SCREEN_ZOOM,P2D);
}

void setup(){  
  simulation = new LBM(GRID_SIZE_X,GRID_SIZE_Y);
  
  simulation.setGlobalForceY(-0.0016f);
  
  boolean cave[][] = generateCaves(300, 300, 0.49, 10);
  
  for(int i=0; i<simulation.getNx() ;i++)
    for(int j=0; j<simulation.getNy() ;j++){
      
           if (j==0)                            simulation.setCell(i, j, new Cell(CELL_TYPE.SOLID, 1.f, 0.f, 0.f, 0.f, false)); 
      else if (j==simulation.getNy()-1)         simulation.setCell(i, j, new Cell(CELL_TYPE.SOLID, 1.f, 0.f, 0.f, 0.f, false));
      else if (j>270)                           simulation.setCell(i, j, new Cell(CELL_TYPE.LIQUID, 1.f, 0.f, 0.f, 1.f, false));
      else if (j>250 || j<50)                   simulation.setCell(i, j, new Cell(CELL_TYPE.GAS, 1.f, 0.f, 0.f, 0.f, false));
      else if (cave[i][j])                      simulation.setCell(i, j, new Cell(CELL_TYPE.SOLID, 1.f, 0.f, 0.f, 0.f, false)); 
      else                                      simulation.setCell(i, j, new Cell(CELL_TYPE.GAS, 1.f, 0.f, 0.f, 0.f, false));

      
      /*
      if (j==0)                                                          simulation.setCell(i, j, new Cell(CELL_TYPE.SOLID, 1.f, 0.f, 0.f, 0.f)); 
      else if (j==simulation.getNy()-1)                                  simulation.setCell(i, j, new Cell(CELL_TYPE.SOLID, 1.f, 0.f, 0.f, 0.f));
      else if (inSphere(i,j,GRID_SIZE_X/15,GRID_SIZE_X/2,GRID_SIZE_Y/2)) simulation.setCell(i, j, new Cell(CELL_TYPE.SOLID, 1.f, 0.f, 0.f, 0.f));
      else if (j<100 || 
               (j>GRID_SIZE_Y-GRID_SIZE_Y/10-GRID_SIZE_Y/30 && 
                j<GRID_SIZE_Y-GRID_SIZE_Y/30 &&
                i>GRID_SIZE_X/6 &&
                i<GRID_SIZE_X-GRID_SIZE_X/6))                            simulation.setCell(i, j, new Cell(CELL_TYPE.L, 1.f, 0.f, 0.f, 1.f));
      else                                                               simulation.setCell(i, j, new Cell(CELL_TYPE.G, 1.f, 0.f, 0.f, 0.f));
      */
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
          pixels[((simulation.getNy()-1-j)*SCREEN_ZOOM+a) * width + (i*SCREEN_ZOOM+b)] = col;
        }
      }
    }
  }
  updatePixels();
  
  float minP = Float.MAX_VALUE;
  float maxP = Float.MIN_VALUE;
  float minU = Float.MAX_VALUE;
  float maxU = Float.MIN_VALUE;
  for(int i=0; i<simulation.getNx() ;i++) {
    for(int j=0; j<simulation.getNy() ;j++) {
      if(simulation.getCell(i,j).getType()==CELL_TYPE.SOLID) continue;
      minP = min(minP, simulation.getCell(i,j).getDensity());
      maxP = max(maxP, simulation.getCell(i,j).getDensity());
      float tmpU = sqrt(sq(simulation.getCell(i,j).getVelocityX())+sq(simulation.getCell(i,j).getVelocityY()));
      minU = min(minU, tmpU);
      maxU = max(maxU, tmpU);
    }
  }
  
  text("Display : "+(colorType==COLOR_TYPE.PRESSURE ? "Density" : (colorType==COLOR_TYPE.VELOCITY ? "Velocity - Magnitude" : (colorType==COLOR_TYPE.TYPE ? "Type" : "Bubbles"))),10,15);
  text("Frame : "+simulation.getT(),10,30);
  text("P : ["+nf(minP,0,3)+", "+nf(maxP,0,3)+"]",10,45); 
  text("U : ["+nf(minU,0,3)+", "+nf(maxU,0,3)+"]",10,60);
  
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
}

void keyPressed(){
  if(keyCode==' ') paused = !paused;
  if(key==TAB) colorType = COLOR_TYPE.values()[(colorType.ordinal() + 1) % COLOR_TYPE.values().length];
}

void mousePressed(){
  int radius = (GRID_SIZE_X+GRID_SIZE_Y)/100;
  int mX = mouseX, mY = height-mouseY;
  
  for(int i=(mX/SCREEN_ZOOM-radius); i<(mX/SCREEN_ZOOM+radius) ;i++)
    for(int j=(mY/SCREEN_ZOOM-radius); j<(mY/SCREEN_ZOOM+radius) ;j++)
      if(i>=0 && i<simulation.getNx() && j>=0 && j<simulation.getNy() && inSphere(i,j,radius,mX/SCREEN_ZOOM,mY/SCREEN_ZOOM) && simulation.getCell(i,j).getType()!=CELL_TYPE.SOLID){
          simulation.getCell(i,j).pushType(CELL_TYPE.LIQUID);
          simulation.getCell(i,j).setEquilibrium(false);
      }
   
  for(int i=0; i<simulation.getNx() ;i++)
    for(int j=0; j<simulation.getNy() ;j++)
      simulation.getCell(i,j).init(i, j, simulation);
  simulation.updateBubbles();
}
