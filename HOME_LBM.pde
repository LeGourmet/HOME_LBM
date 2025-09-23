LBM simulation;

void setup(){
  size(500,200,P2D);
  
  simulation = new LBM(500,200,0.f,0.f,0.001f,1.f);
  
  for(int i=0; i<simulation.getNx() ;i++)
    for(int j=0; j<simulation.getNy() ;j++) {
      if(inSphere(i,j,20,50,100) || j<1 || j>=simulation.getNy()-1){
        simulation.setCell(i,j,new Cell(CELL_TYPE.SOLID, 0.f, 0.f, 1.f, 0.f, 0.f));
      }else if(i==0) {
        simulation.setCell(i,j,new Cell(CELL_TYPE.EQUILIBRIUM, 0.f, 0.0f, 1.f, 0.15f, 0.f));
      } else if(i==simulation.getNx()-1) {
        simulation.setCell(i,j,new Cell(CELL_TYPE.EQUILIBRIUM, 0.f, 0.f, 1.f, 0.f, 0.f));
      } else {
        simulation.setCell(i,j,new Cell(CELL_TYPE.FLUID, 0.f, 0.f, 1.f, 0.f, 0.f));
      }
    }
    
  noStroke();
}

void draw(){
  background(30);
  
  loadPixels();
  for(int i=0; i<simulation.getNx() ;i++) {
    for(int j=0; j<simulation.getNy() ;j++) {
      if(simulation.getCell(i,j) == null) continue;
      pixels[j * width + i] = simulation.getCell(i,j).getColor();
    }
  }
  updatePixels();
    
  simulation.doTimeStep();
}

void keyPressed(){
  if(key=='p' || key=='P') println("framerate :"+str(frameRate));
}

void mousePressed(){
  int radius = 15;
  
  for(int i=(mouseX-radius); i<(mouseX+radius) ;i++)
    for(int j=(mouseY-radius); j<(mouseY+radius) ;j++) {
      if(i>=0 && i<simulation.getNx() && j>=0 && j<simulation.getNy() && inSphere(i,j,radius,mouseX,mouseY))
        //simulation.getCell(i,j).setPressure(1.5f);
        simulation.getCell(i,j).setPressure(simulation.getCell(i,j).getPressure()+1.5f);  
  }
}
