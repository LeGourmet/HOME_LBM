LBM simulation;

void setup(){
  size(400,400,P2D);
  
  simulation = new LBM(400,400,0.f,0.001f);
  
  for(int i=0; i<simulation.getNx() ;i++)
    for(int j=0; j<simulation.getNy() ;j++) {
      if(inSphere(i,j,20,200,50) || i<=1 || i>=simulation.getNx()-2 ){
        simulation.setCell(i,j,new Cell(CELL_TYPE.SOLID, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f));
      } else if(inSphere(i,j,20,200,150)){
          simulation.setCell(i,j,new Cell(CELL_TYPE.FLUID, 0.f, 0.f, 1.f, 0.f, 0.f, 1.f));        
      } else if(j==0|| j==simulation.getNy()-1) {
        simulation.setCell(i,j,new Cell(CELL_TYPE.EQUILIBRIUM, 0.f, 0.0f, 1.f, 0.f, 0.f, 0.f));
      } else {
        simulation.setCell(i,j,new Cell(CELL_TYPE.FLUID, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f));
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
      //pixels[j * width + i] = simulation.getCell(i,j).getColorType();
      pixels[j * width + i] = simulation.getCell(i,j).getColorVelocity();
      //pixels[j * width + i] = simulation.getCell(i,j).getColor();
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
        simulation.getCell(i,j).setPressure(1.5f);
        //simulation.getCell(i,j).setPressure(simulation.getCell(i,j).getPressure()+1.5f);  
  }
}
