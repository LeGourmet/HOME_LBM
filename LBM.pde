public enum CELL_TYPE { SOLID, FLUID, EQUILIBRIUM };
public enum COLOR_TYPE { PRESSURE, VELOCITY, TYPE };

public class LBM {
  // ----------------------------------------------------- ATTRIBUTS -----------------------------------------------------
  private int Nx;
  private int Ny;  
  private int t;

  private float GlobalForceX;
  private float GlobalForceY;
  private float ForceFieldX[][];
  private float ForceFieldY[][];
  
  private CELL_TYPE gridType[][];
  private CellFlow gridFlow[][];
  private CellPhase gridPhase[][];
  
  private static final int PhaseScale = 2;
  
  // --------------------------------------------- DESTRUCTOR / CONSTRUCTOR ----------------------------------------------
  public LBM(int p_Nx, int p_Ny) {
    this.Nx = max(1,p_Nx);
    this.Ny = max(1,p_Ny);
    this.t  = 0;
    
    this.GlobalForceX = 0.f;
    this.GlobalForceY = 0.f;
    this.ForceFieldX = new float[this.Nx][this.Ny];
    this.ForceFieldY = new float[this.Nx][this.Ny];
    for(int i=0; i<this.Nx ;i++){
      for(int j=0; j<this.Ny ;j++){
        this.ForceFieldX[i][j] = 0.f;
        this.ForceFieldY[i][j] = 0.f;    
      }
    }
    
    this.gridType = new CELL_TYPE[this.Nx][this.Ny];
    this.gridFlow = new CellFlow[this.Nx][this.Ny];
    this.gridPhase = new CellPhase[this.Nx*LBM.PhaseScale][this.Ny*LBM.PhaseScale];
  }  
  
  // ------------------------------------------------------ GETTERS ------------------------------------------------------
  public int getNx() { return this.Nx; }
  public int getNy() { return this.Ny; }
  public int getT() { return this.t; }
  public float getForceX(int p_x, int p_y) { return this.GlobalForceX+this.ForceFieldX[p_x][p_y]; }
  public float getForceY(int p_x, int p_y) { return this.GlobalForceY+this.ForceFieldY[p_x][p_y]; }
  public CELL_TYPE getType(int p_x, int p_y) { return this.gridType[p_x][p_y]; }
  public float getPressure(int p_x, int p_y) { return this.gridFlow[p_x][p_y].getPressure(); }
  public float getVelocityX(int p_x, int p_y) { return this.gridFlow[p_x][p_y].getVelocityX(); }
  public float getOldVelocityX(int p_x, int p_y) { return this.gridFlow[p_x][p_y].getOldVelocityX(); }
  public float getVelocityY(int p_x, int p_y) { return this.gridFlow[p_x][p_y].getVelocityY(); }
  public float getOldVelocityY(int p_x, int p_y) { return this.gridFlow[p_x][p_y].getOldVelocityY(); }
  public float getSxx(int p_x, int p_y){ return this.gridFlow[p_x][p_y].getSxx(); }
  public float getSyy(int p_x, int p_y){ return this.gridFlow[p_x][p_y].getSyy(); }
  public float getSxy(int p_x, int p_y){ return this.gridFlow[p_x][p_y].getSxy(); }
    
  // get Velocity => should trilinear for subgrid need
  
  public float getPhi(int p_x, int p_y){ 
    float phi = 0.f;
    for(int i=0; i<LBM.PhaseScale ;i++)
      for(int j=0; j<LBM.PhaseScale ;j++)  
        phi += this.gridPhase[p_x*LBM.PhaseScale+i][p_y*LBM.PhaseScale+j].getPhi();
    return phi/sq(float(LBM.PhaseScale)); 
  }
  
  public float getMicroPhi(int p_x, int p_y){ return gridPhase[p_x][p_y].getPhi(); }
  public float getHi(int p_x, int p_y, int p_i){ return gridPhase[p_x][p_y].getHi(p_i); }
    
  public int getColor(int p_x, int p_y, COLOR_TYPE p_colorType) {
    if(p_x<0 || p_x>=this.Nx || p_y<0 || p_y>=this.Ny) return color(0);
    
    CELL_TYPE type = getType(p_x, p_y);
    if (type == CELL_TYPE.SOLID) return color(0);
    
    int palette[];
    float val;
    
    if(p_colorType == COLOR_TYPE.PRESSURE){
      palette = new int[]{color(68,1,84), color(59,82,139), color(33,145,140), color(94,201,98), color(253,231,37)};
      val = constrain(getPressure(p_x,p_y)*0.5f, 0.f, 1.f);
    }else if(p_colorType == COLOR_TYPE.VELOCITY){
      palette = new int[]{color(70,70,219), color(0,255,91), color(0,128,0), color(255,255,0), color(255,96,0), color(107,0,0), color(223,77,77)};
      val = constrain(sqrt(sq(getVelocityX(p_x,p_y)) + sq(getVelocityY(p_x,p_y)))/cs, 0.f, 1.f);
    }else{
      //palette = new int[]{color(40,40,180), color(150,150,200), color(221,221,221), color(200,150,150), color(180,40,40)};
      //val = (grid[p_x][p_y].phi>0.99f) ? 1.f : (grid[p_x][p_y].phi>0.01f) ? 0.f : 0.5f;
      palette = new int[]{color(40,40,180), color(180,40,40)};
      val = constrain(getPhi(p_x,p_y),0,1);
  }
      
    float x = val*0.999f*(palette.length-1.f);
    int idx = (int)floor(x);
    return lerpColor(palette[idx],palette[idx+1],  x-(float)(idx));
  }
  
  // ------------------------------------------------------ SETTERS ------------------------------------------------------
  public void setGlobalForceX(float p_force) { this.GlobalForceX = p_force; }
  public void setGlobalForceY(float p_force) { this.GlobalForceY = p_force; }
  public void setForceFieldX(int p_x, int p_y, float p_force){ this.ForceFieldX[p_x][p_y] = p_force; }
  public void setForceFieldY(int p_x, int p_y, float p_force){ this.ForceFieldY[p_x][p_y] = p_force; }
  public void setType(int p_x, int p_y, CELL_TYPE p_type) { this.gridType[p_x][p_y] = p_type; }
  public void setPressure(int p_x, int p_y, float p_p) { this.gridFlow[p_x][p_y].setPressure(p_p); }
  public void setVelocityX(int p_x, int p_y, float p_ux) { this.gridFlow[p_x][p_y].setVelocityX(p_ux); }
  public void setVelocityY(int p_x, int p_y, float p_uy) { this.gridFlow[p_x][p_y].setVelocityY(p_uy); }
  
  public void setCell(int p_x, int p_y, CELL_TYPE p_type, float p_p, float p_ux, float p_uy, float p_phi) {
    this.gridType[p_x][p_y] = p_type;
    this.gridFlow[p_x][p_y] = new CellFlow(p_p, p_ux, p_uy);
    for(int i=0; i<LBM.PhaseScale ;i++)
      for(int j=0; j<LBM.PhaseScale ;j++)  
        this.gridPhase[p_x*LBM.PhaseScale+i][p_y*LBM.PhaseScale+j] = new CellPhase(p_phi,p_ux,p_uy);
  } 
  
  // ----------------------------------------------------- FUNCTIONS -----------------------------------------------------
  void doTimeStep(){
    for(int i=0; i<Nx ;i++)
      for(int j=0; j<Ny ;j++)
        gridFlow[i][j].streaming(i, j, this);
    
    for(int i=0; i<Nx ;i++)
      for(int j=0; j<Ny ;j++)
        gridFlow[i][j].collision(i, j, this);
    
    // first update with gi(x,t)
    for(int i=0; i<Nx ;i++)
      for(int j=0; j<Ny ;j++)
        for(int a=0; a<LBM.PhaseScale ;a++)
          for(int b=0; b<LBM.PhaseScale ;b++)
            gridPhase[i*LBM.PhaseScale+a][j*LBM.PhaseScale+b].collision(i*LBM.PhaseScale+a, j*LBM.PhaseScale+b, LBM.PhaseScale, getType(i,j), getOldVelocityX(i,j), getOldVelocityY(i,j), this);
    
    for(int i=0; i<Nx ;i++)
      for(int j=0; j<Ny ;j++)
        for(int a=0; a<LBM.PhaseScale ;a++)
          for(int b=0; b<LBM.PhaseScale ;b++)
            gridPhase[i*LBM.PhaseScale+a][j*LBM.PhaseScale+b].streaming(i*LBM.PhaseScale+a, j*LBM.PhaseScale+b, LBM.PhaseScale, getType(i,j), this);
    
    // should trilinear evaluate u !! 
    // second update with (gi(x,t)+gi(x,t+1))/2
    for(int i=0; i<Nx ;i++)
      for(int j=0; j<Ny ;j++)
        for(int a=0; a<LBM.PhaseScale ;a++)
          for(int b=0; b<LBM.PhaseScale ;b++)
            gridPhase[i*LBM.PhaseScale+a][j*LBM.PhaseScale+b].collision(i*LBM.PhaseScale+a, j*LBM.PhaseScale+b, LBM.PhaseScale, getType(i,j), (getOldVelocityX(i,j)+getVelocityX(i,j))*0.5f, (getOldVelocityY(i,j)+getVelocityY(i,j))*0.5f, this);
    
    // should trilinear evaluate u !! 
    for(int i=0; i<Nx ;i++)
      for(int j=0; j<Ny ;j++)
        for(int a=0; a<LBM.PhaseScale ;a++)
          for(int b=0; b<LBM.PhaseScale ;b++)
            gridPhase[i*LBM.PhaseScale+a][j*LBM.PhaseScale+b].streaming(i*LBM.PhaseScale+a, j*LBM.PhaseScale+b, LBM.PhaseScale, getType(i,j), this);
    
    t++;
  }
}
