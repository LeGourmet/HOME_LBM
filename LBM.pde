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
  
  private Cell grid[][];
  
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
    
    this.grid = new Cell[this.Nx][this.Ny];
  }  
  
  // ------------------------------------------------------ GETTERS ------------------------------------------------------
  public int getNx() { return this.Nx; }
  public int getNy() { return this.Ny; }
  public int getT() { return this.t; }
  
  public float getForceX(int p_x, int p_y) { 
    if(p_x<0 || p_x>=this.Nx || p_y<0 || p_y>=this.Ny) return this.GlobalForceX; 
    return this.GlobalForceX+this.ForceFieldX[p_x][p_y]; 
  }
  
  public float getForceY(int p_x, int p_y) { 
    if(p_x<0 || p_x>=this.Nx || p_y<0 || p_y>=this.Ny) return this.GlobalForceY; 
    return this.GlobalForceY+this.ForceFieldY[p_x][p_y]; 
  }
  
  public Cell getCell(int p_x, int p_y) { 
    if(p_x<0 || p_x>=this.Nx || p_y<0 || p_y>=this.Ny) return null; 
    return this.grid[p_x][p_y];
  }
  
  public int getColor(int p_x, int p_y, COLOR_TYPE p_colorType) {
    if(p_x<0 || p_x>=this.Nx || p_y<0 || p_y>=this.Ny) return color(0);
    
    if (grid[p_x][p_y].type == CELL_TYPE.SOLID) return color(0);
    
    int palette[];
    float val;
    
    if(p_colorType == COLOR_TYPE.PRESSURE){
      palette = new int[]{color(68,1,84), color(59,82,139), color(33,145,140), color(94,201,98), color(253,231,37)};
      val = constrain(grid[p_x][p_y].p*0.5f, 0.f, 1.f);
    }else if(p_colorType == COLOR_TYPE.VELOCITY){
      palette = new int[]{color(70,70,219), color(0,255,91), color(0,128,0), color(255,255,0), color(255,96,0), color(107,0,0), color(223,77,77)};
      val = constrain(sqrt(grid[p_x][p_y].ux*grid[p_x][p_y].ux + grid[p_x][p_y].uy*grid[p_x][p_y].uy)/cs, 0.f, 1.f);
    }else{
      //palette = new int[]{color(40,40,180), color(150,150,200), color(221,221,221), color(200,150,150), color(180,40,40)};
      //val = (grid[p_x][p_y].phi>0.99f) ? 1.f : (grid[p_x][p_y].phi>0.01f) ? 0.f : 0.5f;
      palette = new int[]{color(40,40,180), color(180,40,40)};
      val = constrain(grid[p_x][p_y].phi,0,1);
  }
      
    float x = val*0.999f*(palette.length-1.f);
    int idx = (int)floor(x);
    return lerpColor(palette[idx],palette[idx+1],  x-(float)(idx));
  }
  
  // ------------------------------------------------------ SETTERS ------------------------------------------------------
  public void setGlobalForceX(float p_force) { this.GlobalForceX = p_force; }
  public void setGlobalForceY(float p_force) { this.GlobalForceY = p_force; }
  
  public void setForceFieldX(int p_x, int p_y, float p_force){
    if(p_x>=0 || p_x<this.Nx || p_y>=0 || p_y<this.Ny) this.ForceFieldX[p_x][p_y] = p_force;
  }
  
  public void setForceFieldY(int p_x, int p_y, float p_force){
    if(p_x>=0 || p_x<this.Nx || p_y>=0 || p_y<this.Ny) this.ForceFieldY[p_x][p_y] = p_force;
  }
  
  public void setCell(int p_x, int p_y, Cell p_cell) {
    if(p_x>=0 || p_x<this.Nx || p_y>=0 || p_y<this.Ny) this.grid[p_x][p_y] = p_cell;
  }
  
  // ----------------------------------------------------- FUNCTIONS -----------------------------------------------------
  void doTimeStep(){    
    for(int i=0; i<Nx ;i++)
      for(int j=0; j<Ny ;j++)
        grid[i][j].flowStreaming(i, j, this);
    
    for(int i=0; i<Nx ;i++)
        for(int j=0; j<Ny ;j++)
          grid[i][j].flowCollision(i, j, this);
    
    for(int a=0; a<1 ;a++){
      for(int i=0; i<Nx ;i++)
        for(int j=0; j<Ny ;j++)
          grid[i][j].phaseCollision(i, j, this);
        
      for(int i=0; i<Nx ;i++)
        for(int j=0; j<Ny ;j++)
          grid[i][j].phaseStreaming(i, j, this);
    }
    
    t++;
  }
}
