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
    for(int i=0; i<this.Nx ;i++)
      for(int j=0; j<this.Ny ;j++){
        this.ForceFieldX[i][j] = 0.f;
        this.ForceFieldY[i][j] = 0.f;    
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
        grid[i][j].streaming(i, j, this);
    
    for(int i=0; i<Nx ;i++)
        for(int j=0; j<Ny ;j++)
          grid[i][j].collision(i, j, this);
  
    t++;
  }
}
