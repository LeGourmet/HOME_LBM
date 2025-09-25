public class LBM {
  // ----------------------------------------------------- ATTRIBUTS -----------------------------------------------------
  private int Nx;
  private int Ny;
  
  private float Fx;
  private float Fy;
  
  private int t;
  
  private float nu;
  private float rho;
  
  private Cell grid[][];
  
  // --------------------------------------------- DESTRUCTOR / CONSTRUCTOR ----------------------------------------------
  public LBM(int p_Nx, int p_Ny, float p_Fx, float p_Fy, float p_nu, float p_rho) {
    this.Nx = p_Nx;
    this.Ny = p_Ny;
    
    this.Fx = p_Fx;
    this.Fy = p_Fy;
    
    this.t  = 0;
    
    this.nu = p_nu;
    this.rho = p_rho;
    
    this.grid = new Cell[Nx][Ny];
  }  
  
  // ------------------------------------------------------ GETTERS ------------------------------------------------------
  public int getNx() { return this.Nx; }
  public int getNy() { return this.Ny; }
  public float getForceX() { return this.Fx; }
  public float getForceY() { return this.Fy; }
  public int getT() { return this.t; }
  public float getFluidNu() { return this.nu; }
  public float getFluidRho() { return this.rho; }
  public Cell getCell(int p_x, int p_y) { 
    if(p_x<0 || p_x>=Nx || p_y<0 || p_y>=Ny) return null; 
    return this.grid[p_x][p_y];
  }
  
  // ------------------------------------------------------ SETTERS ------------------------------------------------------
  public void setForceX(float p_force ) { this.Fx = p_force; }
  public void setForceY(float p_force) { this.Fy = p_force; }
  public void setCell(int p_x, int p_y, Cell p_cell) {
    if(p_x>=0 || p_x<Nx || p_y>=0 || p_y<Ny) this.grid[p_x][p_y] = p_cell;
  }
  
  // ----------------------------------------------------- FUNCTIONS -----------------------------------------------------
  void doTimeStep(){
    t++;
    
    for(int i=0; i<Nx ;i++)
      for(int j=0; j<Ny ;j++)
        grid[i][j].streaming(i, j, Nx, Ny, grid);
    
    for(int i=0; i<Nx ;i++)
        for(int j=0; j<Ny ;j++)
          grid[i][j].collision(nu, rho, Fx, Fy);
  }
}
