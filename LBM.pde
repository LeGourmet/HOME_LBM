public enum COLOR_TYPE { PRESSURE, VELOCITY, TYPE, BUBBLE };

public class LBM {
  // ----------------------------------------------------- ATTRIBUTS -----------------------------------------------------
  private int Nx;
  private int Ny;  
  private int B;
  private int t;

  private float GlobalForceX;
  private float GlobalForceY;
  private float ForceFieldX[][];
  private float ForceFieldY[][];
  
  private Cell grid[][];
  private Bubble bubbles[];
  
  // --------------------------------------------- DESTRUCTOR / CONSTRUCTOR ----------------------------------------------
  public LBM(int p_Nx, int p_Ny) {
    this.Nx = max(1,p_Nx);
    this.Ny = max(1,p_Ny);
    this.B = this.Nx*this.Ny+((this.Nx*this.Ny)%2);
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
    
    this.bubbles = new Bubble[this.B];
    for(int id=0; id<this.B ;id++) bubbles[id] = new Bubble();
  }  
  
  // ------------------------------------------------------ GETTERS ------------------------------------------------------
  public int getNx() { return this.Nx; }
  public int getNy() { return this.Ny; }
  public int getB() { return this.B; }
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
  
  public Bubble getBubble(int p_id) {
    if(p_id>=this.B || p_id<0) return null; 
    return this.bubbles[p_id];  
  }
  
  public int getColor(int p_x, int p_y, COLOR_TYPE p_colorType) {
    if(p_x<0 || p_x>=this.Nx || p_y<0 || p_y>=this.Ny) return color(0);
     
    if (p_colorType == COLOR_TYPE.TYPE) {
      if(grid[p_x][p_y].type == CELL_TYPE.GAS) return color(127);
      if(grid[p_x][p_y].type == CELL_TYPE.INTERFACE) return color(0,125,125);
      if(grid[p_x][p_y].type == CELL_TYPE.LIQUID) return color(0,0,200);
      return color(0);
    }
    
    if (grid[p_x][p_y].type == CELL_TYPE.SOLID )return color(0);
    
    if (p_colorType == COLOR_TYPE.BUBBLE) {
      randomSeed(simulation.getCell(p_x,p_y).getBubbleId());
      return color(random(255), random(255), random(255));
    }
        
    if (grid[p_x][p_y].type == CELL_TYPE.GAS)return color(30);
    
    int palette[];
    float val;
    
    if(p_colorType == COLOR_TYPE.PRESSURE){
      palette = new int[]{color(68,1,84), color(59,82,139), color(33,145,140), color(94,201,98), color(253,231,37)};
      val = constrain(grid[p_x][p_y].rho*0.5f, 0.f, 1.f);
    }else{
      palette = new int[]{color(70,70,219), color(0,255,91), color(0,128,0), color(255,255,0), color(255,96,0), color(107,0,0), color(223,77,77)};
      val = constrain(sqrt(sq(grid[p_x][p_y].ux) + sq(grid[p_x][p_y].uy))/cs, 0.f, 1.f);
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
  void init() {
    // ------------------------ Flag Init ------------------------
    for(int i=0; i<Nx ;i++)
      for(int j=0; j<Ny ;j++)
        grid[i][j].init(i, j, this);
    
    updateBubbles();
  }
  
  void doTimeStep(){    
    // ------------------- Flow Stream Collide -------------------
    for(int i=0; i<Nx ;i++)
      for(int j=0; j<Ny ;j++)
        grid[i][j].flowStreamingCollision(i, j, this);
    
    // ----------------------- Moments Swap ----------------------
    for(int i=0; i<Nx ;i++)
      for(int j=0; j<Ny ;j++)
        grid[i][j].flowSwapMoments(i, j, this);
        
    // ------------------------- Surface -------------------------
    for(int i=0; i<Nx ;i++)
      for(int j=0; j<Ny ;j++)
        grid[i][j].surface1(i, j, this);
    
    for(int i=0; i<Nx ;i++)
      for(int j=0; j<Ny ;j++)
        grid[i][j].surface2(i, j, this);
        
    for(int i=0; i<Nx ;i++)
      for(int j=0; j<Ny ;j++)
        grid[i][j].surface3(i, j, this);
        
    // ------------- Bubble Update Id and VolumeInit -------------
    updateBubbles();
    
    t++;
  }
  
  // ----------------------------------------------------- FUNCTIONS BUBBLES -----------------------------------------------------
  void updateBubbles() {
    // --------------------------- graph generation ---------------------------
    for(int i=0; i<Nx ;i++) {
      for(int j=0; j<Ny ;j++) {        
        if (!(grid[i][j].getType()==CELL_TYPE.INTERFACE || grid[i][j].getType()==CELL_TYPE.GAS)) {
          if(grid[i][j].getBubbleId()>=0) bubbles[grid[i][j].getBubbleId()].numberCells--;
          grid[i][j].setBubbleIdTmp(-1);
          continue;
        }
        int minId = i*Ny + j;
        for (int k=1; k<9 ;k+=1) {
          int idI = mod(i+D2Q9_cx[k],Nx);
          int idJ = mod(j+D2Q9_cy[k],Ny);
          int m = idI*Ny + idJ;
          if (grid[idI][idJ].getType()==CELL_TYPE.INTERFACE || grid[idI][idJ].getType()==CELL_TYPE.GAS) minId = min(minId,m);
        }
        grid[i][j].setBubbleIdTmp(minId);
      }
    }
    
    boolean unionContinue = true;
    while(unionContinue) {
      unionContinue = false;
      // --------------------------- graph reduction to local min id ---------------------------
      for(int i=0; i<Nx ;i++) {
        for(int j=0; j<Ny ;j++) {
          if (!(grid[i][j].getType()==CELL_TYPE.INTERFACE || grid[i][j].getType()==CELL_TYPE.GAS)) continue;
          int n = i*Ny + j;
          int m = grid[i][j].getBubbleIdTmp();
          if(m==n) continue;
          while(n!=m){
            n = m;
            m = grid[m/Ny][m%Ny].getBubbleIdTmp(); 
          }
          grid[i][j].setBubbleIdTmp(n);
        }
      }
      
      // --------------------------- graph local min id union ---------------------------
      for(int i=0; i<Nx ;i++) {
        for(int j=0; j<Ny ;j++) {
          if (!(grid[i][j].getType()==CELL_TYPE.INTERFACE || grid[i][j].getType()==CELL_TYPE.GAS)) continue;
          int n = grid[i][j].getBubbleIdTmp();
          for (int k=1; k<9 ;k+=1) {
            int idI = mod(i+D2Q9_cx[k],Nx);
            int idJ = mod(j+D2Q9_cy[k],Ny);
            int m = grid[idI][idJ].getBubbleIdTmp();
            if ((grid[idI][idJ].getType()==CELL_TYPE.INTERFACE || grid[idI][idJ].getType()==CELL_TYPE.GAS) && m<n) {
              grid[n/Ny][n%Ny].setBubbleIdTmp(m);
              unionContinue = true;
            }
          }
        }
      }    
    }
    
    for(int i=0; i<Nx ;i++) {
      for(int j=0; j<Ny ;j++) {
        int idBubbleN = grid[i][j].getBubbleIdTmp();
        int idBubbleN_old = grid[i][j].getBubbleId();
    
        if (!(grid[i][j].getType()==CELL_TYPE.INTERFACE || grid[i][j].getType()==CELL_TYPE.GAS)) {
          grid[i][j].setBubbleId(-1);
          continue;
        }
        
        idBubbleN = ((t%2)==0) ? idBubbleN-(idBubbleN%2) : idBubbleN+((idBubbleN+1)%2);
        
        bubbles[idBubbleN].volume += 1.f-grid[i][j].getPhi();
        if(idBubbleN_old >= 0) bubbles[idBubbleN].volumeInit += bubbles[idBubbleN_old].volumeInit/float(bubbles[idBubbleN_old].numberCells);
        bubbles[idBubbleN].numberCells++;
        
        grid[i][j].setBubbleId(idBubbleN);
      }
    }

    for(int b=0; b<B ;b++){
      if(bubbles[b].numberCells == 0) continue;
      
      if (t%2 == b%2){
        if (bubbles[b].volumeInit == 0.f) bubbles[b].volumeInit = bubbles[b].volume;
        bubbles[b].rho = cs * bubbles[b].volumeInit / bubbles[b].volume;
      } else bubbles[b].reset();
    }
    
  }
  
}
