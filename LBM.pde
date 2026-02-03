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
    this.B = this.Nx*this.Ny;
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
      if(grid[p_x][p_y].type == CELL_TYPE.G) return color(127);
      if(grid[p_x][p_y].type == CELL_TYPE.I) return color(0,125,125);
      if(grid[p_x][p_y].type == CELL_TYPE.L) return color(0,0,200);
      return color(0);
    }
    
    if (grid[p_x][p_y].type == CELL_TYPE.SOLID )return color(0);
    
    if (p_colorType == COLOR_TYPE.BUBBLE) {
      randomSeed(simulation.getCell(p_x,p_y).getBubbleId());
      return color(random(255), random(255), random(255));
    }
        
    if (grid[p_x][p_y].type == CELL_TYPE.G )return color(30);
    
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
    
    // ------- Connect Component Labelling (=> flood fill) -------
    int currentLabel = 0;
    Stack<Integer> coords = new Stack<>();
    
    for(int i=0; i<Nx ;i++) {
      for(int j=0; j<Ny ;j++) {
        if(grid[i][j].getBubbleId() != -1 || (grid[i][j].getType() != CELL_TYPE.I && grid[i][j].getType() != CELL_TYPE.G)) continue;
        coords.clear();
        coords.push(j*Nx+i);
        
        while (!coords.empty()) {
            int coord = coords.pop();
            int coordX = coord%Nx;
            int coordY = coord/Nx;
            
            if(grid[coordX][coordY].getBubbleId() != -1) continue;
            
            grid[coordX][coordY].bubbleId = currentLabel;
            bubbles[currentLabel].numberCells++;
            bubbles[currentLabel].volumeInit += 1.f-grid[coordX][coordY].getPhi();
            
            for (int k=1; k<9 ;k++) {
              int offsetX = mod(coordX+D2Q9_cx[k],Nx);
              int offsetY = mod(coordY+D2Q9_cy[k],Ny);
              if(grid[offsetX][offsetY].getBubbleId() == -1 && (grid[offsetX][offsetY].getType() == CELL_TYPE.I || grid[offsetX][offsetY].getType() == CELL_TYPE.G)) coords.push(offsetY*Nx+offsetX);
            }
        }
        
        //println("Bubble :", currentLabel, " with", bubbles[currentLabel].numberCells, " cells, volume :", bubbles[currentLabel].volumeInit);
        
        currentLabel++;
      }
    }
    
    // ----------------------- Bubble Init -----------------------
    for(int id=0; id<B ;id++) {
      if(bubbles[id].numberCells==0) continue;
      bubbles[id].rho = cs2;
      bubbles[id].volume = bubbles[id].volumeInit;
    }
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
    
    for(int i=0; i<Nx ;i++) {
      for(int j=0; j<Ny ;j++) {
        int idB = grid[i][j].getBubbleId();
        if(!(idB>=0 && grid[i][j].getType()==CELL_TYPE.L)) continue; // if(!split) continue;
                
        bubbles[idB].numberCells--;
        if(bubbles[idB].numberCells<=0) bubbles[idB].reset();
        grid[i][j].setBubbleId(-1);
      }
    }
    
    for(int id=0; id<B ; id++)
      if(bubbles[id].numberCells>0) 
        bubbles[id].deprecated = true;
    
    Stack<Integer> coords = new Stack<>();
        
    for(int i=0; i<Nx ;i++) {
      for(int j=0; j<Ny ;j++) {
        int idN = grid[i][j].getBubbleId(); 
        CELL_TYPE typeN = grid[i][j].getType();
        if((idN<0 && !(typeN==CELL_TYPE.I || typeN==CELL_TYPE.G)) || (idN>=0 && !bubbles[idN].deprecated)) continue; // if(mark || Solid || Liquid)
        
        int currentId = -1;
        for(int id=0; id<B ; id++)
          if(bubbles[id].numberCells==0) 
            { currentId=id; break; }
        // if(currentId==-1) ????
        
        coords.clear();
        coords.push(j*Nx+i);
        
        while (!coords.empty()) {
          int coord = coords.pop();
          int coordX = coord%Nx;
          int coordY = coord/Nx;
          int idC = grid[coordX][coordY].getBubbleId();
          
          if(idC>=0 && !bubbles[idC].deprecated) continue; 
          
          bubbles[currentId].numberCells++;
          bubbles[currentId].volume += 1.f-grid[coordX][coordY].getPhi();
          if(idC >= 0) bubbles[currentId].volumeInit += bubbles[idC].volumeInit/float(bubbles[idC].numberCells); 
          //if(idC >= 0) bubbles[currentId].volumeInit += (1.f-grid[coordX][coordY].getPhi()) * bubbles[idC].rho / cs2;
          grid[coordX][coordY].setBubbleId(currentId);
            
          for (int n=1; n<9 ;n++) {
            int offsetX = mod(coordX+D2Q9_cx[n],Nx);
            int offsetY = mod(coordY+D2Q9_cy[n],Ny);
            int idO = grid[offsetX][offsetY].getBubbleId();
            CELL_TYPE typeO = grid[offsetX][offsetY].getType();
            if((idO<0 && (typeO==CELL_TYPE.I || typeO==CELL_TYPE.G)) || (idO>=0 && bubbles[idO].deprecated)) coords.push(offsetY*Nx+offsetX); // warning, bad if!
          }
        }
        
        if (bubbles[currentId].volumeInit==0.f) bubbles[currentId].volumeInit = bubbles[currentId].volume;
        bubbles[currentId].rho = cs2 * bubbles[currentId].volumeInit/bubbles[currentId].volume;
        println("bubble " + currentId + ", rho="+bubbles[currentId].rho + ", volumeInit=" + bubbles[currentId].volumeInit + ", volume=" + bubbles[currentId].volume);
      }
    }
    
    for(int id=0; id<B ; id++)
      if(bubbles[id].deprecated)
        bubbles[id].reset();
    
    t++;
  }
}
