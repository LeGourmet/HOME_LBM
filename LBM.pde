public enum COLOR_TYPE { PRESSURE, VELOCITY, TYPE, BUBBLE };

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
  private Bubble bubbles[];
  
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
    
    this.bubbles = new Bubble[this.Nx*this.Ny];
    for(int i=0; i<this.Nx*this.Ny ;i++)
      bubbles[i] = new Bubble();
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
  
  public Bubble getBubble(int p_id) {
    if(p_id>=(Nx*Ny) || p_id<0) return null; 
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
    
    for(int i=0; i<Nx ;i++)
      for(int j=0; j<Ny ;j++) {
        if(grid[i][j].getBubbleId() != -1 || (grid[i][j].getType() != CELL_TYPE.I && grid[i][j].getType() != CELL_TYPE.G)) continue;
        coords.clear();
        coords.push(j*Nx+i);
        
        while (!coords.empty()) {
            int coord = coords.pop();
            int coordX = coord%Nx;
            int coordY = coord/Nx;
            
            grid[coordX][coordY].bubbleId = currentLabel;
            bubbles[currentLabel].numberCells++;
            bubbles[currentLabel].volumeInit += 1.f-grid[coordX][coordY].getPhi();
            
            for (int k=1; k<9 ;k++) {
              int offsetX = mod(coordX+D2Q9_cx[k],Nx);
              int offsetY = mod(coordY+D2Q9_cy[k],Ny);
              if(grid[offsetX][offsetY].getBubbleId() == -1 && (grid[offsetX][offsetY].getType() == CELL_TYPE.I || grid[offsetX][offsetY].getType() == CELL_TYPE.G)) coords.push(offsetY*Nx+offsetX);
            }
        }
        
        currentLabel++;
      }
    
    // ----------------------- Bubble Init -----------------------
    for(int i=0; i<Nx*Ny ;i++) {
      if(bubbles[i].numberCells==0) continue;
      bubbles[i].rho = cs2;
      bubbles[i].volume = bubbles[i].volumeInit;
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
    for(int i=0; i<Nx ;i++)
      for(int j=0; j<Ny ;j++) {
        grid[i][j].setBubbleIdPrev(grid[i][j].getBubbleId());
        grid[i][j].setBubbleId(-1);
      }
      
    Stack<Integer> coords = new Stack<>();
    IntList labels = new IntList();
    
    // --- split bubbles ---
    for(int id=0; id<Nx*Ny ;id++){
      if(bubbles[id].volumeInit == 0.f) continue;
      
      int start = labels.size();
      
      for(int i=0; i<Nx ;i++)
        for(int j=0; j<Ny ;j++) {
          if(grid[i][j].getBubbleId() >= 0 || grid[i][j].getBubbleIdPrev() != id || !(grid[i][j].getType()==CELL_TYPE.I || grid[i][j].getType()==CELL_TYPE.G)) continue;
            
          coords.clear();
          coords.push(j*Nx+i);
          labels.append(0);
            
          while (!coords.empty()) {
            int coord = coords.pop();
            int coordX = coord%Nx;
            int coordY = coord/Nx;
              
            int currentLabel = labels.size()-1;
            labels.increment(currentLabel);
            grid[coordX][coordY].setBubbleId(currentLabel);
                
            for (int n=1; n<9 ;n++) {
              int offsetX = mod(coordX+D2Q9_cx[n],Nx);
              int offsetY = mod(coordY+D2Q9_cy[n],Ny);
              if(grid[offsetX][offsetY].getBubbleId() < 0 && grid[offsetX][offsetY].getBubbleIdPrev() == id && (grid[offsetX][offsetY].getType() == CELL_TYPE.I || grid[offsetX][offsetY].getType() == CELL_TYPE.G)) coords.push(offsetY*Nx+offsetX);
            }
          }
        }
      
      if((labels.size()-start)==0) { bubbles[id].volumeInit=0.f; }
      else if((labels.size()-start)==1) {  }
      else {
        int ttNB = 0;
        for(int l=start ; l<labels.size() ;l++) ttNB += labels.get(l);
        for(int l=start ; l<labels.size() ;l++) {
          // split here
          // new id (write prev)
          // vo=V0*NB/ttNB
          // delete old bubble id (v0=0)
        }
      }
    }
    
    for(int i=0; i<Nx ;i++)
      for(int j=0; j<Ny ;j++)
        grid[i][j].setBubbleId(-1);
    
    // --- merge bubbles ---    
    for(int i=0; i<Nx ;i++)
      for(int j=0; j<Ny ;j++) {
        if(!grid[i][j].getBubbleMerge()) continue;  
        labels.clear();
        coords.clear();
        coords.push(j*Nx+i);
        
        while (!coords.empty()) {
          int coord = coords.pop();
          int coordX = coord%Nx;
          int coordY = coord/Nx;
            
          grid[coordX][coordY].setBubbleId(-2);
          grid[coordX][coordY].setBubbleMergeFalse();
          
          int id = grid[coordX][coordY].getBubbleIdPrev();
          if(id >= 0 && !labels.hasValue(id)) labels.append(id);
            
          for (int n=1; n<9 ;n++) {
            int offsetX = mod(coordX+D2Q9_cx[n],Nx);
            int offsetY = mod(coordY+D2Q9_cy[n],Ny);
            if(grid[offsetX][offsetY].getBubbleId() != -2 && (grid[offsetX][offsetY].getType() == CELL_TYPE.I || grid[offsetX][offsetY].getType() == CELL_TYPE.G)) coords.push(offsetY*Nx+offsetX);
          }
        }
        
        int currentBubbleId = -1;
        
        // create bubble
        if(labels.size() == 0) {
          for(int id=0; id<Nx*Ny ;id++)
            if(bubbles[id].volumeInit==0.f) 
              { currentBubbleId = id; break; }

          for(int x=0; x<Nx ;x++)
            for(int y=0; y<Ny ;y++)
              if(grid[x][y].getBubbleId()==-2) 
                bubbles[currentBubbleId].volumeInit += 1.f-grid[x][y].getPhi();
        } 
        
        // enlarge a bubble
        else if(labels.size() == 1) {
          currentBubbleId = labels.get(0);
        } 
        
        // merge bubbles
        else {
          currentBubbleId = labels.min();
          for(int l : labels)
            if(l != currentBubbleId)
              bubbles[currentBubbleId].volumeInit += bubbles[l].volumeInit;
        }
        
        // update id
        for(int x=0; x<Nx ;x++)
          for(int y=0; y<Ny ;y++)
            if(grid[x][y].getBubbleId()==-2) 
              grid[x][y].setBubbleId(currentBubbleId);
      }
    
    
    // ------------------- Bubble Update Volume ------------------
    for(int i=0; i<Nx*Ny ;i++) bubbles[i].volume = 0.f;
    
    for(int i=0; i<Nx ;i++)
      for(int j=0; j<Ny ;j++) {
        int id = grid[i][j].getBubbleId();
        if(id>=0) bubbles[id].volume += 1.f-grid[i][j].getPhi();
      }
    
    // --------------------- Bubble Update Rho -------------------
    for(int i=0; i<Nx*Ny ;i++) 
      if(bubbles[i].volume==0.f) {
        bubbles[i].rho = 0.f;
        bubbles[i].volumeInit = 0.f;
      } else 
        bubbles[i].rho = cs2 * bubbles[i].volumeInit/bubbles[i].volume;
      
    t++;
  }
}
