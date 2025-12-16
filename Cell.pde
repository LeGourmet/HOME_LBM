public enum CELL_TYPE { SOLID, FLUID, INTERFACE, GAS, EQUILIBRIUM, INTERFACE_TO_FLUID, INTERFACE_TO_GAS, GAS_TO_INTERFACE };

public class Cell {
  // ----------------------------------------------------- ATTRIBUTS -----------------------------------------------------
  private CELL_TYPE type;
    
  private float rho;
  private float ux;
  private float uy;
  private float Sxx;
  private float Syy;
  private float Sxy;
  
  private float _rho;
  private float _ux;
  private float _uy;
  private float _Sxx;
  private float _Syy;
  private float _Sxy;
  
  private float phi;
  private float mass;
  private float massex;
  
  // --------------------------------------------- DESTRUCTOR / CONSTRUCTOR ----------------------------------------------  
  public Cell(CELL_TYPE p_type, float p_rho, float p_ux, float p_uy, float p_phi) {
    this.type = p_type;
    
    this.rho = max(p_rho,0.f);
    this.ux  = p_ux;
    this.uy  = p_uy;
    this.Sxx = 0.f;
    this.Syy = 0.f;
    this.Sxy = 0.f;
    
    if(this.type == CELL_TYPE.SOLID) {
      this.phi = 0.f;  
    } else if(this.type == CELL_TYPE.EQUILIBRIUM) {
      this.phi = constrain(p_phi, 0.f, 1.f);
      this.Sxx = (this.ux*this.ux - cs2);
      this.Syy = (this.uy*this.uy - cs2);
      this.Sxy = (this.ux*this.uy); // Syx = Sxy; // symetric tensor
    } else {
      this.phi = (p_phi<0.5f) ? 0.f : 1.f;
      if(this.phi==1.f)
        this.type = CELL_TYPE.FLUID;
      else {
        this.type = CELL_TYPE.GAS;
        this.rho = 1.f;
        this.ux = 0.f;
        this.uy = 0.f;
      }
    }
    
    this._rho   = this.rho;
    this._ux  = this.ux;
    this._uy  = this.uy;
    this._Sxx = this.Sxx;
    this._Syy = this.Syy;
    this._Sxy = this.Sxy;  // Syx = Sxy; // symetric tensor
    
    this.mass = this.phi*this.rho;
    this.massex = 0.f;
  } 
  
  // ------------------------------------------------------ GETTERS ------------------------------------------------------
  public CELL_TYPE getType()   { return this.type; }
  public float getDensity()    { return this.rho; }
  public float get_Density()   { return this._rho; }
  public float getVelocityX()  { return this.ux; }
  public float get_VelocityX() { return this._ux; }
  public float getVelocityY()  { return this.uy; }
  public float get_VelocityY() { return this._uy; }
  public float getSxx()        { return this.Sxx; }
  public float getSyy()        { return this.Syy; }
  public float getSxy()        { return this.Sxy; }
  public float getPhi()        { return this.phi; }
  public float getMass()       { return this.mass; }
  public float getMassex()     { return this.massex; }
  
  // ------------------------------------------------------ SETTERS ------------------------------------------------------
  public void setType(CELL_TYPE p_type)      { this.type = p_type; }
  public void setDensity(float p_density)    { this.rho = max(0.f,p_density); }
  public void setVelocityX(float p_velocity) { this.ux = p_velocity; }
  public void setVelocityY(float p_velocity) { this.uy = p_velocity; }
  
  // -------------------------------------------------- FUNCTIONS DDFs ---------------------------------------------------  
  
  // without ddf shifting => feq(i) - w(i) 
  private float computeFiEq(int p_i, float p_rho, float p_ux, float p_uy) {
    return D2Q9_w[p_i] * p_rho * (1.f + (D2Q9_cx[p_i]*p_ux + D2Q9_cy[p_i]*p_uy)/cs2 + 0.5f*sq(D2Q9_cx[p_i]*p_ux + D2Q9_cy[p_i]*p_uy)/cs4 - 0.5f*(p_ux*p_ux + p_uy*p_uy)/cs2);
  }
  
  // without ddf shifting => feq(i) - w(i)
  private float computeFi(int p_i, float p_rho, float p_ux, float p_uy, float p_Sxx, float p_Syy, float p_Sxy) {
    return D2Q9_w[p_i] * p_rho * (1.f +
                                   (D2Q9_cx[p_i]*p_ux + D2Q9_cy[p_i]*p_uy)/cs2 +
                                   0.5f*( 2.f*p_Sxy*D2Q9_cx[p_i]*D2Q9_cy[p_i] + p_Sxx*(D2Q9_cx[p_i]*D2Q9_cx[p_i]-cs2) + p_Syy*(D2Q9_cy[p_i]*D2Q9_cy[p_i]-cs2))/cs4 +
                                   0.5f*( (D2Q9_cx[p_i]*D2Q9_cx[p_i]*D2Q9_cy[p_i]-D2Q9_cy[p_i]*cs2) * (p_Sxx*p_uy+2.f*p_Sxy*p_ux-2.f*p_ux*p_ux*p_uy) +
                                          (D2Q9_cx[p_i]*D2Q9_cy[p_i]*D2Q9_cy[p_i]-D2Q9_cx[p_i]*cs2) * (p_Syy*p_ux+2.f*p_Sxy*p_uy-2.f*p_ux*p_uy*p_uy))/cs6 );
  }
  
  // -------------------------------------------------- FUNCTIONS FLOW ---------------------------------------------------  
  public void flowStreamingCollision(int p_x, int p_y, LBM p_simulation) {    
    if(type==CELL_TYPE.SOLID || type==CELL_TYPE.EQUILIBRIUM || type==CELL_TYPE.GAS) return;
    
    // external forces = body forces : Newton second law of motion : F = M*G = M*A => g = A
    float fx = p_simulation.getForceX(p_x, p_y);
    float fy = p_simulation.getForceY(p_x, p_y);
    
    float nu = nu_fluid;
    float st = st_fluid;
    
    // ---------------------------------- STREAMING ----------------------------------
    float[] fin = new float[9];
    fin[0] = computeFi(0, rho, ux, uy, Sxx, Syy, Sxy);
    float[] fon = new float[9];
    fon[0] = fin[0];
    
    for(int i=1; i<9 ;i++){
      int j = (i%2==0) ? i-1 : i+1;
      
      int idNx = mod(p_x+D2Q9_cx[j],p_simulation.getNx());
      int idNy = mod(p_y+D2Q9_cy[j],p_simulation.getNy());

      CELL_TYPE typeN = p_simulation.getCell(idNx,idNy).getType();
      float rhoN = p_simulation.getCell(idNx,idNy).getDensity();
      float uxN  = p_simulation.getCell(idNx,idNy).getVelocityX();
      float uyN  = p_simulation.getCell(idNx,idNy).getVelocityY();
      float SxxN = p_simulation.getCell(idNx,idNy).getSxx();
      float SyyN = p_simulation.getCell(idNx,idNy).getSyy();
      float SxyN = p_simulation.getCell(idNx,idNy).getSxy();

      if(typeN==CELL_TYPE.SOLID)
        fin[i] = computeFi(i, rho, uxN, uyN, Sxx+uxN*uxN-ux*ux, Syy+uyN*uyN-uy*uy, Sxy+uxN*uyN-ux*uy);
      else if(typeN==CELL_TYPE.EQUILIBRIUM)
        fin[i] = computeFiEq(i, rhoN, uxN, uyN);
      else
        fin[i] = computeFi(i, rhoN, uxN, uyN, SxxN, SyyN, SxyN);
      
      fon[i] = computeFi(i, rho, ux, uy, Sxx, Syy, Sxy);
    }

    // ---------------------------------- SURFACE 0 ---------------------------------- 
    for(int i=1; i<9 ;i++)
      mass += p_simulation.getCell(mod(p_x+D2Q9_cx[i],p_simulation.getNx()), mod(p_y+D2Q9_cy[i],p_simulation.getNy())).getMassex();
    
    if(type == CELL_TYPE.FLUID)
      for(int i=1; i<9 ;i++)
        mass += fin[i] - fon[i];
    else if(type == CELL_TYPE.INTERFACE){
      float phiN[] = new float[9]; // cache fill level of neighbor lattice points
      phiN[0] = rho>0.f ? constrain(mass/rho, 0.f, 1.f) : 0.5f; // don't load phi[n] from memory, instead recalculate it with mass corrected by excess mass
      for(int i=1; i<9 ;i++) {
        int idNx = mod(p_x+D2Q9_cx[i],p_simulation.getNx());
        int idNy = mod(p_y+D2Q9_cy[i],p_simulation.getNy());
        phiN[i] = p_simulation.getCell(idNx, idNy).getPhi(); // cache fill level of neighbor lattice points
        
        // do some other stuff if solid !
      }
      
      // do some buble compute
      
      //float rho_laplace = 6.f* st * calculate_curvature(phiN);
      float rho_laplace = 6.f * st;
      
      // limit for stability purpose (simulation can't exceed mach 1)
      float uxTmp = ux;
      float uyTmp = uy;
      float uTmpNorm = sqrt(sq(uxTmp+0.5f*fx)+sq(uyTmp+0.5f*fy));
      if(uTmpNorm>cs){
        uxTmp *= cs/uTmpNorm;
        uyTmp *= cs/uTmpNorm;
      }
  
      float[] feq = new float[9]; // reconstruct f from neighbor gas lattice points
      for(int i=0; i<9 ;i++)
        feq[i] = computeFiEq(i, 1.f - rho_laplace, uxTmp, uyTmp); // do some buble stuffs
      
      for(int i=1; i<9 ;i++) {
        int j = (i%2==0) ? i-1 : i+1;        
        CELL_TYPE typeN = p_simulation.getCell(mod(p_x+D2Q9_cx[i],p_simulation.getNx()), mod(p_y+D2Q9_cy[i],p_simulation.getNy())).getType();
             if (typeN == CELL_TYPE.FLUID) mass += fin[j] - fon[i];
        else if (typeN == CELL_TYPE.INTERFACE || typeN == CELL_TYPE.INTERFACE_TO_FLUID || typeN == CELL_TYPE.INTERFACE_TO_GAS || typeN == CELL_TYPE.GAS_TO_INTERFACE) mass += 0.5f * (phiN[i] + phiN[0]) * (fin[j] - fon[i]);
      }
  
      for(int i=1; i<9 ;i++) {
        int j = (i%2==0) ? i-1 : i+1;
        if(p_simulation.getCell(mod(p_x+D2Q9_cx[i],p_simulation.getNx()), mod(p_y+D2Q9_cy[i],p_simulation.getNy())).getType() == CELL_TYPE.GAS) fin[i] = feq[j] - fon[j] + feq[i];
      }
    }

    // ---------------------------------- RECONSTRUCTION ----------------------------------    
    float rhoT = fin[0] + fin[1] + fin[2] + fin[3] + fin[4] + fin[5] + fin[6] + fin[7] + fin[8];
    float uxT  = (fin[1] - fin[2] + fin[5] - fin[6] + fin[7] - fin[8]) / rhoT;
    float uyT  = (fin[3] - fin[4] + fin[5] - fin[6] + fin[8] - fin[7]) / rhoT;    
    float SxxT = (fin[1] + fin[2] + fin[5] + fin[6] + fin[7] + fin[8]) / rhoT - cs2;
    float SyyT = (fin[3] + fin[4] + fin[5] + fin[6] + fin[8] + fin[7]) / rhoT - cs2;
    float SxyT = (fin[5] + fin[6] - fin[7] - fin[8]) / rhoT;
    
    // ---------------------------------- SURFACE 0.5 ---------------------------------- 
    
    if (type == CELL_TYPE.INTERFACE) {
      boolean TYPE_NO_F = true, TYPE_NO_G = true; // temporary flags for no fluid or gas neighbors
      for (int i=1; i<9 ;i++) {
        CELL_TYPE typeN = p_simulation.getCell(mod(p_x+D2Q9_cx[i],p_simulation.getNx()), mod(p_y+D2Q9_cy[i],p_simulation.getNy())).getType();
        TYPE_NO_F = TYPE_NO_F && typeN != CELL_TYPE.FLUID;
        TYPE_NO_G = TYPE_NO_G && typeN != CELL_TYPE.GAS;
      }
           if (mass > rhoT || TYPE_NO_G) type = CELL_TYPE.INTERFACE_TO_FLUID; // set flag interface->fluid
      else if (mass < 0.f || TYPE_NO_F) type = CELL_TYPE.INTERFACE_TO_GAS; // set flag interface->gas
    }
    // do some stuff with bubble and Stress tensor
    
    // ---------------------------------- COLLISION ----------------------------------
    
    float uTNorm = sqrt(sq(uxT+0.5f*fx)+sq(uyT+0.5f*fy));
    if(uTNorm>cs){
      uxT *= cs/uTNorm;
      uyT *= cs/uTNorm;
    }
    
    float tau = nu/cs2 + 0.5f;
    
    // compute moment at time t+1
    _rho = rhoT; 
    _ux = uxT + 0.5f*fx; // (Guo forcing, Krueger p.233f)
    _uy = uyT + 0.5f*fy; // (Guo forcing, Krueger p.233f)
    _Sxx = (tau-1.f)/(2.f*tau)*(SxxT-SyyT+uyT*uyT+fx*uxT-fy*uyT) + (tau+1.f)/(2.f*tau)*uxT*uxT + fx*uxT;
    _Syy = (tau-1.f)/(2.f*tau)*(SyyT-SxxT+uxT*uxT+fy*uyT-fx*uxT) + (tau+1.f)/(2.f*tau)*uyT*uyT + fy*uyT;
    _Sxy = (1.f-1.f/tau)*SxyT + uxT*uyT/tau + (2.f*tau-1.f)/(2.f*tau)*(fx*uyT+fy*uxT);
  }
  
  public void swapMoments() {     
    this.rho = this._rho;
    this.ux = this._ux;
    this.uy = this._uy;
    this.Sxx = this._Sxx;
    this.Syy = this._Syy;
    this.Sxy = this._Sxy;
  }
  
  // -------------------------------------------------- FUNCTIONS VOF ---------------------------------------------------  
  
  public void surfaceVOF_init(int p_x, int p_y, LBM p_simulation){
    if(type != CELL_TYPE.GAS) return;
    
    float rhot = 0.f, uxt = 0.f, uyt = 0.f, counter = 0.f; // average over all fluid/interface neighbors
    for (int i=1; i<9 ;i++){
      int idNx = mod(p_x+D2Q9_cx[i],p_simulation.getNx());
      int idNy = mod(p_y+D2Q9_cy[i],p_simulation.getNy());
  
      if (p_simulation.getCell(idNx,idNy).getType() == CELL_TYPE.FLUID) {
        counter += 1.f;
        rhot += p_simulation.getCell(idNx,idNy).getDensity();
        uxt  += p_simulation.getCell(idNx,idNy).getVelocityX();
        uyt  += p_simulation.getCell(idNx,idNy).getVelocityY();
      }
    }
    
    if (counter == 0.f) return;
    
    type = CELL_TYPE.INTERFACE;
    
    rho = rhot / counter;
    ux  = uxt / counter;
    uy  = uyt / counter;
    
   
    Sxx = (ux*ux - cs2);
    Syy = (uy*uy - cs2);
    Sxy = (ux*uy);

    _rho = rho;
    _ux  = ux;
    _uy  = uy;
    _Sxx = Sxx;
    _Syy = Syy;
    _Sxy = Sxy;  // Syx = Sxy; // symetric tensor
    
    phi = 0.5f; 
    mass = phi*rho;
    massex = 0.f; 
  }
  
  public void sufaceVOF_1(int p_x, int p_y, LBM p_simulation) {
    if (type == CELL_TYPE.INTERFACE_TO_FLUID) {
      for(int i=1; i<9; i++) {
        int idNx = mod(p_x+D2Q9_cx[i],p_simulation.getNx());
        int idNy = mod(p_y+D2Q9_cy[i],p_simulation.getNy());
        
        CELL_TYPE typeN = p_simulation.getCell(idNx,idNy).getType();
             if (typeN == CELL_TYPE.INTERFACE_TO_GAS) p_simulation.getCell(idNx,idNy).setType(CELL_TYPE.INTERFACE); // prevent interface neighbor cells from becoming gas
        else if (typeN == CELL_TYPE.GAS             ) p_simulation.getCell(idNx,idNy).setType(CELL_TYPE.GAS_TO_INTERFACE); // neighbor cell was gas and must change to interface
      }
    }
  }
  
  public void sufaceVOF_2(int p_x, int p_y, LBM p_simulation) {
    if (type == CELL_TYPE.GAS_TO_INTERFACE) { // initialize the fi of gas cells that should become interface
      float rhot = 0.f, uxt = 0.f, uyt = 0.f, counter = 0.f; // average over all fluid/interface neighbors
      for (int i=1; i<9 ;i++){
        int idNx = mod(p_x+D2Q9_cx[i],p_simulation.getNx());
        int idNy = mod(p_y+D2Q9_cy[i],p_simulation.getNy());
        
        CELL_TYPE typeN = p_simulation.getCell(idNx,idNy).getType(); 
        if (typeN == CELL_TYPE.FLUID || typeN == CELL_TYPE.INTERFACE || typeN == CELL_TYPE.INTERFACE_TO_FLUID) {
          counter += 1.f;
          rhot += p_simulation.getCell(idNx,idNy).get_Density();
          uxt  += p_simulation.getCell(idNx,idNy).get_VelocityX();
          uyt  += p_simulation.getCell(idNx,idNy).get_VelocityY();
        }
      }
      
      _rho = counter > 0.f ? rhot  / counter : 1.f;
      _ux  = counter > 0.f ? uxt / counter : 0.f;
      _uy  = counter > 0.f ? uyt / counter : 0.f;
      // different in FREE-HOME-LBM => use eq without ddf shifting ...
      _Sxx = (_ux*_ux - cs2);
      _Syy = (_uy*_uy - cs2);
      _Sxy = (_ux*_uy);
    } else if (type == CELL_TYPE.INTERFACE_TO_GAS) { // flag interface->gas is set
      for(int i=1; i<9; i++) {
        int idNx = mod(p_x+D2Q9_cx[i],p_simulation.getNx());
        int idNy = mod(p_y+D2Q9_cy[i],p_simulation.getNy());
        
        CELL_TYPE typeN = p_simulation.getCell(idNx,idNy).getType();
        if (typeN == CELL_TYPE.FLUID || typeN == CELL_TYPE.INTERFACE_TO_FLUID)
          type = CELL_TYPE.INTERFACE; // in fluidx and base code => prevent neighboors become fluid ?
          //p_simulation.getCell(idNx,idNy).setType(CELL_TYPE.INTERFACE); // prevent fluid or interface neighbors that turn to fluid from being/becoming fluid
          // + buble merge detector
      }
    }
  }

  public void sufaceVOF_3(int p_x, int p_y, LBM p_simulation) {
    //+  mlflow[0].previous_tag[curind] = mlflow[0].tag_matrix[curind];
    //+   mlflow[0].tag_matrix[curind] = -1;
    
    if (type == CELL_TYPE.FLUID || type == CELL_TYPE.INTERFACE_TO_FLUID) {
      type = CELL_TYPE.FLUID;
      massex = mass-_rho; // dump mass-rho difference into excess mass
      mass = _rho; // fluid cell mass has to equal rho
      phi = 1.f;
      // mlflow[0].previous_tag[curind] = mlflow[0].tag_matrix[curind];
      // mlflow[0].tag_matrix[curind] = -1;
      // if IL => if (mlflow[0].previous_tag[curind] > 0) atomicExch(&mlflow[0].split_flag, 1);
    } else if (type == CELL_TYPE.INTERFACE || type == CELL_TYPE.GAS_TO_INTERFACE) {
      type = CELL_TYPE.INTERFACE;
      massex = mass > _rho ? mass-_rho : mass<0.f ? mass : 0.f; // allow interface cells with mass>rho or mass<0
      mass = constrain(mass, 0.f, _rho);
      phi = _rho>0.f ? mass/_rho : 0.5f ; // calculate fill level for next step (only necessary for interface cells)
    } else if (type == CELL_TYPE.GAS || type == CELL_TYPE.INTERFACE_TO_GAS) {
      type = CELL_TYPE.GAS;
      massex = mass; // dump remaining mass into excess mass
      mass = 0.f;
      phi = 0.f;
      //_rho = 1.f;
      //_ux = 0.f;
      //_uy = 0.f;
      //_Sxx = -cs2;
      //_Syy = -cs2;
      //_Sxy = 0.f;
    } else return;
  
    float counter = 0.f; // count (fluid|interface) neighbors
    for (int i=1; i<9 ;i++) {
      int idNx = mod(p_x+D2Q9_cx[i],p_simulation.getNx());
      int idNy = mod(p_y+D2Q9_cy[i],p_simulation.getNy());
        
      CELL_TYPE typeN = p_simulation.getCell(idNx,idNy).getType();
      if(typeN==CELL_TYPE.FLUID || typeN==CELL_TYPE.INTERFACE || typeN==CELL_TYPE.INTERFACE_TO_FLUID || typeN==CELL_TYPE.GAS_TO_INTERFACE) 
        counter += 1.f;
    }
    
    mass += (counter>0.f) ? 0.f : massex; // if excess mass can't be distributed to neighboring interface or fluid cells, add it to local mass (ensure mass conservation)
    massex = (counter>0.f) ? massex/counter : 0.f; // divide excess mass up for all interface or fluid neighbors
  }
}
