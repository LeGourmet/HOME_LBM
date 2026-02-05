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
    this.rho = max(0.f, p_rho);
    this.ux = p_ux;
    this.uy = p_uy;
    
    if(p_type == CELL_TYPE.SOLID) {
      this.type = CELL_TYPE.SOLID;
      this.phi = 0.f;
    } else {
      this.phi = (p_phi<0.5f) ? 0.f : 1.f;
      if(this.phi == 0.f){
        this.type = CELL_TYPE.GAS;
        this.rho = 1.f;
        this.ux = 0.f;
        this.uy = 0.f;
      } else this.type = CELL_TYPE.FLUID;
    }
    
    this.Sxx = (ux*ux - cs2);
    this.Syy = (uy*uy - cs2);
    this.Sxy = (ux*uy);  // Syx = Sxy; // symetric tensor
    
    this.mass = phi * rho;
    this.massex = 0.f;
  } 
  
  // ------------------------------------------------------ GETTERS ------------------------------------------------------
  public CELL_TYPE getType()   { return this.type; }
  public float getDensity()    { return this.rho; }
  public float getVelocityX()  { return this.ux; }
  public float getVelocityY()  { return this.uy; }
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
  
  private float computeFiEq(int p_i, float p_rho, float p_ux, float p_uy) {
    return D2Q9_w[p_i] * ( p_rho * (1.f + (D2Q9_cx[p_i]*p_ux + D2Q9_cy[p_i]*p_uy)/cs2 + 0.5f*sq(D2Q9_cx[p_i]*p_ux + D2Q9_cy[p_i]*p_uy)/cs4 - 0.5f*(p_ux*p_ux + p_uy*p_uy)/cs2) - 1.f);
  }
    
  private float computeFi(int p_i, float p_rho, float p_ux, float p_uy, float p_Sxx, float p_Syy, float p_Sxy) {
    return D2Q9_w[p_i] * ( p_rho * ( 1.f +
                                     (D2Q9_cx[p_i]*p_ux + D2Q9_cy[p_i]*p_uy)/cs2 +
                                     0.5f*( 2.f*p_Sxy*D2Q9_cx[p_i]*D2Q9_cy[p_i] + p_Sxx*(D2Q9_cx[p_i]*D2Q9_cx[p_i]-cs2) + p_Syy*(D2Q9_cy[p_i]*D2Q9_cy[p_i]-cs2))/cs4 +
                                     0.5f*( (D2Q9_cx[p_i]*D2Q9_cx[p_i]*D2Q9_cy[p_i]-D2Q9_cy[p_i]*cs2) * (p_Sxx*p_uy+2.f*p_Sxy*p_ux-2.f*p_ux*p_ux*p_uy) +
                                            (D2Q9_cx[p_i]*D2Q9_cy[p_i]*D2Q9_cy[p_i]-D2Q9_cx[p_i]*cs2) * (p_Syy*p_ux+2.f*p_Sxy*p_uy-2.f*p_uy*p_uy*p_ux))/cs6) 
                          - 1.f);
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
    for (int i=1; i<9 ;i++){
      int idNx = mod(p_x+D2Q9_cx[i],p_simulation.getNx());
      int idNy = mod(p_y+D2Q9_cy[i],p_simulation.getNy());
      mass += p_simulation.getCell(idNx,idNy).getMassex();
    }
    
    if (type == CELL_TYPE.FLUID) {
      for (int i=1; i<9 ;i++)
        mass += fin[i] - fon[i]; // neighbor is fluid or interface cell
    }
    else if(type == CELL_TYPE.INTERFACE){
      float[] phiN = new float[9]; // cache fill level of neighbor lattice points
      phiN[0] = rho>0.f ? constrain(mass/rho, 0.f, 1.f) : 0.5f; // don't load phi[n] from memory, instead recalculate it with mass corrected by excess mass
      
      // cache fill level of neighbor lattice points
      for (int i=1; i<9 ;i++){
        int idNx = mod(p_x+D2Q9_cx[i],p_simulation.getNx());
        int idNy = mod(p_y+D2Q9_cy[i],p_simulation.getNy());
        phiN[i] = (p_simulation.getCell(idNx,idNy).getType()==CELL_TYPE.SOLID) ? 1.f : p_simulation.getCell(idNx,idNy).getPhi();
      }
      
      //float rho_laplace = 6.f* st * calculate_curvature(phiN);
      float rho_k = 1.f;
      float def_6_sigma_k = 2.f * st / cs2;
      float rho_laplace = (def_6_sigma_k==0.f) ? 0.f :  def_6_sigma_k * calculate_curvature(phiN);
      
      // limit for stability purpose (simulation can't exceed mach 1)
      float tmpUx = ux + 0.5f * fx;
      float tmpUy = uy + 0.5f * fy;
      float tmpUNorm = sqrt(sq(tmpUx)+sq(tmpUy));
      if(tmpUNorm>cs){
        tmpUx *= cs/tmpUNorm;
        tmpUy *= cs/tmpUNorm;
      }
      
      float[] feq = new float[9]; // reconstruct f from neighbor gas lattice points
      for(int i=0; i<9 ;i++)
        feq[i] = computeFiEq(i, rho_k - rho_laplace, tmpUx, tmpUy); // do some buble stuffs
      
      for (int i=1; i<9 ;i++) {
        int j = (i%2==0) ? i-1 : i+1;
        int idNx = mod(p_x+D2Q9_cx[j],p_simulation.getNx());
        int idNy = mod(p_y+D2Q9_cy[j],p_simulation.getNy());
        
        CELL_TYPE typeN = p_simulation.getCell(idNx,idNy).getType();
        if (typeN == CELL_TYPE.FLUID) mass += fin[i] - fon[j];
        else if (typeN == CELL_TYPE.INTERFACE || typeN == CELL_TYPE.INTERFACE_TO_FLUID || typeN == CELL_TYPE.INTERFACE_TO_GAS || typeN == CELL_TYPE.GAS_TO_INTERFACE) mass += 0.5f * (phiN[j] + phiN[0]) * (fin[i] - fon[j]);
      }
  
      for (int i=1; i<9 ;i++) {
        int j = (i%2==0) ? i-1 : i+1;
        int idNx = mod(p_x+D2Q9_cx[j],p_simulation.getNx());
        int idNy = mod(p_y+D2Q9_cy[j],p_simulation.getNy());
        
        CELL_TYPE typeN = p_simulation.getCell(idNx,idNy).getType();
        if (typeN == CELL_TYPE.GAS) fin[i] = feq[j] - fon[j] + feq[i];
      }
    }

    // ---------------------------------- RECONSTRUCTION ----------------------------------    
     for(int i=0; i<9 ;i++)
      fin[i] += D2Q9_w[i];
    
    float rhoT = fin[0] + fin[1] + fin[2] + fin[3] + fin[4] + fin[5] + fin[6] + fin[7] + fin[8];  // sum(fi);
    float uxT  = (fin[1] - fin[2] + fin[5] - fin[6] + fin[7] - fin[8]) / rhoT + 0.5f * fx;         // sum(fi * cix);
    float uyT  = (fin[3] - fin[4] + fin[5] - fin[6] + fin[8] - fin[7]) / rhoT + 0.5f * fy;         // sum(fi * ciy);
    float SxxT = (fin[1] + fin[2] + fin[5] + fin[6] + fin[7] + fin[8]) / rhoT - cs2;               // sum(fi * (cix*cix-cs2));
    float SyyT = (fin[3] + fin[4] + fin[5] + fin[6] + fin[8] + fin[7]) / rhoT - cs2;               // sum(fi * (ciy*ciy-cs2));
    float SxyT = (fin[5] + fin[6] - fin[7] - fin[8]) / rhoT;                                       // sum(fi * (cix*ciy));
    
    // ---------------------------------- SURFACE TAG UPDATE  ---------------------------------- 
    
    if (type==CELL_TYPE.INTERFACE) {
      boolean TYPE_NO_L = true, TYPE_NO_G = true; // temporary flags for no fluid or gas neighbors
      for (int i=1; i<9 ;i++) {
        int idNx = mod(p_x+D2Q9_cx[i],p_simulation.getNx());
        int idNy = mod(p_y+D2Q9_cy[i],p_simulation.getNy());
        CELL_TYPE typeN = p_simulation.getCell(idNx,idNy).getType();
        TYPE_NO_L = TYPE_NO_L && (typeN != CELL_TYPE.FLUID);
        TYPE_NO_G = TYPE_NO_G && (typeN != CELL_TYPE.GAS);
      }
      if (mass > rhoT || TYPE_NO_G) type = CELL_TYPE.INTERFACE_TO_FLUID; // set flag interface->fluid
      else if (mass < 0.f || TYPE_NO_L) type = CELL_TYPE.INTERFACE_TO_GAS; // set flag interface->gas
    }
    
    // ---------------------------------- COLLISION ----------------------------------
    
    float uTNorm = sqrt(sq(uxT)+sq(uyT));
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
  
  public void surfaceInit(int p_x, int p_y, LBM p_simulation){
     if(type!=CELL_TYPE.GAS) return;
    
    float rhon=0.f, uxn=0.f, uyn=0.f, counter=0.f; // average over all fluid/interface neighbors

    for(int i=1; i<9 ;i++) {
      int idNx = mod(p_x+D2Q9_cx[i],p_simulation.getNx());
      int idNy = mod(p_y+D2Q9_cy[i],p_simulation.getNy());

      CELL_TYPE typeN = p_simulation.getCell(idNx,idNy).getType();
      if (typeN == CELL_TYPE.FLUID) {
        counter += 1.f;
        rhon += p_simulation.getCell(idNx,idNy).getDensity();
        uxn += p_simulation.getCell(idNx,idNy).getVelocityX();
        uyn += p_simulation.getCell(idNx,idNy).getVelocityY();
      }
    }
      
    if(counter > 0.f){
      type = CELL_TYPE.INTERFACE;
      
      rho = rhon / counter;
      ux = uxn / counter;
      uy = uyn / counter;
      Sxx = (ux*ux - cs2);
      Syy = (uy*uy - cs2);
      Sxy = (ux*uy);  // Syx = Sxy; // symetric tensor
      
      phi = 0.5f;
      mass = phi * rho;
    }
  }
  
  public void surface1(int p_x, int p_y, LBM p_simulation) { 
    if(type == CELL_TYPE.INTERFACE_TO_FLUID){
      for(int i=1; i<9 ;i++){
        int idNx = mod(p_x+D2Q9_cx[i],p_simulation.getNx());
        int idNy = mod(p_y+D2Q9_cy[i],p_simulation.getNy());

        CELL_TYPE typeN = p_simulation.getCell(idNx,idNy).getType();
        if(typeN == CELL_TYPE.INTERFACE_TO_GAS) p_simulation.getCell(idNx,idNy).setType(CELL_TYPE.INTERFACE);
        if(typeN == CELL_TYPE.GAS) p_simulation.getCell(idNx,idNy).setType(CELL_TYPE.GAS_TO_INTERFACE);
      }
    }
  }
  
  public void surface2(int p_x, int p_y, LBM p_simulation) { 
    if(type == CELL_TYPE.GAS_TO_INTERFACE){
      float rhon=0.f, uxn=0.f, uyn=0.f, counter=0.f; // average over all fluid/interface neighbors

      for(int i=1; i<9 ;i++) {
        int idNx = mod(p_x+D2Q9_cx[i],p_simulation.getNx());
        int idNy = mod(p_y+D2Q9_cy[i],p_simulation.getNy());

        CELL_TYPE typeN = p_simulation.getCell(idNx,idNy).getType();
        if (typeN == CELL_TYPE.FLUID || typeN == CELL_TYPE.INTERFACE || typeN == CELL_TYPE.INTERFACE_TO_FLUID) { // fluid or interface or (interface->fluid) neighbor
          counter += 1.f;
          rhon += p_simulation.getCell(idNx,idNy).getDensity();
          uxn += p_simulation.getCell(idNx,idNy).getVelocityX();
          uyn += p_simulation.getCell(idNx,idNy).getVelocityY();
        }
      }
      
      rhon = counter > 0.f ? rhon / counter : 1.f;
      uxn = counter > 0.f ? uxn / counter : 0.f;
      uyn = counter > 0.f ? uyn / counter : 0.f;
      
      float[] feq = new float[9];
      for(int i=0; i<9 ;i++)
        feq[i] = computeFiEq(i, rhon, uxn, uyn) + D2Q9_w[i]; // unshift
      
      rho = rhon;
      ux = uxn;
      uy = uyn;
      Sxx = (feq[1] + feq[2] + feq[5] + feq[6] + feq[7] + feq[8]) / rhon - cs2;
      Syy = (feq[3] + feq[4] + feq[5] + feq[6] + feq[8] + feq[7]) / rhon - cs2;
      Sxy = (feq[5] + feq[6] - feq[7] - feq[8]) / rhon;
    } else if(type == CELL_TYPE.INTERFACE_TO_GAS) {
      for(int i=1; i<9 ;i++){
        int idNx = mod(p_x+D2Q9_cx[i],p_simulation.getNx());
        int idNy = mod(p_y+D2Q9_cy[i],p_simulation.getNy());

        CELL_TYPE typeN = p_simulation.getCell(idNx,idNy).getType();
        if(typeN == CELL_TYPE.FLUID || typeN == CELL_TYPE.INTERFACE_TO_FLUID)
          p_simulation.getCell(idNx,idNy).setType(CELL_TYPE.INTERFACE);
      }
    }
  }
  
  public void surface3(int p_x, int p_y, LBM p_simulation) { 
    if(type==CELL_TYPE.SOLID) return;
    
    float rhon = rho;
    float massn = mass;
    float massexn = 0.f;
    
    if (type==CELL_TYPE.FLUID || type==CELL_TYPE.INTERFACE_TO_FLUID) {
      if (type==CELL_TYPE.INTERFACE_TO_FLUID) type = CELL_TYPE.FLUID;
      massexn = massn - rhon; // dump mass-rho difference into excess mass
      massn = rhon; // fluid cell mass has to equal rho
      phi = 1.f;
    }
    else if (type==CELL_TYPE.INTERFACE || type==CELL_TYPE.GAS_TO_INTERFACE) {
      if (type==CELL_TYPE.GAS_TO_INTERFACE)  type = CELL_TYPE.INTERFACE;
      massexn = massn > rhon ? massn - rhon : massn < 0.f ? massn : 0.f; // allow interface cells with mass>rho or mass<0
      massn = constrain(massn, 0.f, rhon);
      phi = rhon > 0.f ? massn / rhon : 0.5f; // calculate fill level for next step (only necessary for interface cells)
    }
    else if (type==CELL_TYPE.GAS || type==CELL_TYPE.INTERFACE_TO_GAS) {
      if (type==CELL_TYPE.INTERFACE_TO_GAS) {
        type = CELL_TYPE.GAS;
        rho = 1.f;
        ux = 0.f;
        uy = 0.f;
      }
      massexn = massn; // dump remaining mass into excess mass
      massn = 0.f;
      phi = 0.f;
    }
    
    int counter = 0; // count (fluid|interface) neighbors
    for (int i=1; i<9 ;i++) {
      int idNx = mod(p_x+D2Q9_cx[i],p_simulation.getNx());
      int idNy = mod(p_y+D2Q9_cy[i],p_simulation.getNy());
      
      CELL_TYPE typeN = p_simulation.getCell(idNx,idNy).getType();
      if(typeN == CELL_TYPE.FLUID || typeN == CELL_TYPE.INTERFACE || typeN == CELL_TYPE.INTERFACE_TO_FLUID || typeN == CELL_TYPE.GAS_TO_INTERFACE) counter++;
    }
    
    mass = massn + ((counter > 0) ? 0.f : massexn); // if excess mass can't be distributed to neighboring interface or fluid cells, add it to local mass (ensure mass conservation)
    massex = ((counter > 0) ? massexn / (float)counter : 0.f); // divide excess mass up for all interface or fluid neighbors
  }
}
