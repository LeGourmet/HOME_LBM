public enum CELL_TYPE { L, I, G, IL, IG, GI, EQUILIBRIUM, SOLID };

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
    if(p_type == CELL_TYPE.SOLID) {
      this.type = CELL_TYPE.SOLID;
      this.phi = 0.f;
      this.rho = p_rho;
      this.ux = p_ux;
      this.uy = p_uy;
    } else {
      this.phi = (p_phi<0.5f) ? 0.f : 1.f;
      if(this.phi == 0.f){
        this.type = CELL_TYPE.G;
        this.rho = 1.f;
        this.ux = 0.f;
        this.uy = 0.f;
      } else {
        this.type = CELL_TYPE.L;
        this.rho = p_rho;
        this.ux = p_ux;
        this.uy = p_uy;
      }
    }
    
    this.Sxx = (ux*ux - cs2);
    this.Syy = (uy*uy - cs2);
    this.Sxy = (ux*uy);  // Syx = Sxy; // symetric tensor
    
    this.mass = phi * rho;
    this.massex = 0.f;
  } 
  
  // ------------------------------------------------------ GETTERS ------------------------------------------------------
  public CELL_TYPE getType() { return this.type; }
  public float getDensity() { return this.rho; }
  public float getVelocityX() { return this.ux; }
  public float getVelocityY() { return this.uy; }
  public float getSxx() { return this.Sxx; }
  public float getSyy() { return this.Syy; }
  public float getSxy() { return this.Sxy; }
  public float getPhi() { return this.phi; }
  public float getMass() { return this.mass; }
  public float getMassex() { return this.massex; }
  
  // ------------------------------------------------------ SETTERS ------------------------------------------------------
  public void setType(CELL_TYPE p_type) { this.type = p_type; }
  public void setDensity(float p_rho) { this.rho = max(0.f,p_rho); }
  public void setVelocityX(float p_velocity) { this.ux = p_velocity; }
  public void setVelocityY(float p_velocity) { this.uy = p_velocity; }
  
  // -------------------------------------------------- INIT ---------------------------------------------------  
  public void init(int p_x, int p_y, LBM p_simulation) {
    if(type!=CELL_TYPE.G) return;
    
    float rhon=0.f, uxn=0.f, uyn=0.f, counter=0.f; // average over all fluid/interface neighbors

    for(int i=1; i<9 ;i++) {
      int idNx = mod(p_x+D2Q9_cx[i],p_simulation.getNx());
      int idNy = mod(p_y+D2Q9_cy[i],p_simulation.getNy());

      CELL_TYPE typeN = p_simulation.getCell(idNx,idNy).getType();
      if (typeN == CELL_TYPE.L) {
        counter += 1.f;
        rhon += p_simulation.getCell(idNx,idNy).getDensity();
        uxn += p_simulation.getCell(idNx,idNy).getVelocityX();
        uyn += p_simulation.getCell(idNx,idNy).getVelocityY();
      }
    }
      
    if(counter > 0.f){
      type = CELL_TYPE.I;
      
      rho = rhon / counter;
      ux = uxn / counter;
      uy = uyn / counter;
      Sxx = (ux*ux - cs2);
      Syy = (uy*uy - cs2);
      Sxy = (ux*uy);  // Syx = Sxy; // symetric tensor
      
      mass = phi * rho;
    }
  }
  
  // -------------------------------------------------- FUNCTIONS DDFs ---------------------------------------------------  
  private float computeFiEq(int p_i, float p_rho, float p_ux, float p_uy) {
    return D2Q9_w[p_i] * ( p_rho * (1.f + (D2Q9_cx[p_i]*p_ux + D2Q9_cy[p_i]*p_uy)/cs2 + 0.5f*sq(D2Q9_cx[p_i]*p_ux + D2Q9_cy[p_i]*p_uy)/cs4 - 0.5f*(p_ux*p_ux + p_uy*p_uy)/cs2) - 1.f);
  }
  
  private float computeFi(int p_i, float p_rho, float p_ux, float p_uy, float p_Sxx, float p_Syy, float p_Sxy) {
    return D2Q9_w[p_i] * ( p_rho * ( 1.f +
                                     (D2Q9_cx[p_i]*p_ux + D2Q9_cy[p_i]*p_uy)/cs2 +
                                     0.5f*( 2.f*p_Sxy*D2Q9_cx[p_i]*D2Q9_cy[p_i] + p_Sxx*(D2Q9_cx[p_i]*D2Q9_cx[p_i]-cs2) + p_Syy*(D2Q9_cy[p_i]*D2Q9_cy[p_i]-cs2))/cs4 +
                                     0.5f*( (D2Q9_cx[p_i]*D2Q9_cx[p_i]*D2Q9_cy[p_i]-D2Q9_cy[p_i]*cs2) * (p_Sxx*p_uy+2.f*p_Sxy*p_ux-2.f*p_ux*p_ux*p_uy) +
                                            (D2Q9_cx[p_i]*D2Q9_cy[p_i]*D2Q9_cy[p_i]-D2Q9_cx[p_i]*cs2) * (p_Syy*p_ux+2.f*p_Sxy*p_uy-2.f*p_ux*p_uy*p_uy))/cs6) 
                          - 1.f);
  }
  
  // -------------------------------------------------- FUNCTIONS FLOW ---------------------------------------------------  
  public void flowStreamingCollision(int p_x, int p_y, LBM p_simulation) {
    if(type==CELL_TYPE.SOLID || type==CELL_TYPE.EQUILIBRIUM || type==CELL_TYPE.G ) return;
     
    float tau = 0.5f + nu/cs2;
     
    // external forces = body forces : Newton second law of motion : F = M*G = M*A => g = A
    float fx = p_simulation.getForceX(p_x, p_y);
    float fy = p_simulation.getForceY(p_x, p_y);
    
    // ---------------------------------- STREAMING ----------------------------------
    float[] fin = new float[9];
    float[] fout = new float[9];
    
    float fi = computeFi(0, rho, ux, uy, Sxx, Syy, Sxy);
    fin[0] = fi;
    fout[0] = fi;
    
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
        
      fout[i] = computeFi(i, rho, ux, uy, Sxx, Syy, Sxy);
    }
    
    // ---------------------------------- SURFACE 0 ----------------------------------    
    float massn = mass;
    
    for (int i=1; i<9 ;i++){
      int idNx = mod(p_x+D2Q9_cx[i],p_simulation.getNx());
      int idNy = mod(p_y+D2Q9_cy[i],p_simulation.getNy());
      massn += p_simulation.getCell(idNx,idNy).getMassex();
    }
    
    if (type == CELL_TYPE.L) {
      for (int i=1; i<9 ;i++)
        massn += fin[i] - fout[i]; // neighbor is fluid or interface cell
    }
    else if (type == CELL_TYPE.I) {
      float[] phiM = new float[9]; // cache fill level of neighbor lattice points
      phiM[0] = rho>0.f ? constrain(massn/rho, 0.f, 1.f) : 0.5f; // don't load phi[n] from memory, instead recalculate it with mass corrected by excess mass
      
      // cache fill level of neighbor lattice points
      for (int i=1; i<9 ;i++){
        int idNx = mod(p_x+D2Q9_cx[i],p_simulation.getNx());
        int idNy = mod(p_y+D2Q9_cy[i],p_simulation.getNy());
        phiM[i] = p_simulation.getCell(idNx,idNy).getPhi();
      }
      
      float def_6_sigma_k = 2.f * ca / cs2;
      float rho_laplace = def_6_sigma_k * 1.f; //*calculate_curvature(phiM);
      float rho_k = 1.f;
      
      /*long bubbleId = bubbleTag[n];
      if (bubbleId >= 0ll) {
        rho_k = bubbleInfo[(ulong)(bubbleId)];                                    // for bubble pressure
        if (bubbleInfo[2ul*def_bubbles_N + (ulong)(bubbleId)] > 5000000.f) def_6_sigma_k = 1e-6f;          // for air layer surface tension
        if ((def_6_sigma_k>1e-3f) && (bubbleInfo[def_bubbles_N + (ulong)(bubbleId)] < 64.f)) def_6_sigma_k = 2e-4f;  // for small bubble surface tension
      }*/
    
      // limit for stability purpose (simulation can't exceed mach 1)
      float tmpUx = ux + 0.5f * fx;
      float tmpUy = uy + 0.5f * fy;
      float tmpUNorm = sqrt(sq(tmpUx)+sq(tmpUy));
      if(tmpUNorm>cs){
        tmpUx *= cs/tmpUNorm;
        tmpUy *= cs/tmpUNorm;
      }
    
      // calculate gas equilibrium DDFs with constant ambient pressure => rho = (p(g) - 2*gamma*k(x))/def_cs2
      float[] feq = new float[9]; // reconstruct f from neighbor gas lattice points
      for (int i=0; i<9 ;i++)
        feq[i] = computeFiEq(i, rho_k - rho_laplace, tmpUx, tmpUy);
    
      for (int i=1; i<9 ;i++) {
        int j = (i%2==0) ? i-1 : i+1;
        int idNx = mod(p_x+D2Q9_cx[j],p_simulation.getNx());
        int idNy = mod(p_y+D2Q9_cy[j],p_simulation.getNy());
        
        CELL_TYPE typeN = p_simulation.getCell(idNx,idNy).getType();
        if (typeN == CELL_TYPE.L) massn += fin[i] - fout[j];
        else if (typeN == CELL_TYPE.I || typeN == CELL_TYPE.IL || typeN == CELL_TYPE.IG || typeN == CELL_TYPE.GI) massn += 0.5f * (phiM[j] + phiM[0]) * (fin[i] - fout[j]);
      }
    
      for (int i=1; i<9 ;i++) {
        int j = (i%2==0) ? i-1 : i+1;
        int idNx = mod(p_x+D2Q9_cx[j],p_simulation.getNx());
        int idNy = mod(p_y+D2Q9_cy[j],p_simulation.getNy());
        
        CELL_TYPE typeN = p_simulation.getCell(idNx,idNy).getType();
        if (typeN == CELL_TYPE.G) fin[i] = feq[j] - fout[j] + feq[i];
      }
    }
    mass = massn;
    
    // ---------------------------------- RECONSTRUCTION ----------------------------------
    for(int i=0; i<9 ;i++)
      fin[i] += D2Q9_w[i];
    
    float rhoT =  fin[0] + fin[1] + fin[2] + fin[3] + fin[4] + fin[5] + fin[6] + fin[7] + fin[8];  // sum(fi);
    float uxT  = (fin[1] - fin[2] + fin[5] - fin[6] + fin[7] - fin[8]) / rhoT + 0.5f * fx;         // sum(fi * cix);
    float uyT  = (fin[3] - fin[4] + fin[5] - fin[6] + fin[8] - fin[7]) / rhoT + 0.5f * fy;         // sum(fi * ciy);
    float SxxT = (fin[1] + fin[2] + fin[5] + fin[6] + fin[7] + fin[8] + cs2) / rhoT - cs2;         // sum(fi * (cix*cix-cs2));
    float SyyT = (fin[3] + fin[4] + fin[5] + fin[6] + fin[8] + fin[7] + cs2) / rhoT - cs2;         // sum(fi * (ciy*ciy-cs2));
    float SxyT = (fin[5] + fin[6] - fin[7] - fin[8]) / rhoT;                                       // sum(fi * (cix*ciy));
    
    // ---------------------------------- SURFACE TAG UPDATE ----------------------------------
    
    if (type==CELL_TYPE.I) {
      boolean TYPE_NO_L = true, TYPE_NO_G = true; // temporary flags for no fluid or gas neighbors
      for (int i=1; i<9 ;i++) {
        int idNx = mod(p_x+D2Q9_cx[i],p_simulation.getNx());
        int idNy = mod(p_y+D2Q9_cy[i],p_simulation.getNy());
        CELL_TYPE typeN = p_simulation.getCell(idNx,idNy).getType();
        TYPE_NO_L = TYPE_NO_L && (typeN != CELL_TYPE.L);
        TYPE_NO_G = TYPE_NO_G && (typeN != CELL_TYPE.G);
      }
      if (mass > rhoT || TYPE_NO_G) type = CELL_TYPE.IL; // set flag interface->fluid
      else if (mass < 0.f || TYPE_NO_L) type = CELL_TYPE.IG; // set flag interface->gas
    }
    
    /*uint x, y, z;
    coordinates(n, &x, &y, &z);
    
    // LES model for neighborhood of bubbles (add eddy viscosity)
    for (int a=-3; a<3 ;a++)
      for (int b=-3; b<3 ;b++)
        for (int c=-3; c<3 ;c++) {
          ulong m = index(((int)(x)+a+(int)(def_Nx))%def_Nx, ((int)(y)+b+(int)(def_Ny))%def_Ny, ((int)(z)+c+(int)(def_Nz))%def_Nz);
          long bubbleId = bubbleTag[m];
          if ((bubbleId>=0) && (bubbleInfo[def_bubbles_N + (ulong)(bubbleId)] < 5000000.f)) {
            tau = (_nu + 4.f * sqrt(sq(_Sxx) + sq(_Syy) + sq(_Szz) + 2.f * (sq(_Sxy) + sq(_Sxz) + sq(_Syz)))) / def_cs2 + 0.5f;
            break;
          }
        }*/
    
    // ---------------------------------- COLLISION ----------------------------------
    float uTNorm = sqrt(sq(uxT)+sq(uyT));
    if(uTNorm>cs){
      uxT *= cs/uTNorm;
      uyT *= cs/uTNorm;
    }
    
    // compute moment at time t+1
    _rho = rhoT; 
    _ux = uxT + 0.5f*fx; // (Guo forcing, Krueger p.233f)
    _uy = uyT + 0.5f*fy; // (Guo forcing, Krueger p.233f)
    _Sxx = (tau-1.f)/(2.f*tau)*(SxxT-SyyT+uyT*uyT+fx*uxT-fy*uyT) + (tau+1.f)/(2.f*tau)*uxT*uxT + fx*uxT;
    _Syy = (tau-1.f)/(2.f*tau)*(SyyT-SxxT+uxT*uxT+fy*uyT-fx*uxT) + (tau+1.f)/(2.f*tau)*uyT*uyT + fy*uyT;
    _Sxy = (1.f-1.f/tau)*SxyT + uxT*uyT/tau + (2.f*tau-1.f)/(2.f*tau)*(fx*uyT+fy*uxT);
  }  
  
  public void flowSwapMoments(int p_x, int p_y, LBM p_simulation) { 
    this.rho = this._rho;
    this.ux = this._ux;
    this.uy = this._uy;
    this.Sxx = this._Sxx;
    this.Syy = this._Syy;
    this.Sxy = this._Sxy;
  }
  
  public void surface1(int p_x, int p_y, LBM p_simulation) { 
    if(type == CELL_TYPE.IL){
      for(int i=1; i<9 ;i++){
        int idNx = mod(p_x+D2Q9_cx[i],p_simulation.getNx());
        int idNy = mod(p_y+D2Q9_cy[i],p_simulation.getNy());

        CELL_TYPE typeN = p_simulation.getCell(idNx,idNy).getType();
        if(typeN == CELL_TYPE.IG) p_simulation.getCell(idNx,idNy).setType(CELL_TYPE.I);
        if(typeN == CELL_TYPE.G) p_simulation.getCell(idNx,idNy).setType(CELL_TYPE.GI);
      }
    }
  }
  
  public void surface2(int p_x, int p_y, LBM p_simulation) { 
    if(type == CELL_TYPE.GI){
      float rhon=0.f, uxn=0.f, uyn=0.f, counter=0.f; // average over all fluid/interface neighbors

      for(int i=1; i<9 ;i++) {
        int idNx = mod(p_x+D2Q9_cx[i],p_simulation.getNx());
        int idNy = mod(p_y+D2Q9_cy[i],p_simulation.getNy());

        CELL_TYPE typeN = p_simulation.getCell(idNx,idNy).getType();
        if (typeN == CELL_TYPE.L || typeN == CELL_TYPE.I || typeN == CELL_TYPE.IL) { // fluid or interface or (interface->fluid) neighbor
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
      Sxx = (feq[1] + feq[2] + feq[5] + feq[6] + feq[7] + feq[8] + cs2) / rhon - cs2;
      Syy = (feq[3] + feq[4] + feq[5] + feq[6] + feq[8] + feq[7] + cs2) / rhon - cs2;
      Sxy = (feq[5] + feq[6] - feq[7] - feq[8]) / rhon;
    } else if(type == CELL_TYPE.IG) {
      for(int i=1; i<9 ;i++){
        int idNx = mod(p_x+D2Q9_cx[i],p_simulation.getNx());
        int idNy = mod(p_y+D2Q9_cy[i],p_simulation.getNy());

        CELL_TYPE typeN = p_simulation.getCell(idNx,idNy).getType();
        if(typeN == CELL_TYPE.L) {
          p_simulation.getCell(idNx,idNy).setType(CELL_TYPE.I); // + bubble merge
        } else if(typeN == CELL_TYPE.IL) {
          p_simulation.getCell(idNx,idNy).setType(CELL_TYPE.I);
        }
      }
    }
  }
  
  public void surface3(int p_x, int p_y, LBM p_simulation) { 
    if(type==CELL_TYPE.SOLID) return;
    
    float rhon = rho;
    float massn = mass;
    float massexn = 0.f;
    float phin = 0.f;
    
    if (type==CELL_TYPE.L || type==CELL_TYPE.IL) {
      if (type==CELL_TYPE.IL) type = CELL_TYPE.L; // + bubble split
      massexn = massn - rhon; // dump mass-rho difference into excess mass
      massn = rhon; // fluid cell mass has to equal rho
      phin = 1.f;
    
      // use previous tag to record for the bubble volume change computation
      //mlflow[0].previous_tag[curind] = mlflow[0].tag_matrix[curind];
      //mlflow[0].tag_matrix[curind] = -1;
    }
    else if (type==CELL_TYPE.I || type==CELL_TYPE.GI) {
      if (type==CELL_TYPE.GI)  type = CELL_TYPE.I;
      massexn = massn > rhon ? massn - rhon : massn < 0.f ? massn : 0.f; // allow interface cells with mass>rho or mass<0
      massn = constrain(massn, 0.f, rhon);
      phin = rhon > 0.f ? massn / rhon : 0.5f; // calculate fill level for next step (only necessary for interface cells)
    }
    else if (type==CELL_TYPE.G || type==CELL_TYPE.IG) {
      if (type==CELL_TYPE.IG) {
        type = CELL_TYPE.G;
        rho = 1.f;
        ux = 0.f;
        uy = 0.f;
      }
      massexn = massn; // dump remaining mass into excess mass
      massn = 0.f;
      phin = 0.f;
    }
    
    int counter = 0; // count (fluid|interface) neighbors
    for (int i=1; i<9 ;i++) {
      int idNx = mod(p_x+D2Q9_cx[i],p_simulation.getNx());
      int idNy = mod(p_y+D2Q9_cy[i],p_simulation.getNy());
      
      CELL_TYPE typeN = p_simulation.getCell(idNx,idNy).getType();
      if(typeN == CELL_TYPE.L || typeN == CELL_TYPE.I || typeN == CELL_TYPE.IL || typeN == CELL_TYPE.GI) counter++;
    }
    
    mass = massn + ((counter > 0) ? 0.f : massexn); // if excess mass can't be distributed to neighboring interface or fluid cells, add it to local mass (ensure mass conservation)
    massex = ((counter > 0) ? massexn / (float)counter : 0.f); // divide excess mass up for all interface or fluid neighbors
    phi = phin;
  }
}
