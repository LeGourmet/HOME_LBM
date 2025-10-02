// should use 4 time finer discreate field for phi (2 per directions)
// should update hi at half time fi(t) => fi(t+1) ; hi(t) => hi(t+0.5) => hi(t+1)
// should update solid for moving solid boundary treatment

public enum CELL_TYPE { SOLID, FLUID, EQUILIBRIUM };

public class Cell {
  // ----------------------------------------------------- ATTRIBUTS -----------------------------------------------------
  private CELL_TYPE type;
    
  private float p;
  private float ux;
  private float uy;
  private float Sxx;
  private float Syy;
  private float Sxy;
  private float phi;
  
  private float _p = 0.f;
  private float _ux = 0.f;
  private float _uy = 0.f;
  private float _Sxx = 0.f;
  private float _Syy = 0.f;
  private float _Sxy = 0.f;
  
  protected float[] hi;
  
  // --------------------------------------------- DESTRUCTOR / CONSTRUCTOR ----------------------------------------------  
  public Cell(CELL_TYPE p_type, float p_p, float p_ux, float p_uy, float p_phi) {
    this.type = p_type;
    
    this.p = p_p;
    this.ux = p_ux;
    this.uy = p_uy;
    this.Sxx = (p_ux*p_ux - cs2);
    this.Syy = (p_uy*p_uy - cs2);
    this.Sxy = (p_ux*p_uy);  // Syx = Sxy; // symetric tensor
    this.phi = p_phi;
    
    this.hi = new float[5];
    for(int i=0; i<5 ;i++) this.hi[i] = computeHiEq(i,p_phi,p_ux,p_uy);
  } 
  
  // ------------------------------------------------------ GETTERS ------------------------------------------------------
  public CELL_TYPE getType() { return this.type; }
  public float getPressure() { return this.p; }
  public float getVelocityX() { return this.ux; }
  public float getVelocityY() { return this.uy; }
  public float getSxx() { return this.Sxx; }
  public float getSyy() { return this.Syy; }
  public float getSxy() { return this.Sxy; }
  public float getPhi() { return this.phi; }
  
  // ------------------------------------------------------ SETTERS ------------------------------------------------------
  public void setType(CELL_TYPE p_type) { this.type = p_type; }
  public void setPressure(float p_pressure) { this.p = max(0.f,p_pressure); }
  public void setVelocityX(float p_velocity) { this.ux = p_velocity; }
  public void setVelocityY(float p_velocity) { this.uy = p_velocity; }
  public void setPhi(float p_phi) { this.phi = constrain(p_phi,0.f,1.f); }
  
  // -------------------------------------------------- FUNCTIONS DDFs ---------------------------------------------------  
  // DDF shifting : fi_eq_shift = fi_eq-wi 
  private float computeFiEq(int p_i, float p_p, float p_ux, float p_uy) {
    return D2Q9_w[p_i] * (p_p + (D2Q9_cx[p_i]*p_ux + D2Q9_cy[p_i]*p_uy)/cs2 + 0.5f*sq(D2Q9_cx[p_i]*p_ux + D2Q9_cy[p_i]*p_uy)/cs4 - 0.5f*(p_ux*p_ux + p_uy*p_uy)/cs2 - 1.f);
  }
  
  // DDF shifting : fi_shift = fi-wi 
  private float computeFi(int p_i, float p_p, float p_ux, float p_uy, float p_Sxx, float p_Syy, float p_Sxy) {
    return D2Q9_w[p_i] * (p_p +
                          (D2Q9_cx[p_i]*p_ux + D2Q9_cy[p_i]*p_uy)/cs2 +
                          0.5f*( 2.f*p_Sxy*D2Q9_cx[p_i]*D2Q9_cy[p_i] + p_Sxx*(D2Q9_cx[p_i]*D2Q9_cx[p_i]-1.f/3.f) + p_Syy*(D2Q9_cy[p_i]*D2Q9_cy[p_i]-1.f/3.f))/cs4 +
                          0.5f*( (D2Q9_cx[p_i]*D2Q9_cx[p_i]*D2Q9_cy[p_i]-D2Q9_cy[p_i]*1.f/3.f) * (p_Sxx*p_uy+2.f*p_Sxy*p_ux-2.f*p_ux*p_ux*p_uy) +
                                 (D2Q9_cx[p_i]*D2Q9_cy[p_i]*D2Q9_cy[p_i]-D2Q9_cx[p_i]*1.f/3.f) * (p_Syy*p_ux+2.f*p_Sxy*p_uy-2.f*p_ux*p_uy*p_uy))/cs6 
                          - 1.f);
  }
  
  private float computeHiEq(int p_i, float p_phi, float p_ux, float p_uy) {
    return D2Q5_w[p_i] * p_phi * (1.f + (D2Q5_cx[p_i]*p_ux + D2Q5_cy[p_i]*p_uy)/cs2);
  }
  
  private float computeHi(int p_i, float p_phi, float p_ux, float p_uy) {
    return computeHiEq(p_i, p_phi, p_ux, p_uy);
  }
  
  // -------------------------------------------------- FUNCTIONS FLOW ---------------------------------------------------  
  // // DDF shifting : p = sum(fi)+1 
  public void flowStreaming(int p_idX, int p_idY, LBM p_simulation) {
    if(type==CELL_TYPE.SOLID || type==CELL_TYPE.EQUILIBRIUM) return;
        
    _p = 1.f;
    _ux = 0.f;
    _uy = 0.f;
    _Sxx = 0.f;
    _Syy = 0.f;
    _Sxy = 0.f;
    
    for(int i=0; i<9 ;i++){
      int j = (i==0) ? i : ((i%2==0) ? i-1 : i+1);
      
      int idNx = mod(p_idX+D2Q9_cx[j],p_simulation.getNx());
      int idNy = mod(p_idY+D2Q9_cy[j],p_simulation.getNy());

      CELL_TYPE typeN = p_simulation.getCell(idNx,idNy).getType();
      float pN = p_simulation.getCell(idNx,idNy).getPressure();
      float uxN = p_simulation.getCell(idNx,idNy).getVelocityX();
      float uyN = p_simulation.getCell(idNx,idNy).getVelocityY();
      float SxxN = p_simulation.getCell(idNx,idNy).getSxx();
      float SyyN = p_simulation.getCell(idNx,idNy).getSyy();
      float SxyN = p_simulation.getCell(idNx,idNy).getSxy();

      float fi;
      if(typeN==CELL_TYPE.SOLID){
        fi = computeFi(i, p, uxN, uyN, Sxx+uxN*uxN-ux*ux, Syy+uyN*uyN-uy*uy, Sxy+uxN*uyN-ux*uy);
      } else if(typeN==CELL_TYPE.EQUILIBRIUM){
        fi = computeFiEq(i, pN, uxN, uyN);
      } else { 
        fi = computeFi(i, pN, uxN, uyN, SxxN, SyyN, SxyN);
      }
      
      _p += fi;
      _ux += fi * D2Q9_cx[i];
      _uy += fi * D2Q9_cy[i];
      _Sxx += fi * (D2Q9_cx[i]*D2Q9_cx[i]-1.f/3.f);
      _Syy += fi * (D2Q9_cy[i]*D2Q9_cy[i]-1.f/3.f);
      _Sxy += fi * (D2Q9_cx[i]*D2Q9_cy[i]);
    }
  }
  
  public void flowCollision(int p_idX, int p_idY, LBM p_simulation) {    
    if(type==CELL_TYPE.SOLID || type==CELL_TYPE.EQUILIBRIUM) return;
    
    float rho = max(1e-5f,rho_air + phi * (rho_fluid - rho_air));
    float nu = 1.f/( (1.f-phi)/nu_air + phi/nu_fluid );
    float tau = 0.5f + nu/cs2;
    
    // unroll => wi*ci_alpha*phi(x+cix,y+ciy)/cs2
    float GphiX = (D2Q5_w[1]*p_simulation.getCell(mod(p_idX+1,p_simulation.getNx()),p_idY).getPhi() - D2Q5_w[2]*p_simulation.getCell(mod(p_idX-1,p_simulation.getNx()),p_idY).getPhi())/cs2;
    float GphiY = (D2Q5_w[3]*p_simulation.getCell(p_idX,mod(p_idY+1,p_simulation.getNy())).getPhi() - D2Q5_w[4]*p_simulation.getCell(p_idX,mod(p_idY-1,p_simulation.getNy())).getPhi())/cs2;
    float GrhoX = (rho_fluid - rho_air) * GphiX;
    float GrhoY = (rho_fluid - rho_air) * GphiY;
    // unroll => 2*wi*(phi(x+cix,y+ciy)/cs2
    float GphiSQ= 2.f/cs * ( D2Q5_w[1]*(p_simulation.getCell(mod(p_idX+1,p_simulation.getNx()),p_idY).getPhi()-phi) + 
                             D2Q5_w[2]*(p_simulation.getCell(mod(p_idX-1,p_simulation.getNx()),p_idY).getPhi()-phi) + 
                             D2Q5_w[3]*(p_simulation.getCell(p_idX,mod(p_idY+1,p_simulation.getNy())).getPhi()-phi) + 
                             D2Q5_w[4]*(p_simulation.getCell(p_idX,mod(p_idY-1,p_simulation.getNy())).getPhi()-phi) );
                             
    // external forces
    float fx = 0.f;
    float fy = 0.f;
    
    // body forces
    fx += p_simulation.getForceX(p_idX, p_idY);
    fy += p_simulation.getForceY(p_idX, p_idY);
    
    // pressure force
    fx += -p * GrhoX; //fx += -_p * GrhoX;
    fy += -p * GrhoY; //fy += -_p * GrhoY;
    
    // viscosity force
    fx += (ux*ux-Sxx) * GrhoX + (ux*uy-Sxy) * GrhoY; //fx += (_ux*_ux-_Sxx) * GrhoX + (_ux*_uy-_Sxy) * GrhoY;
    fy += (uy*ux-Sxy) * GrhoX + (uy*uy-Syy) * GrhoY; //fy += (_uy*_ux-_Sxy) * GrhoX + (_uy*_uy-_Syy) * GrhoY;
    
    // surface tension force
    float surfaceTensionForceTmp = (ca_air+ca_fluid) * (24.f/interfacial_thickness * phi * (1.f-phi) * (1.f-2.f*phi) - 3.f*interfacial_thickness/2.f * GphiSQ);
    fx += surfaceTensionForceTmp * GphiX;
    fy += surfaceTensionForceTmp * GphiY;
    
    // Sum_forces = Mass * A
    fx /= rho; 
    fy /= rho;
    
    // float F = D2Q9_w[i] * (1.f-1.f/(2.f*tau)) * (((D2Q9_cx[i]-ux)/cs2 + ((D2Q9_cx[i]*ux+D2Q9_cy[i]*uy)*D2Q9_cx[i])/cs4)*f.x + ((D2Q9_cy[i]-uy)/cs2 + ((D2Q9_cx[i]*ux+D2Q9_cy[i]*uy)*D2Q9_cy[i])/cs4)*f.y) ;
        
    // ajoute de la pression jusqu'au buffer overflow ! => trouver mieux
    float _uNorm = sqrt(sq(_ux+0.5f*fx)+sq(_uy+0.5f*fy));
    if(_uNorm>cs){
      _ux *= cs/_uNorm;
      _uy *= cs/_uNorm;
    }
    
    p = _p; 
    ux = _ux + 0.5f*fx; // (Guo forcing, Krueger p.233f) => manque la moitier
    uy = _uy + 0.5f*fy; // (Guo forcing, Krueger p.233f) => manque la moitier
    Sxx = (tau-1.f)/(2.f*tau)*(_Sxx-_Syy+_uy*_uy+fx*_ux-fy*_uy) + (tau+1.f)/(2.f*tau)*_ux*_ux + fx*_ux;
    Syy = (tau-1.f)/(2.f*tau)*(_Syy-_Sxx+_ux*_ux+fy*_uy-fx*_ux) + (tau+1.f)/(2.f*tau)*_uy*_uy + fy*_uy;
    Sxy = (1.f-1.f/tau)*_Sxy + _ux*_uy/tau + (2.f*tau-1.f)/(2.f*tau)*(fx*_uy+fy*_ux);
  }
  
  // -------------------------------------------------- FUNCTIONS PHASE -------------------------------------------------- 
  public void phaseStreaming(int p_idX, int p_idY, LBM p_simulation){
    /*for(int i=0; i<5 ;i++){
      int j = (i==0) ? i : ((i%2==0) ? i-1 : i+1);
      
      int idNx = (p_idX+D2Q5_cx[j]+p_simulation.getNx())%p_simulation.getNx();
      int idNy = (p_idY+D2Q5_cy[j]+p_simulation.getNy())%p_simulation.getNy();
      
      CELL_TYPE typeN = p_simulation.getCell(idNx,idNy).getType();
      float phiN = p_simulation.getCell(idNx,idNy).getPhi();
      float uxN = p_simulation.getCell(idNx,idNy).getVelocityX();
      float uyN = p_simulation.getCell(idNx,idNy).getVelocityY();
      
      if(typeN==CELL_TYPE.SOLID){
        _phi +=  computeHi(i, phi, uxN, uyN);
      } else if(typeN==CELL_TYPE.EQUILIBRIUM){
        _phi += computeHiEq(i, phiN, uxN, uyN);
      } else { 
        _phi += computeHi(i, phiN, uxN, uyN);
      }
    }*/
  }
  
  public void phaseCollision(int p_idX, int p_idY, LBM p_simulation){
    phi = 0.f;
    
    for(int i=0; i<5 ;i++){
      int j = (i==0) ? i : ((i%2==0) ? i-1 : i+1);
      int idNx = mod(p_idX+D2Q5_cx[j],p_simulation.getNx());
      int idNy = mod(p_idY+D2Q5_cy[j],p_simulation.getNy());
      //phi += hi[]; // => it or that 
    }
  }
}
