public enum CELL_TYPE { SOLID, LIQUID, INTERFACE, GAS };

public class Cell {
  // ----------------------------------------------------- ATTRIBUTS -----------------------------------------------------
  private CELL_TYPE type;
  private boolean eq;
    
  private float p;
  private float ux;
  private float uy;
  private float Sxx;
  private float Syy;
  private float Sxy;
  
  private float _p;
  private float _ux;
  private float _uy;
  private float _Sxx;
  private float _Syy;
  private float _Sxy;
  
  private float phi;
  protected float[] hi;
  
  // --------------------------------------------- DESTRUCTOR / CONSTRUCTOR ----------------------------------------------  
  public Cell(CELL_TYPE p_type, float p_p, float p_ux, float p_uy, float p_phi, boolean p_eq) {
    this.type = p_type;
    this.eq = p_eq;
    
    this.p = max(p_p, 0.f);
    this.ux = p_ux;
    this.uy = p_uy;
    this.Sxx = (p_ux*p_ux - cs2);
    this.Syy = (p_uy*p_uy - cs2);
    this.Sxy = (p_ux*p_uy);  // Syx = Sxy; // symetric tensor
    
    this._p   = this.p;
    this._ux  = this.ux;
    this._uy  = this.uy;
    this._Sxx = this.Sxx;
    this._Syy = this.Syy;
    this._Sxy = this.Sxy;  // Syx = Sxy; // symetric tensor
    
    this.phi = constrain(p_phi,0.f,1.f);
    
    this.hi = new float[5];
    for(int i=0; i<5 ;i++) this.hi[i] = computeHiEq(i, this.phi, this.ux, this.uy);
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
  public boolean isEQ() { return this.eq; }
  
  // ------------------------------------------------------ SETTERS ------------------------------------------------------
  public void setType(CELL_TYPE p_type) { this.type = p_type; }
  public void setPressure(float p_pressure) { this.p = max(0.f,p_pressure); }
  public void setVelocityX(float p_velocity) { this.ux = p_velocity; }
  public void setVelocityY(float p_velocity) { this.uy = p_velocity; }
  public void setPhi(float p_phi) {
    this.phi = constrain(p_phi,0.f,1.f);
    this.hi = new float[5];
    for(int i=0; i<5 ;i++) this.hi[i] = computeHiEq(i, this.phi, this.ux, this.uy);
  }
  public void setEq(boolean p_eq) { this.eq = p_eq; }
  
  // -------------------------------------------------- FUNCTIONS DDFs ---------------------------------------------------  
  private float computeFiEq(int p_i, float p_p, float p_ux, float p_uy) {
    return D2Q9_w[p_i] * (p_p + (D2Q9_cx[p_i]*p_ux + D2Q9_cy[p_i]*p_uy)/cs2 + 0.5f*sq(D2Q9_cx[p_i]*p_ux + D2Q9_cy[p_i]*p_uy)/cs4 - 0.5f*(p_ux*p_ux + p_uy*p_uy)/cs2 - 1.f);
  }
  
  private float computeFi(int p_i, float p_p, float p_ux, float p_uy, float p_Sxx, float p_Syy, float p_Sxy) {
    return D2Q9_w[p_i] * (p_p +
                          (D2Q9_cx[p_i]*p_ux + D2Q9_cy[p_i]*p_uy)/cs2 +
                          0.5f*( 2.f*p_Sxy*D2Q9_cx[p_i]*D2Q9_cy[p_i] + p_Sxx*(D2Q9_cx[p_i]*D2Q9_cx[p_i]-cs2) + p_Syy*(D2Q9_cy[p_i]*D2Q9_cy[p_i]-cs2))/cs4 +
                          0.5f*( (D2Q9_cx[p_i]*D2Q9_cx[p_i]*D2Q9_cy[p_i]-D2Q9_cy[p_i]*cs2) * (p_Sxx*p_uy+2.f*p_Sxy*p_ux-2.f*p_ux*p_ux*p_uy) +
                                 (D2Q9_cy[p_i]*D2Q9_cy[p_i]*D2Q9_cx[p_i]-D2Q9_cx[p_i]*cs2) * (p_Syy*p_ux+2.f*p_Sxy*p_uy-2.f*p_uy*p_uy*p_ux))/cs6
                          - 1.f);
  }
  
  private float computeHiEq(int p_i, float p_phi, float p_ux, float p_uy) {
    return D2Q5_w[p_i] * p_phi * (1.f + (D2Q5_cx[p_i]*p_ux + D2Q5_cy[p_i]*p_uy)/cs2);
  }
  
  // -------------------------------------------------- FUNCTIONS FLOW ---------------------------------------------------  
  public void flowStreamingCollision(int p_x, int p_y, LBM p_simulation) {
    if(type==CELL_TYPE.SOLID || eq) return;
        
    float rho = phi*rho_fluid + (1.f-phi)*rho_air;
    float nu = 1.f/(phi/nu_fluid + (1.f-phi)/nu_air);
    float tau = 0.5f + nu/cs2;
        
    // unroll => wi*ci_alpha*phi(x+cix,y+ciy)/cs2
    float GphiX = (D2Q5_w[1]*p_simulation.getCell(mod(p_x+1,p_simulation.getNx()),p_y).getPhi() - D2Q5_w[2]*p_simulation.getCell(mod(p_x-1,p_simulation.getNx()),p_y).getPhi())/cs2;
    float GphiY = (D2Q5_w[3]*p_simulation.getCell(p_x,mod(p_y+1,p_simulation.getNy())).getPhi() - D2Q5_w[4]*p_simulation.getCell(p_x,mod(p_y-1,p_simulation.getNy())).getPhi())/cs2;
    float GrhoX = (rho_fluid - rho_air) * GphiX;
    float GrhoY = (rho_fluid - rho_air) * GphiY;
    // unroll => 2*wi*(phi(x+cix,y+ciy)/cs2
    float GphiSQ= 2.f/cs * ( D2Q5_w[1]*(p_simulation.getCell(mod(p_x+1,p_simulation.getNx()),p_y).getPhi()-phi) + 
                             D2Q5_w[2]*(p_simulation.getCell(mod(p_x-1,p_simulation.getNx()),p_y).getPhi()-phi) + 
                             D2Q5_w[3]*(p_simulation.getCell(p_x,mod(p_y+1,p_simulation.getNy())).getPhi()-phi) + 
                             D2Q5_w[4]*(p_simulation.getCell(p_x,mod(p_y-1,p_simulation.getNy())).getPhi()-phi) );
    
    // body forces
    float fbx = p_simulation.getForceX(p_x, p_y);
    float fby = p_simulation.getForceY(p_x, p_y);
    
    // pressure force : -P * cs2 * gradient(rho) :|: -P * gradient(rho)
    float fpx = 0.1f * -p * cs2 * GrhoX;
    float fpy = 0.1f * -p * cs2 * GrhoY;
    
    // viscosity force : nu * [gradient(U) + transpose(grandient(u))] * gradient(rho)
    float fvx = ( (ux*ux-Sxx) * GrhoX + (ux*uy-Sxy) * GrhoY);
    float fvy = ( (uy*ux-Sxy) * GrhoX + (uy*uy-Syy) * GrhoY);
    
    // surface tension force : abs(ca_fluid - ca_air) < cs/100.f
    float fsx = 0.1f * (ca_fluid + ca_air) * GphiX * (24.f/interfacial_thickness * sq(missibility) * phi * (1.f-phi) * (1.f-2.f*phi) - 3.f*interfacial_thickness/2.f * GphiSQ);
    float fsy = 0.1f * (ca_fluid + ca_air) * GphiY * (24.f/interfacial_thickness * sq(missibility) * phi * (1.f-phi) * (1.f-2.f*phi) - 3.f*interfacial_thickness/2.f * GphiSQ);
        
    // external forces : Newton second law of motion
    float fx = fbx + (fpx+fvx+fsx)/rho; 
    float fy = fby + (fpy+fvy+fsy)/rho;
        
    // ---------------------------------- STREAMING ----------------------------------
    float[] fin = new float[9];
    fin[0] = computeFi(0, p, ux, uy, Sxx, Syy, Sxy);
    
    for(int i=1; i<9 ;i++){
      int j = (i%2==0) ? i-1 : i+1;
      
      int idNx = mod(p_x+D2Q9_cx[j],p_simulation.getNx());
      int idNy = mod(p_y+D2Q9_cy[j],p_simulation.getNy());

      CELL_TYPE typeN = p_simulation.getCell(idNx,idNy).getType();
      boolean eqN = p_simulation.getCell(idNx,idNy).isEQ();
      float pN   = p_simulation.getCell(idNx,idNy).getPressure();
      float uxN  = p_simulation.getCell(idNx,idNy).getVelocityX();
      float uyN  = p_simulation.getCell(idNx,idNy).getVelocityY();
      float SxxN = p_simulation.getCell(idNx,idNy).getSxx();
      float SyyN = p_simulation.getCell(idNx,idNy).getSyy();
      float SxyN = p_simulation.getCell(idNx,idNy).getSxy();

      if(typeN==CELL_TYPE.SOLID)
        fin[i] = computeFi(i, p, uxN, uyN, Sxx+uxN*uxN-ux*ux, Syy+uyN*uyN-uy*uy, Sxy+uxN*uyN-ux*uy);
      else if(eqN)
        fin[i] = computeFiEq(i, pN, uxN, uyN);
      else
        fin[i] = computeFi(i, pN, uxN, uyN, SxxN, SyyN, SxyN);
    }
    
    // ---------------------------------- RECONSTRUCTION ----------------------------------
    float pT   =  1.f + fin[0] + fin[1] + fin[2] + fin[3] + fin[4] + fin[5] + fin[6] + fin[7] + fin[8];  // 1 + sum(fi);
    float uxT  = (fin[1] - fin[2] + fin[5] - fin[6] + fin[7] - fin[8]) + 0.5f * fx;                      //     sum(fi * cix);
    float uyT  = (fin[3] - fin[4] + fin[5] - fin[6] + fin[8] - fin[7]) + 0.5f * fy;                      //     sum(fi * ciy);
    float SxxT = (fin[1] + fin[2] + fin[5] + fin[6] + fin[7] + fin[8]) - (pT-1.f)*cs2;                   //     sum(fi * (cix*cix-cs2));
    float SyyT = (fin[3] + fin[4] + fin[5] + fin[6] + fin[8] + fin[7]) - (pT-1.f)*cs2;                   //     sum(fi * (ciy*ciy-cs2));
    float SxyT = (fin[5] + fin[6] - fin[7] - fin[8]);                                                    //     sum(fi * (cix*ciy));
    
    // ---------------------------------- COLLISION ----------------------------------
    float uTNorm = sqrt(sq(uxT)+sq(uyT));
    if(uTNorm>cs){
      uxT *= cs/uTNorm;
      uyT *= cs/uTNorm;
    }
    
    // compute moment at time t+1
    _p = pT; 
    _ux = uxT + 0.5f*fx; // (Guo forcing, Krueger p.233f)
    _uy = uyT + 0.5f*fy; // (Guo forcing, Krueger p.233f)
    _Sxx = (tau-1.f)/(2.f*tau)*(SxxT-SyyT+uyT*uyT+fx*uxT-fy*uyT) + (tau+1.f)/(2.f*tau)*uxT*uxT + fx*uxT;
    _Syy = (tau-1.f)/(2.f*tau)*(SyyT-SxxT+uxT*uxT+fy*uyT-fx*uxT) + (tau+1.f)/(2.f*tau)*uyT*uyT + fy*uyT;
    _Sxy = (1.f-1.f/tau)*SxyT + uxT*uyT/tau + (2.f*tau-1.f)/(2.f*tau)*(fx*uyT+fy*uxT);
  }
  
   public void swapMoments(int p_x, int p_y, LBM p_simulation) { 
    this.p = this._p;
    this.ux = this._ux;
    this.uy = this._uy;
    this.Sxx = this._Sxx;
    this.Syy = this._Syy;
    this.Sxy = this._Sxy;
  } 
  
  // -------------------------------------------------- FUNCTIONS PHASE ---------------------------------------------------  
  public void phaseCollision(int p_x, int p_y, LBM p_simulation){
    if(type==CELL_TYPE.SOLID || eq) return;
    
    float mo = 1.f/(phi/mo_fluid + (1.f-phi)/mo_air);
    float tau = 0.5f + mo/cs2;
        
    float GphiX = (D2Q5_w[1]*p_simulation.getCell(mod(p_x+1,p_simulation.getNx()),p_y).phi - D2Q5_w[2]*p_simulation.getCell(mod(p_x-1,p_simulation.getNx()),p_y).phi)/cs2;
    float GphiY = (D2Q5_w[3]*p_simulation.getCell(p_x,mod(p_y+1,p_simulation.getNy())).phi - D2Q5_w[4]*p_simulation.getCell(p_x,mod(p_y-1,p_simulation.getNy())).phi)/cs2;
    float nGphi = sqrt(sq(GphiX)+sq(GphiY));
    float nGphiX = (nGphi==0.f) ? 0.f : GphiX/nGphi;
    float nGphiY = (nGphi==0.f) ? 0.f : GphiY/nGphi;
        
    hi[0] += (computeHiEq(0,phi,ux,uy)-hi[0])/tau;
    hi[1] += (computeHiEq(1,phi,ux,uy)-hi[1])/tau + D2Q5_w[1] * 4.f/interfacial_thickness * phi*(1.f-phi) * missibility * nGphiX;
    hi[2] += (computeHiEq(2,phi,ux,uy)-hi[2])/tau - D2Q5_w[2] * 4.f/interfacial_thickness * phi*(1.f-phi) * missibility * nGphiX;
    hi[3] += (computeHiEq(3,phi,ux,uy)-hi[3])/tau + D2Q5_w[3] * 4.f/interfacial_thickness * phi*(1.f-phi) * missibility * nGphiY;
    hi[4] += (computeHiEq(4,phi,ux,uy)-hi[4])/tau - D2Q5_w[4] * 4.f/interfacial_thickness * phi*(1.f-phi) * missibility * nGphiY;
  }
  
  public void phaseStreaming(int p_x, int p_y, LBM p_simulation){
    if(type==CELL_TYPE.SOLID || eq) return;
    
    phi = 0.f;
    for(int i=0; i<5 ;i++){
      int j = (i==0) ? i : ((i%2==0) ? i-1 : i+1);
      Cell cell = p_simulation.getCell(mod(p_x+D2Q5_cx[j],p_simulation.getNx()),mod(p_y+D2Q5_cy[j],p_simulation.getNy())); 
      phi += (cell.type == CELL_TYPE.SOLID) ? hi[j] : cell.hi[i];
    }
  }
}
