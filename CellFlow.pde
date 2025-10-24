public class CellFlow {
  // ----------------------------------------------------- ATTRIBUTS -----------------------------------------------------
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
  
  // --------------------------------------------- DESTRUCTOR / CONSTRUCTOR ----------------------------------------------  
  public CellFlow(float p_p, float p_ux, float p_uy) {
    this.p = p_p;
    this.ux = p_ux;
    this.uy = p_uy;
    this.Sxx = (p_ux*p_ux - cs2);
    this.Syy = (p_uy*p_uy - cs2);
    this.Sxy = (p_ux*p_uy);  // Syx = Sxy; // symetric tensor
  
    this._p = this.p;
    this._ux = this.ux;
    this._uy = this.uy;
    this._Sxx = this.Sxx;
    this._Syy = this.Syy;
    this._Sxy = this.Sxy;  // Syx = Sxy; // symetric tensor
  } 
  
  // ------------------------------------------------------ GETTERS ------------------------------------------------------
  public float getPressure() { return this.p; }
  public float getVelocityX() { return this.ux; }
  public float getOldVelocityX() { return this._ux; }
  public float getVelocityY() { return this.uy; }
  public float getOldVelocityY() { return this._uy; }
  public float getSxx() { return this.Sxx; }
  public float getSyy() { return this.Syy; }
  public float getSxy() { return this.Sxy; }
  
  // ------------------------------------------------------ SETTERS ------------------------------------------------------
  public void setPressure(float p_pressure) { this.p = max(0.f,p_pressure); }
  public void setVelocityX(float p_velocity) { this.ux = p_velocity; }
  public void setVelocityY(float p_velocity) { this.uy = p_velocity; }
  
  // -------------------------------------------------- FUNCTIONS VBDFs --------------------------------------------------  
  // DDF shifting : fi_eq_shift = fi_eq-wi 
  private float computeFiEq(int p_i, float p_p, float p_ux, float p_uy) {
    return D2Q9_w[p_i] * (p_p + (D2Q9_cx[p_i]*p_ux + D2Q9_cy[p_i]*p_uy)/cs2 + 0.5f*sq(D2Q9_cx[p_i]*p_ux + D2Q9_cy[p_i]*p_uy)/cs4 - 0.5f*(p_ux*p_ux + p_uy*p_uy)/cs2 - 1.f);
  }
  
  // DDF shifting : fi_shift = fi-wi 
  private float computeFi(int p_i, float p_p, float p_ux, float p_uy, float p_Sxx, float p_Syy, float p_Sxy) {
    return D2Q9_w[p_i] * (p_p +
                          (D2Q9_cx[p_i]*p_ux + D2Q9_cy[p_i]*p_uy)/cs2 +
                          0.5f*( 2.f*p_Sxy*D2Q9_cx[p_i]*D2Q9_cy[p_i] + p_Sxx*(D2Q9_cx[p_i]*D2Q9_cx[p_i]-cs2) + p_Syy*(D2Q9_cy[p_i]*D2Q9_cy[p_i]-cs2))/cs4 +
                          0.5f*( (D2Q9_cx[p_i]*D2Q9_cx[p_i]*D2Q9_cy[p_i]-D2Q9_cy[p_i]*cs2) * (p_Sxx*p_uy+2.f*p_Sxy*p_ux-2.f*p_ux*p_ux*p_uy) +
                                 (D2Q9_cx[p_i]*D2Q9_cy[p_i]*D2Q9_cy[p_i]-D2Q9_cx[p_i]*cs2) * (p_Syy*p_ux+2.f*p_Sxy*p_uy-2.f*p_ux*p_uy*p_uy))/cs6 
                          - 1.f);
  }
  
  // -------------------------------------------------- FUNCTIONS FLOW ---------------------------------------------------  
  // // DDF shifting : p = sum(fi)+1 
  public void streaming(int p_x, int p_y, LBM p_simulation) {
    CELL_TYPE type = p_simulation.getType(p_x, p_y);
    if(type==CELL_TYPE.SOLID || type==CELL_TYPE.EQUILIBRIUM) return;
        
    _p = 1.f;
    _ux = 0.f;
    _uy = 0.f;
    _Sxx = 0.f;
    _Syy = 0.f;
    _Sxy = 0.f;
    
    for(int i=0; i<9 ;i++){
      int j = (i==0) ? i : ((i%2==0) ? i-1 : i+1);
      
      int idNx = mod(p_x+D2Q9_cx[j],p_simulation.getNx());
      int idNy = mod(p_y+D2Q9_cy[j],p_simulation.getNy());

      CELL_TYPE typeN = p_simulation.getType(idNx,idNy);
      float pN = p_simulation.getPressure(idNx,idNy);
      float uxN = p_simulation.getVelocityX(idNx,idNy);
      float uyN = p_simulation.getVelocityY(idNx,idNy);
      float SxxN = p_simulation.getSxx(idNx,idNy);
      float SyyN = p_simulation.getSyy(idNx,idNy);
      float SxyN = p_simulation.getSxy(idNx,idNy);

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
      _Sxx += fi * (D2Q9_cx[i]*D2Q9_cx[i]-cs2);
      _Syy += fi * (D2Q9_cy[i]*D2Q9_cy[i]-cs2);
      _Sxy += fi * (D2Q9_cx[i]*D2Q9_cy[i]);
    }
  }
  
  public void collision(int p_x, int p_y, LBM p_simulation) { 
    CELL_TYPE type = p_simulation.getType(p_x, p_y);
    if(type==CELL_TYPE.SOLID || type==CELL_TYPE.EQUILIBRIUM) return;
    
    float phi = simulation.getPhi(p_x, p_y);
        
    float rho = (1-phi)*rho_air + phi*rho_fluid;         // max 1e-5f
    float nu = 1.f/( (1.f-phi)/nu_air + phi/nu_fluid );  // max 1e-5f
    float tau = 0.5f + nu/cs2;
        
    // unroll => wi*ci_alpha*phi(x+cix,y+ciy)/cs2
    float GphiX = (D2Q5_w[1]*p_simulation.getPhi(mod(p_x+1,p_simulation.getNx()),p_y) - D2Q5_w[2]*p_simulation.getPhi(mod(p_x-1,p_simulation.getNx()),p_y))/cs2;
    float GphiY = (D2Q5_w[3]*p_simulation.getPhi(p_x,mod(p_y+1,p_simulation.getNy())) - D2Q5_w[4]*p_simulation.getPhi(p_x,mod(p_y-1,p_simulation.getNy())))/cs2;
    /*float nGphi = sqrt(sq(GphiX)+sq(GphiY));
    if(nGphi > (0.6f * 4.f*phi*(1.f-phi)/interfacial_thickness)){
      GphiX = GphiX/nGphi * 0.6f;
      GphiY = GphiY/nGphi * 0.6f;
    }*/
    float GrhoX = (rho_fluid - rho_air) * GphiX;
    float GrhoY = (rho_fluid - rho_air) * GphiY;
    
    // unroll => 2*wi*(phi(x+cix,y+ciy)/cs2
    float GphiSQ= 2.f/cs * ( D2Q5_w[1]*(p_simulation.getPhi(mod(p_x+1,p_simulation.getNx()),p_y)-phi) + 
                             D2Q5_w[2]*(p_simulation.getPhi(mod(p_x-1,p_simulation.getNx()),p_y)-phi) + 
                             D2Q5_w[3]*(p_simulation.getPhi(p_x,mod(p_y+1,p_simulation.getNy()))-phi) + 
                             D2Q5_w[4]*(p_simulation.getPhi(p_x,mod(p_y-1,p_simulation.getNy()))-phi) );
   
    // body forces
    float fbx = p_simulation.getForceX(p_x, p_y);
    float fby = p_simulation.getForceY(p_x, p_y);
    
    // pressure force : -P * cs2 * gradient(rho) : -P * gradient(rho)
    float fpx = 0.0001f* -p * cs2 * GrhoX;
    float fpy = 0.0001f* -p * cs2 * GrhoY;
    
    // viscosity force : nu * [gradient(U) + transpose(grandient(u))] * gradient(rho)
    float fvx = 0.0001f* ( (ux*ux-Sxx) * GrhoX + (ux*uy-Sxy) * GrhoY);
    float fvy = 0.0001f* ( (uy*ux-Sxy) * GrhoX + (uy*uy-Syy) * GrhoY);
    
    // surface tension force : abs(ca_fluid - ca_air) < cs/100.f
    float fsx = 0.0001f * (ca_fluid + ca_air) * GphiX * (24.f/interfacial_thickness * sq(missibility) * phi * (1.f-phi) * (1.f-2.f*phi) - 3.f*interfacial_thickness/2.f * GphiSQ);
    float fsy = 0.0001f * (ca_fluid + ca_air) * GphiY * (24.f/interfacial_thickness * sq(missibility) * phi * (1.f-phi) * (1.f-2.f*phi) - 3.f*interfacial_thickness/2.f * GphiSQ);   
 
    // external forces : Newton second law of motion
    float fx = fbx + (fpx+fvx+fsx)/rho; 
    float fy = fby + (fpy+fvy+fsy)/rho;
    
    // ajoute de la pression jusqu'au buffer overflow ! => trouver mieux
    float _uNorm = sqrt(sq(_ux+0.5f*fx)+sq(_uy+0.5f*fy));
    if(_uNorm>cs){
      _ux *= cs/_uNorm;
      _uy *= cs/_uNorm;
    }
    
    float p_tmp = p;
    float ux_tmp = ux;
    float uy_tmp = uy;
    float Sxx_tmp = Sxx;
    float Syy_tmp = Syy;
    float Sxy_tmp = Sxy;
    
    // compute moment at time t+1
    p = _p; 
    ux = _ux + 0.5f*fx; // (Guo forcing, Krueger p.233f)
    uy = _uy + 0.5f*fy; // (Guo forcing, Krueger p.233f)
    Sxx = (tau-1.f)/(2.f*tau)*(_Sxx-_Syy+_uy*_uy+fx*_ux-fy*_uy) + (tau+1.f)/(2.f*tau)*_ux*_ux + fx*_ux;
    Syy = (tau-1.f)/(2.f*tau)*(_Syy-_Sxx+_ux*_ux+fy*_uy-fx*_ux) + (tau+1.f)/(2.f*tau)*_uy*_uy + fy*_uy;
    Sxy = (1.f-1.f/tau)*_Sxy + _ux*_uy/tau + (2.f*tau-1.f)/(2.f*tau)*(fx*_uy+fy*_ux);
    
    // save moments at time t
    _p = p_tmp;
    _ux = ux_tmp;
    _uy = uy_tmp;
    _Sxx = Sxx_tmp;
    _Syy = Syy_tmp;
    _Sxy = Sxy_tmp;
  }
}
