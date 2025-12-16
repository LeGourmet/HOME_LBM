public enum CELL_TYPE { SOLID, FLUID, EQUILIBRIUM };

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
  
  // --------------------------------------------- DESTRUCTOR / CONSTRUCTOR ----------------------------------------------  
  public Cell(CELL_TYPE p_type, float p_rho, float p_ux, float p_uy) {
    this.type = p_type;
    
    this.rho = max(p_rho, 0.f);
    this.ux = p_ux;
    this.uy = p_uy;
    this.Sxx = (p_ux*p_ux - cs2);
    this.Syy = (p_uy*p_uy - cs2);
    this.Sxy = (p_ux*p_uy);  // Syx = Sxy; // symetric tensor
    
    this._rho   = this.rho;
    this._ux  = this.ux;
    this._uy  = this.uy;
    this._Sxx = this.Sxx;
    this._Syy = this.Syy;
    this._Sxy = this.Sxy;  // Syx = Sxy; // symetric tensor
  } 
  
  // ------------------------------------------------------ GETTERS ------------------------------------------------------
  public CELL_TYPE getType() { return this.type; }
  public float getDensity() { return this.rho; }
  public float getVelocityX() { return this.ux; }
  public float getVelocityY() { return this.uy; }
  public float getSxx() { return this.Sxx; }
  public float getSyy() { return this.Syy; }
  public float getSxy() { return this.Sxy; }
  
  // ------------------------------------------------------ SETTERS ------------------------------------------------------
  public void setType(CELL_TYPE p_type) { this.type = p_type; }
  public void setDensity(float p_density) { this.rho = max(0.f,p_density); }
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
                                            (D2Q9_cx[p_i]*D2Q9_cy[p_i]*D2Q9_cy[p_i]-D2Q9_cx[p_i]*cs2) * (p_Syy*p_ux+2.f*p_Sxy*p_uy-2.f*p_ux*p_uy*p_uy))/cs6) 
                          - 1.f);
  }
  
  // -------------------------------------------------- FUNCTIONS FLOW ---------------------------------------------------  
  public void flowStreamingCollision(int p_x, int p_y, LBM p_simulation) {
    if(type==CELL_TYPE.SOLID || type==CELL_TYPE.EQUILIBRIUM) return;
     
    // external forces = body forces : Newton second law of motion : F = M*G = M*A => g = A
    float fx = p_simulation.getForceX(p_x, p_y);
    float fy = p_simulation.getForceY(p_x, p_y);
        
    float nu = nu_fluid;
    
    // ---------------------------------- STREAMING ----------------------------------
    float[] fin = new float[9];
    fin[0] = computeFi(0, rho, ux, uy, Sxx, Syy, Sxy);
    
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
    }
    
    // ---------------------------------- RECONSTRUCTION ----------------------------------
    float rhoT =  1.f + fin[0] + fin[1] + fin[2] + fin[3] + fin[4] + fin[5] + fin[6] + fin[7] + fin[8];  // 1 + sum(fi);
    float uxT  = (fin[1] - fin[2] + fin[5] - fin[6] + fin[7] - fin[8]) / rhoT + 0.5f * fx;               //     sum(fi * cix);
    float uyT  = (fin[3] - fin[4] + fin[5] - fin[6] + fin[8] - fin[7]) / rhoT + 0.5f * fy;               //     sum(fi * ciy);
    float SxxT = (fin[1] + fin[2] + fin[5] + fin[6] + fin[7] + fin[8] + cs2) / rhoT - cs2;               //     sum(fi * (cix*cix-cs2));
    float SyyT = (fin[3] + fin[4] + fin[5] + fin[6] + fin[8] + fin[7] + cs2) / rhoT - cs2;               //     sum(fi * (ciy*ciy-cs2));
    float SxyT = (fin[5] + fin[6] - fin[7] - fin[8]) / rhoT;                                             //     sum(fi * (cix*ciy));
    
    // ---------------------------------- COLLISION ----------------------------------
    float uTNorm = sqrt(sq(uxT)+sq(uyT));
    if(uTNorm>cs){
      uxT *= cs/uTNorm;
      uyT *= cs/uTNorm;
    }
    
    float tau = 0.5f + nu/cs2;
    
    // compute moment at time t+1
    _rho = rhoT; 
    _ux = uxT + 0.5f*fx; // (Guo forcing, Krueger p.233f)
    _uy = uyT + 0.5f*fy; // (Guo forcing, Krueger p.233f)
    _Sxx = (tau-1.f)/(2.f*tau)*(SxxT-SyyT+uyT*uyT+fx*uxT-fy*uyT) + (tau+1.f)/(2.f*tau)*uxT*uxT + fx*uxT;
    _Syy = (tau-1.f)/(2.f*tau)*(SyyT-SxxT+uxT*uxT+fy*uyT-fx*uxT) + (tau+1.f)/(2.f*tau)*uyT*uyT + fy*uyT;
    _Sxy = (1.f-1.f/tau)*SxyT + uxT*uyT/tau + (2.f*tau-1.f)/(2.f*tau)*(fx*uyT+fy*uxT);
  }  
  
  public void flowCollision(int p_x, int p_y, LBM p_simulation) { 
    this.rho = this._rho;
    this.ux = this._ux;
    this.uy = this._uy;
    this.Sxx = this._Sxx;
    this.Syy = this._Syy;
    this.Sxy = this._Sxy;
  } 
}
