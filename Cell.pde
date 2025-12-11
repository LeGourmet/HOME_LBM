public enum CELL_TYPE { SOLID, FLUID, INTERFACE, GAS, EQUILIBRIUM };

public class Cell {
  // ----------------------------------------------------- ATTRIBUTS -----------------------------------------------------
  private CELL_TYPE type;
    
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
  private float mass;
  private float massex;
  
  // --------------------------------------------- DESTRUCTOR / CONSTRUCTOR ----------------------------------------------  
  public Cell(CELL_TYPE p_type, float p_p, float p_ux, float p_uy, float p_phi) {
    this.type = p_type;
    
    this.p   = max(p_p,0.f);
    this.ux  = p_ux;
    this.uy  = p_uy;
    this.Sxx = (p_ux*p_ux - cs2);
    this.Syy = (p_uy*p_uy - cs2);
    this.Sxy = (p_ux*p_uy);  // Syx = Sxy; // symetric tensor
    
    this._p   = this.p;
    this._ux  = this.ux;
    this._uy  = this.uy;
    this._Sxx = this.Sxx;
    this._Syy = this.Syy;
    this._Sxy = this.Sxy;  // Syx = Sxy; // symetric tensor
    
    this.phi = (p_phi<0.5f) ? 0.f : 1.f;
  } 
  
  // ------------------------------------------------------ GETTERS ------------------------------------------------------
  public CELL_TYPE getType() { return this.type; }
  public float getDensity() { return this.p; }
  public float getVelocityX() { return this.ux; }
  public float getVelocityY() { return this.uy; }
  public float getSxx() { return this.Sxx; }
  public float getSyy() { return this.Syy; }
  public float getSxy() { return this.Sxy; }
  public float getPhi() { return this.phi; }
  
  // ------------------------------------------------------ SETTERS ------------------------------------------------------
  public void setType(CELL_TYPE p_type) { this.type = p_type; }
  public void setDensity(float p_density) { this.p = max(0.f,p_density); }
  public void setVelocityX(float p_velocity) { this.ux = p_velocity; }
  public void setVelocityY(float p_velocity) { this.uy = p_velocity; }
  
  // -------------------------------------------------- FUNCTIONS DDFs ---------------------------------------------------  
  private float computeFiEq(int p_i, float p_p, float p_ux, float p_uy) {
    return p_p * D2Q9_w[p_i] * (1.f + (D2Q9_cx[p_i]*p_ux + D2Q9_cy[p_i]*p_uy)/cs2 + 0.5f*sq(D2Q9_cx[p_i]*p_ux + D2Q9_cy[p_i]*p_uy)/cs4 - 0.5f*(p_ux*p_ux + p_uy*p_uy)/cs2) - D2Q9_w[p_i];
  }
  
  private float computeFi(int p_i, float p_p, float p_ux, float p_uy, float p_Sxx, float p_Syy, float p_Sxy) {
    return p_p * D2Q9_w[p_i] * (1.f +
                                (D2Q9_cx[p_i]*p_ux + D2Q9_cy[p_i]*p_uy)/cs2 +
                                0.5f*( 2.f*p_Sxy*D2Q9_cx[p_i]*D2Q9_cy[p_i] + p_Sxx*(D2Q9_cx[p_i]*D2Q9_cx[p_i]-cs2) + p_Syy*(D2Q9_cy[p_i]*D2Q9_cy[p_i]-cs2))/cs4 +
                                0.5f*( (D2Q9_cx[p_i]*D2Q9_cx[p_i]*D2Q9_cy[p_i]-D2Q9_cy[p_i]*cs2) * (p_Sxx*p_uy+2.f*p_Sxy*p_ux-2.f*p_ux*p_ux*p_uy) +
                                       (D2Q9_cx[p_i]*D2Q9_cy[p_i]*D2Q9_cy[p_i]-D2Q9_cx[p_i]*cs2) * (p_Syy*p_ux+2.f*p_Sxy*p_uy-2.f*p_ux*p_uy*p_uy))/cs6) - D2Q9_w[p_i];
  }
  
  public void initialize(){
    if(type == CELL_TYPE.SOLID) 
      phi = 0.f;
    else if(phi == 1.f) 
      type = CELL_TYPE.FLUID;
     else
        type = CELL_TYPE.GAS;
  }
  
  // -------------------------------------------------- FUNCTIONS FLOW ---------------------------------------------------  
  public void flowStreamingCollision(int p_x, int p_y, LBM p_simulation) {
    CELL_TYPE type = p_simulation.getCell(p_x, p_y).getType();
    
    if(type==CELL_TYPE.SOLID || type==CELL_TYPE.EQUILIBRIUM) return;
    
    // ---------------------------------- STREAMING ----------------------------------
    
    float[] fin = new float[9];
    fin[0] = computeFi(0, p, ux, uy, Sxx, Syy, Sxy);
    
    for(int i=1; i<9 ;i++){
      int j = (i%2==0) ? i-1 : i+1;
      
      int idNx = mod(p_x+D2Q9_cx[j],p_simulation.getNx());
      int idNy = mod(p_y+D2Q9_cy[j],p_simulation.getNy());

      CELL_TYPE typeN = p_simulation.getCell(idNx,idNy).getType();
      float pN   = p_simulation.getCell(idNx,idNy).getDensity();
      float uxN  = p_simulation.getCell(idNx,idNy).getVelocityX();
      float uyN  = p_simulation.getCell(idNx,idNy).getVelocityY();
      float SxxN = p_simulation.getCell(idNx,idNy).getSxx();
      float SyyN = p_simulation.getCell(idNx,idNy).getSyy();
      float SxyN = p_simulation.getCell(idNx,idNy).getSxy();

      if(typeN==CELL_TYPE.SOLID)
        fin[i] = computeFi(i, p, uxN, uyN, Sxx+uxN*uxN-ux*ux, Syy+uyN*uyN-uy*uy, Sxy+uxN*uyN-ux*uy);
      else if(typeN==CELL_TYPE.EQUILIBRIUM)
        fin[i] = computeFiEq(i, pN, uxN, uyN);
      else
        fin[i] = computeFi(i, pN, uxN, uyN, SxxN, SyyN, SxyN);
    }

    // ---------------------------------- RECONSTRUCTION ----------------------------------
    float pT   = 1.f + fin[0] + fin[1] + fin[2] + fin[3] + fin[4] + fin[5] + fin[6] + fin[7] + fin[8];
    float uxT  = (fin[1] - fin[2] + fin[5] - fin[6] + fin[7] - fin[8]) / pT;
    float uyT  = (fin[3] - fin[4] + fin[5] - fin[6] + fin[8] - fin[7]) / pT;
    float SxxT = (fin[1] + fin[2] + fin[5] + fin[6] + fin[7] + fin[8] - cs2*(pT-1.f)) / pT;
    float SyyT = (fin[3] + fin[4] + fin[5] + fin[6] + fin[8] + fin[7] - cs2*(pT-1.f)) / pT;
    float SxyT = (fin[5] + fin[6] - fin[7] - fin[8]) / pT;
    
    // ---------------------------------- COLLISION ----------------------------------
    float nu = nu_fluid ;
    float tau = nu/cs2 + 0.5f;
    
    // external forces = body forces : Newton second law of motion : F = M*G = M*A => g = A
    float fx = p_simulation.getForceX(p_x, p_y);
    float fy = p_simulation.getForceY(p_x, p_y);
        
    // ajoute de la pression jusqu'au buffer overflow ! => trouver mieux
    float uTNorm = sqrt(sq(uxT+0.5f*fx)+sq(uyT+0.5f*fy));
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
  
  public void swapMoments() {     
    this.p = this._p;
    this.ux = this._ux;
    this.uy = this._uy;
    this.Sxx = this._Sxx;
    this.Syy = this._Syy;
    this.Sxy = this._Sxy;
  } 
}
