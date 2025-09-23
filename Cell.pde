public enum CELL_TYPE { SOLID, FLUID, EQUILIBRIUM };

public class Cell {
  // ----------------------------------------------------- ATTRIBUTS -----------------------------------------------------
  private CELL_TYPE flag;
  
  private float Fx;
  private float Fy;
  
  private float p;
  private float ux;
  private float uy;
  private float Sxx;
  private float Syy;
  private float Sxy;
  
  private float[] fi = new float[9]; 
  
  // --------------------------------------------- DESTRUCTOR / CONSTRUCTOR ----------------------------------------------  
  public Cell(CELL_TYPE p_flag, float p_Fx, float p_Fy, float p_p, float p_ux, float p_uy) {
    this.flag = p_flag;
    
    if(flag == CELL_TYPE.SOLID) return;
    
    this.Fx = p_Fx;
    this.Fy = p_Fy;
    
    this.p = p_p;
    this.ux = p_ux;
    this.uy = p_uy;
    this.Sxx = (ux*ux - cs2);
    this.Syy = (uy*uy - cs2);
    this.Sxy = (ux*uy);  // Syx = Sxy; // symetric tensor
    
    for(int i=0; i<9 ;i++)
      this.fi[i] = computeFiEq(i, p, ux, uy);
  } 
  
  // ------------------------------------------------------ GETTERS ------------------------------------------------------
  public CELL_TYPE getType() { return this.flag; } 
  public float getForceX() { return this.Fx; }
  public float getForceY() { return this.Fy; }
  public float getPressure() { return this.p; }
  public float getVelocityX() { return this.ux; }
  public float getVelocityY() { return this.uy; }
  
  public int getColor() {
    if (flag == CELL_TYPE.SOLID) return color(0);
  
    float val = constrain(sqrt(ux*ux + uy*uy)*2.5f, 0.f, 1.f);
    //float val = constrain(p/3.f, 0.f, 1.f);
        
    int palette[] = { 
      color(70, 70, 219),
      color(0, 255, 91), 
      color(0, 128, 0),
      color(255, 255, 0),
      color(255, 96, 0),
      color(107, 0, 0),
      color(223, 77, 77)
    };
    
    float x = constrain(val,0.00001f,0.99999f)*(palette.length-1.f);
    int idx = (int)floor(x);
    return lerpColor(palette[idx],palette[idx+1],  x-(float)(idx));
  }
  
  protected float getFi(int p_i) {
    if(p_i<0 || p_i>=9) return 0.f; 
    return this.fi[p_i]; 
  }
  
  // ------------------------------------------------------ SETTERS ------------------------------------------------------
  public void setType(CELL_TYPE p_type) { this.flag = p_type; } 
  public void setForceX(float p_force ) { this.Fx = p_force; }
  public void setForceY(float p_force) { this.Fy = p_force; }
  public void setPressure(float p_pressure) { this.p = p_pressure; }
  public void setVelocityX(float p_velocity) { this.ux = p_velocity; }
  public void setVelocityY(float p_velocity) { this.uy = p_velocity; }
  
  // ----------------------------------------------------- FUNCTIONS -----------------------------------------------------
  private float computeFiEq(int p_i, float p_p, float p_ux, float p_uy) {
    return w[p_i] * (p_p + (cx[p_i]*p_ux + cy[p_i]*p_uy)/cs2 + 0.5f*sq(cx[p_i]*p_ux + cy[p_i]*p_uy)/cs4 - 0.5f*(p_ux*p_ux + p_uy*p_uy)/cs2);
  }
  
  private float computeFi(int p_i, float p_p, float p_ux, float p_uy, float p_Sxx, float p_Syy, float p_Sxy) {
    return w[p_i] * (p_p +
                     (cx[p_i]*p_ux + cy[p_i]*p_uy)/cs2 +
                     0.5f*( 2.f*p_Sxy*cx[p_i]*cy[p_i] + p_Sxx*(cx[p_i]*cx[p_i]-1.f/3.f) + p_Syy*(cy[p_i]*cy[p_i]-1.f/3.f))/cs4 +
                     0.5f*( (cx[p_i]*cx[p_i]*cy[p_i]-cy[p_i]*1.f/3.f) * (p_Sxx*p_uy+2.f*p_Sxy*p_ux-2.f*p_ux*p_ux*p_uy) +
                            (cx[p_i]*cy[p_i]*cy[p_i]-cx[p_i]*1.f/3.f) * (p_Syy*p_ux+2.f*p_Sxy*p_uy-2.f*p_ux*p_uy*p_uy))/cs6);
  }
  
  public void reconstructDDFs() {
    if(flag != CELL_TYPE.FLUID) return;
    
    for(int i=0; i<9 ;i++)
      fi[i] = computeFi(i, p, ux, uy, Sxx, Syy, Sxy);
  }
  
  public void streamCollide(int idX, int idY, int Nx, int Ny, Cell[][] cells, float p_nu, float p_rho, float p_fx, float p_fy) {    
    if(flag != CELL_TYPE.FLUID) return;
    
    float _fx=p_fx+Fx, _fy=p_fy+Fy;
    // pressure force = 0
    // viscosity force = 0
    // surface tension force = 0
    
    float _p = 0.f; // temporary pressure
    float _ux = 0.f, _uy = 0.f; // temporary velocity
    float _Sxx = 0.f, _Syy = 0.f, _Sxy = 0.f; // temporary stress tensor
    float _Hxx = 0.f, _Hyy = 0.f, _Hxy = 0.f; // non-equilibrium stress tensor
   
    for(int i=0; i<9 ;i++){
      int j = (i==0) ? i : ((i%2==0) ? i-1 : i+1);
      
      int idNx = (idX+cx[j]+Nx)%Nx;
      int idNy = (idY+cy[j]+Ny)%Ny;

      float _fi, _fiEq;
      if(cells[idNx][idNy].getType()==CELL_TYPE.SOLID){
        float uxN = cells[idNx][idNy].getVelocityX();
        float uyN = cells[idNx][idNy].getVelocityY();
        _fi = computeFi(i, p, uxN, uyN, Sxx+uxN*uxN-ux*ux, Syy+uyN*uyN-uy*uy, Sxy+uxN*uyN-ux*uy); //or _fi = fi[j]; for bounce-back
        _fiEq = computeFiEq(i, p, uxN, uyN);
      } else { 
        _fi = cells[idNx][idNy].getFi(i);
        _fiEq = computeFiEq(i, p, ux, uy);  
      }
      
      _p += _fi;
      _ux += _fi * cx[i];
      _uy += _fi * cy[i];
      _Sxx += _fi * (cx[i]*cx[i]-1.f/3.f);
      _Syy += _fi * (cy[i]*cy[i]-1.f/3.f);
      _Sxy += _fi * (cx[i]*cy[i]);
      
      _Hxx += (_fi-_fiEq) * (cx[i]*cx[i]-1.f/3.f);
      _Hyy += (_fi-_fiEq) * (cy[i]*cy[i]-1.f/3.f);
      _Hxy += (_fi-_fiEq) * (cx[i]*cy[i]);
    }
    
    // (Guo forcing, Krueger p.233f) for volume force
    _fx /= max(1e-3f,_p);
    _fy /= max(1e-3f,_p);
    _ux += 0.5f*_fx;
    _uy += 0.5f*_fy;
    
    // moving boundary => // apply Dirichlet velocity boundaries if necessary (Krueger p.180, rho_solid=1)
    // => calculate_forcing_terms(uxn, uyn, uzn, fxn, fyn, fzn, Fin); // calculate volume force terms Fin from velocity field (Guo forcing, Krueger p.233f)
    // fiT = 9.f * w[i] * ( (cx[i]*Fx+cy[i]*Fy) * (cx[i]*ux+cy[i]*uy)*cs2 - (ux*Fx+uy*Fy)*cs2 );
        
    float _uNorm = sqrt(_ux*_ux+_uy*_uy);
    if(_uNorm>cs){
      _ux *= cs/_uNorm;
      _uy *= cs/_uNorm;
    }
    
    // Smagorinsky-Lilly subgrid turbulence model, source: https://arxiv.org/pdf/comp-gas/9401004.pdf : 0.09f = 8/(PI*PI*27*cs2) => sqrt(8) / (PI*PI*sqrt(27)*sqrt(cb(Ck))*cs2) with Ck=1.5 
    float tau = 0.5f + p_nu/cs2 + 0.09f * sqrt(_Hxx*_Hxx + _Hyy*_Hyy + 2.f*_Hxy*_Hxy) / ((cs2+2.f*p_nu)*max(1e-3f,_p));
    
    p = _p;
    ux = _ux + 0.5f*_fx;
    uy = _uy + 0.5f*_fy;
    Sxx = (tau-1.f)/(2.f*tau)*(_Sxx-_Syy+_uy*_uy+_fx*_ux-_fy*_uy) + (tau+1.f)/(2.f*tau)*_ux*_ux + _fx*_ux;
    Syy = (tau-1.f)/(2.f*tau)*(_Syy-_Sxx+_ux*_ux+_fy*_uy-_fx*_ux) + (tau+1.f)/(2.f*tau)*_uy*_uy + _fy*_uy;
    Sxy = (1.f-1.f/tau)*_Sxy + _ux*_uy/tau + (2.f*tau-1.f)/(2.f*tau)*(_fx*_uy+_fy*_ux);
  }
}
