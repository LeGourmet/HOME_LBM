// should use 4 time finer discreate field for phi (2 per directions)
// should update hi at half time fi(t) => fi(t+1) ; hi(t) => hi(t+0.5) => hi(t+1)

public enum CELL_TYPE { SOLID, FLUID, AIR, EQUILIBRIUM };

public class Cell {
  // ----------------------------------------------------- ATTRIBUTS -----------------------------------------------------
  private CELL_TYPE type;
  
  private float Fx;
  private float Fy;
  
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
  private float _Hxx = 0.f;
  private float _Hyy = 0.f;
  private float _Hxy = 0.f;
  private float _phi = 0.f;
  
  // --------------------------------------------- DESTRUCTOR / CONSTRUCTOR ----------------------------------------------  
  public Cell(CELL_TYPE p_type, float p_Fx, float p_Fy, float p_p, float p_ux, float p_uy, float p_phi) {
    this.type = p_type;
    
    this.Fx = p_Fx;
    this.Fy = p_Fy;
    
    this.p = p_p;
    this.ux = p_ux;
    this.uy = p_uy;
    this.Sxx = (p_ux*p_ux - cs2);
    this.Syy = (p_uy*p_uy - cs2);
    this.Sxy = (p_ux*p_uy);  // Syx = Sxy; // symetric tensor
    this.phi = p_phi;
  } 
  
  // ------------------------------------------------------ GETTERS ------------------------------------------------------
  public CELL_TYPE getType() { return this.type; } 
  public float getForceX() { return this.Fx; }
  public float getForceY() { return this.Fy; }
  public float getPressure() { return this.p; }
  public float getVelocityX() { return this.ux; }
  public float getVelocityY() { return this.uy; }
  public float getSxx() { return this.Sxx; }
  public float getSyy() { return this.Syy; }
  public float getSxy() { return this.Sxy; }
  public float getPhi() { return this.phi; }
    
  public int getColor() {
    if (type == CELL_TYPE.SOLID) return color(0);
  
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
  
  // ------------------------------------------------------ SETTERS ------------------------------------------------------
  public void setType(CELL_TYPE p_type) { this.type = p_type; } 
  public void setForceX(float p_force ) { this.Fx = p_force; }
  public void setForceY(float p_force) { this.Fy = p_force; }
  public void setPressure(float p_pressure) { this.p = max(1e-3,p_pressure); }
  public void setVelocityX(float p_velocity) { this.ux = p_velocity; }
  public void setVelocityY(float p_velocity) { this.uy = p_velocity; }
  public void setPhi(float p_phi) { this.phi = constrain(p_phi,0.f,1.f); }
  
  // ----------------------------------------------------- FUNCTIONS -----------------------------------------------------
  private float computeFiEq(int p_i, float p_p, float p_ux, float p_uy) {
    return D2Q9_w[p_i] * (p_p + (D2Q9_cx[p_i]*p_ux + D2Q9_cy[p_i]*p_uy)/cs2 + 0.5f*sq(D2Q9_cx[p_i]*p_ux + D2Q9_cy[p_i]*p_uy)/cs4 - 0.5f*(p_ux*p_ux + p_uy*p_uy)/cs2);
  }
  
  private float computeFi(int p_i, float p_p, float p_ux, float p_uy, float p_Sxx, float p_Syy, float p_Sxy) {
    return D2Q9_w[p_i] * (p_p +
                          (D2Q9_cx[p_i]*p_ux + D2Q9_cy[p_i]*p_uy)/cs2 +
                          0.5f*( 2.f*p_Sxy*D2Q9_cx[p_i]*D2Q9_cy[p_i] + p_Sxx*(D2Q9_cx[p_i]*D2Q9_cx[p_i]-1.f/3.f) + p_Syy*(D2Q9_cy[p_i]*D2Q9_cy[p_i]-1.f/3.f))/cs4 +
                          0.5f*( (D2Q9_cx[p_i]*D2Q9_cx[p_i]*D2Q9_cy[p_i]-D2Q9_cy[p_i]*1.f/3.f) * (p_Sxx*p_uy+2.f*p_Sxy*p_ux-2.f*p_ux*p_ux*p_uy) +
                                 (D2Q9_cx[p_i]*D2Q9_cy[p_i]*D2Q9_cy[p_i]-D2Q9_cx[p_i]*1.f/3.f) * (p_Syy*p_ux+2.f*p_Sxy*p_uy-2.f*p_ux*p_uy*p_uy))/cs6);
  }
  
  private float computeHiEq(int p_i, float p_phi, float p_ux, float p_uy) {
    return D2Q5_w[p_i] * p_phi * (1.f + (D2Q5_cx[p_i]*p_ux + D2Q5_cy[p_i]*p_uy)/cs2);
  }
  
  private float computeHi(int p_i, float p_phi, float p_ux, float p_uy) {
    return computeFiEq(p_i, p_phi, p_ux, p_uy);
  }
  
  public void streaming(int p_idX, int p_idY, int p_Nx, int p_Ny, Cell[][] p_cells) {
    if(type==CELL_TYPE.SOLID || type==CELL_TYPE.EQUILIBRIUM) return;
    
    _p = 0.f;
    _ux = 0.f;
    _uy = 0.f;
    _Sxx = 0.f;
    _Syy = 0.f;
    _Sxy = 0.f;
    _Hxx = 0.f;
    _Hyy = 0.f;
    _Hxy = 0.f;
   
    for(int i=0; i<9 ;i++){
      int j = (i==0) ? i : ((i%2==0) ? i-1 : i+1);
      
      int idNx = (p_idX+D2Q9_cx[j]+p_Nx)%p_Nx;
      int idNy = (p_idY+D2Q9_cy[j]+p_Ny)%p_Ny;

      CELL_TYPE typeN = p_cells[idNx][idNy].getType();
      float pN = p_cells[idNx][idNy].getPressure();
      float uxN = p_cells[idNx][idNy].getVelocityX();
      float uyN = p_cells[idNx][idNy].getVelocityY();
      float SxxN = p_cells[idNx][idNy].getSxx();
      float SyyN = p_cells[idNx][idNy].getSyy();
      float SxyN = p_cells[idNx][idNy].getSxy();

      float _fi, _NeqFi;
      if(typeN==CELL_TYPE.SOLID){
        _fi = computeFi(i, p, uxN, uyN, Sxx+uxN*uxN-ux*ux, Syy+uyN*uyN-uy*uy, Sxy+uxN*uyN-ux*uy);
        _NeqFi = _fi - computeFiEq(i, p, uxN, uyN);
      } else if(typeN==CELL_TYPE.EQUILIBRIUM){
        _fi = computeFiEq(i, pN, uxN, uyN);
        _NeqFi = 0.f;
      } else { 
        _fi = computeFi(i, pN, uxN, uyN, SxxN, SyyN, SxyN);
        _NeqFi = _fi - computeFiEq(i, pN, uxN, uyN);
      }
      
      _p += _fi;
      _ux += _fi * D2Q9_cx[i];
      _uy += _fi * D2Q9_cy[i];
      _Sxx += _fi * (D2Q9_cx[i]*D2Q9_cx[i]-1.f/3.f);
      _Syy += _fi * (D2Q9_cy[i]*D2Q9_cy[i]-1.f/3.f);
      _Sxy += _fi * (D2Q9_cx[i]*D2Q9_cy[i]);
      
      _Hxx += _NeqFi * (D2Q9_cx[i]*D2Q9_cx[i]-1.f/3.f);
      _Hyy += _NeqFi * (D2Q9_cy[i]*D2Q9_cy[i]-1.f/3.f);
      _Hxy += _NeqFi * (D2Q9_cx[i]*D2Q9_cy[i]);
    }
    
    _phi = 0.f;
    //for(int i=0; i<5 ;i++)
      //hi[i] = computeHi(i, phi, ux, uy);
  }
  
  public void collision(float p_nu, float p_rho, float p_fx, float p_fy) {    
    if(type==CELL_TYPE.SOLID || type==CELL_TYPE.EQUILIBRIUM) return;
    
    float _fx=0.f, _fy=0.f;
    
    // body forces
    _fx += p_fx+Fx;
    _fy += p_fy+Fy;
    
    //float GphiX = 3 * D2Q5_w[i]*D2Q5_cx[i]*cellsphi[] ;
    //float GphiY = ;
    
    //float GrhoX = (p_rho1 - p_rho2) * GphiX;
    //float GrhoY = (p_rho1 - p_rho2) * GphiY;
    
    // pressure force
    //_fx += -_p * GrhoX;
    //_fy += -_p * GrhoY;
    
    // viscosity force
    _fx += 0.f;
    _fy += 0.f;
    
    // surface tension force
    _fx += 0.f;
    _fy += 0.f;
    
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
