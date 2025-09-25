// should use 4 time finer discreate field for phi (2 per directions)
// should update hi at half time fi(t) => fi(t+1) ; hi(t) => hi(t+0.5) => hi(t+1)
// should update solid for solid boundary treatment
// 0.15 velocity eq == 0.001 force
// rho utility ??

public enum CELL_TYPE { SOLID, FLUID, EQUILIBRIUM };

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
  private float _NeqS_norm = 0.f;
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
    
  /*public int getColorPressure() {
    float val = constrain(p/3.f, 0.f, 1.f);
  }*/
    
  public int getColorType() {
    if (type == CELL_TYPE.SOLID) return color(127);
    
    int palette[] = { 
      color(0, 0, 255),
      color(255, 0, 0)
    };
  
    float val = constrain(phi, 0.f, 1.f);
  
    float x = constrain(val,0.00001f,0.99999f)*(palette.length-1.f);
    int idx = (int)floor(x);
    return lerpColor(palette[idx],palette[idx+1],  x-(float)(idx));
  }
  
  public int getColorVelocity() {
    if (type == CELL_TYPE.SOLID) return color(0);
      
    int palette[] = { 
      color(70, 70, 219),
      color(0, 255, 91), 
      color(0, 128, 0),
      color(255, 255, 0),
      color(255, 96, 0),
      color(107, 0, 0),
      color(223, 77, 77)
    };
    
    float val = constrain(sqrt(ux*ux + uy*uy)*2.5f, 0.f, 1.f);    
    
    float x = constrain(val,0.00001f,0.99999f)*(palette.length-1.f);
    int idx = (int)floor(x);
    return lerpColor(palette[idx],palette[idx+1],  x-(float)(idx));
  }
  
  public int getColor(){
    if (type == CELL_TYPE.SOLID) return color(127);
    
    int palette[] = { 
      color(10, 10, 10),
      color(60, 60, 200)
    };
  
    float val = constrain(phi, 0.f, 1.f);
  
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
    return computeHiEq(p_i, p_phi, p_ux, p_uy);
  }
  
  public void streaming(int p_idX, int p_idY, int p_Nx, int p_Ny, Cell[][] p_cells) {
    if(type==CELL_TYPE.SOLID || type==CELL_TYPE.EQUILIBRIUM) return;
    
    _p = 0.f;
    _ux = 0.f;
    _uy = 0.f;
    _Sxx = 0.f;
    _Syy = 0.f;
    _Sxy = 0.f;
    _phi = 0.f;
        
    float _Hxx = 0.f;
    float _Hyy = 0.f;
    float _Hxy = 0.f;
   
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

      float fi, NeqFi;
      if(typeN==CELL_TYPE.SOLID){
        fi = computeFi(i, p, uxN, uyN, Sxx+uxN*uxN-ux*ux, Syy+uyN*uyN-uy*uy, Sxy+uxN*uyN-ux*uy);
        NeqFi = fi - computeFiEq(i, p, uxN, uyN);
      } else if(typeN==CELL_TYPE.EQUILIBRIUM){
        fi = computeFiEq(i, pN, uxN, uyN);
        NeqFi = 0.f;
      } else { 
        fi = computeFi(i, pN, uxN, uyN, SxxN, SyyN, SxyN);
        NeqFi = fi - computeFiEq(i, pN, uxN, uyN);
      }
      
      _p += fi;
      _ux += fi * D2Q9_cx[i];
      _uy += fi * D2Q9_cy[i];
      _Sxx += fi * (D2Q9_cx[i]*D2Q9_cx[i]-1.f/3.f);
      _Syy += fi * (D2Q9_cy[i]*D2Q9_cy[i]-1.f/3.f);
      _Sxy += fi * (D2Q9_cx[i]*D2Q9_cy[i]);
      
      _Hxx += NeqFi * (D2Q9_cx[i]*D2Q9_cx[i]-1.f/3.f);
      _Hyy += NeqFi * (D2Q9_cy[i]*D2Q9_cy[i]-1.f/3.f);
      _Hxy += NeqFi * (D2Q9_cx[i]*D2Q9_cy[i]);
    }
    
    _NeqS_norm = sqrt(_Hxx*_Hxx + _Hyy*_Hyy + 2.f*_Hxy*_Hxy);
    
    for(int i=0; i<5 ;i++){
      int j = (i==0) ? i : ((i%2==0) ? i-1 : i+1);
      
      int idNx = (p_idX+D2Q5_cx[j]+p_Nx)%p_Nx;
      int idNy = (p_idY+D2Q5_cy[j]+p_Ny)%p_Ny;
      
      CELL_TYPE typeN = p_cells[idNx][idNy].getType();
      float phiN = p_cells[idNx][idNy].getPhi();
      float uxN = p_cells[idNx][idNy].getVelocityX();
      float uyN = p_cells[idNx][idNy].getVelocityY();
      
      if(typeN==CELL_TYPE.SOLID){
        _phi +=  computeHi(i, phi, uxN, uyN);
      } else if(typeN==CELL_TYPE.EQUILIBRIUM){
        _phi += computeHiEq(i, phiN, uxN, uyN);
      } else { 
        _phi += computeHi(i, phiN, uxN, uyN);
      }
    }
  }
  
  public void collision(int p_idX, int p_idY, int p_Nx, int p_Ny, Cell[][] p_cells, float p_fx, float p_fy) {    
    if(type==CELL_TYPE.SOLID || type==CELL_TYPE.EQUILIBRIUM) return;   
    
    float nu = 1.f/( (1.f-phi)/nu_air + phi/nu_fluid );
    float rho = max(1e-3f,rho_air + phi * (rho_fluid - rho_air));
    
    float GphiX=0.f, GphiY=0.f, GphiSQ=0.f;
    for(int i=0; i<5 ;i++){
      int idNx = (p_idX+D2Q5_cx[i]+p_Nx)%p_Nx;
      int idNy = (p_idY+D2Q5_cy[i]+p_Ny)%p_Ny;
      GphiX += D2Q5_w[i] * D2Q5_cx[i] * p_cells[idNx][idNy].getPhi() / cs2;
      GphiY += D2Q5_w[i] * D2Q5_cy[i] * p_cells[idNx][idNy].getPhi() / cs2;
      GphiSQ += D2Q5_w[i] * (p_cells[idNx][idNy].getPhi() - phi) * 2.f/cs2;
    }
    
    float GrhoX = (rho_fluid - rho_air) * GphiX;
    float GrhoY = (rho_fluid - rho_air) * GphiY;
    
    float _fx=0.f, _fy=0.f;
    
    // body forces
    _fx += p_fx+Fx;
    _fy += p_fy+Fy;
    
    // pressure force
    _fx += -p * GrhoX;
    _fy += -p * GrhoY;
    
    // viscosity force
    _fx += (_ux*_ux-_Sxx) * GrhoX + (_ux*_uy-_Sxy) * GrhoY;
    _fy += (_uy*_ux-_Sxy) * GrhoX + (_uy*_uy-_Syy) * GrhoY;
    
    // surface tension force
    float a=0.005f, b=0.0005f;
    float surfaceTensionForceTmp = (a+b) * (24.f/interfacial_thickness * phi * (1.f-phi) * (1.f-2.f*phi) - 3.f*interfacial_thickness/2.f * GphiSQ);
    _fx += surfaceTensionForceTmp * GphiX;
    _fy += surfaceTensionForceTmp * GphiY;
    
    // (Guo forcing, Krueger p.233f) for volume force
    _p = max(1e-3,_p);
    _fx /= _p;
    _fy /= _p;
    _ux += 0.5f*_fx;
    _uy += 0.5f*_fy;
    
    float _u_norm = sqrt(_ux*_ux+_uy*_uy);
    if(_u_norm>cs){
      _ux *= cs/_u_norm;
      _uy *= cs/_u_norm;
    }
    
    // Smagorinsky-Lilly subgrid turbulence model, source: https://arxiv.org/pdf/comp-gas/9401004.pdf : 0.09f = 8/(PI*PI*27*cs2) => sqrt(8) / (PI*PI*sqrt(27)*sqrt(cb(Ck))*cs2) with Ck=1.5 
    float tau = 0.5f + nu/cs2 + 0.09f * _NeqS_norm / ((cs2+2.f*nu)*_p);
    
    p = _p;
    ux = _ux;
    uy = _uy;
    Sxx = (tau-1.f)/(2.f*tau)*(_Sxx-_Syy+_uy*_uy+_fx*_ux-_fy*_uy) + (tau+1.f)/(2.f*tau)*_ux*_ux + _fx*_ux;
    Syy = (tau-1.f)/(2.f*tau)*(_Syy-_Sxx+_ux*_ux+_fy*_uy-_fx*_ux) + (tau+1.f)/(2.f*tau)*_uy*_uy + _fy*_uy;
    Sxy = (1.f-1.f/tau)*_Sxy + _ux*_uy/tau + (2.f*tau-1.f)/(2.f*tau)*(_fx*_uy+_fy*_ux);
    phi = _phi;
  }
}
