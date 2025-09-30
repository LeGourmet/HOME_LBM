// should use 4 time finer discreate field for phi (2 per directions)
// should update hi at half time fi(t) => fi(t+1) ; hi(t) => hi(t+0.5) => hi(t+1)
// should update solid for solid boundary treatment
// 0.15 velocity eq == 0.001 force
// rho utility ??
// ddf shifting and optimisation for digit extinction : fieq, fi, rho usage, ...

public enum CELL_TYPE { SOLID, FLUID, EQUILIBRIUM };
public enum COLOR_TYPE { PRESSURE, VELOCITY, TYPE };

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
  private float _H = 0.f;
  private float _phi = 0.f;
  
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
  public float getTmpPhi() { return this._phi; }
  
  public int getColor(COLOR_TYPE p_colorType) {
    if (type == CELL_TYPE.SOLID) return color(0);
    
    int palette[];
    float val;
    
    if(p_colorType == COLOR_TYPE.PRESSURE){
      palette = new int[]{color(68,1,84), color(59,82,139), color(33,145,140), color(94,201,98), color(253,231,37)};
      val = constrain(p*0.5f, 0.f, 1.f);
    }else if(p_colorType == COLOR_TYPE.VELOCITY){
      palette = new int[]{color(70,70,219), color(0,255,91), color(0,128,0), color(255,255,0), color(255,96,0), color(107,0,0), color(223,77,77)};
      val = constrain(sqrt(ux*ux + uy*uy)/cs, 0.f, 1.f);
    }else{
      palette = new int[]{color(40,40,180), color(150,150,200), color(221,221,221), color(200,150,150), color(180,40,40)};
      val = (phi<0.5f) ? 0.f : 1.f;
    }
      
    float x = val*0.999f*(palette.length-1.f);
    int idx = (int)floor(x);
    return lerpColor(palette[idx],palette[idx+1],  x-(float)(idx));
  }
  
  // ------------------------------------------------------ SETTERS ------------------------------------------------------
  public void setType(CELL_TYPE p_type) { this.type = p_type; }
  public void setPressure(float p_pressure) { this.p = max(1e-3f,p_pressure); }
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
  
  private PVector computeForces(float p_p, float p_ux, float p_uy, float p_Sxx, float p_Syy, float p_Sxy, float p_phi, float p_Fbx, float p_Fby, float p_GphiX, float p_GphiY, float p_GphiSQ, float p_GrhoX, float p_GrhoY){    
    // body forces
    float fx = p_Fbx;
    float fy = p_Fby;
    
    // pressure force
    fx += -p_p * p_GrhoX;
    fy += -p_p * p_GrhoY;
    
    // viscosity force
    fx += (p_ux*p_ux-p_Sxx) * p_GrhoX + (p_ux*p_uy-p_Sxy) * p_GrhoY;
    fy += (p_uy*p_ux-p_Sxy) * p_GrhoX + (p_uy*p_uy-p_Syy) * p_GrhoY;
    
    // surface tension force
    float a=0.005f, b=0.0005f;
    float surfaceTensionForceTmp = (a+b) * (24.f/interfacial_thickness * p_phi * (1.f-p_phi) * (1.f-2.f*p_phi) - 3.f*interfacial_thickness/2.f * p_GphiSQ);
    fx += surfaceTensionForceTmp * p_GphiX;
    fy += surfaceTensionForceTmp * p_GphiY;
    
    return new PVector(fx,fy);
  }
  
  public void streaming(int p_idX, int p_idY, LBM p_simulation) {
    if(type==CELL_TYPE.SOLID || type==CELL_TYPE.EQUILIBRIUM) return;
    
    float nu = 1.f/( (1.f-_phi)/nu_air + _phi/nu_fluid );
    float rho = max(1e-3f,rho_air + _phi * (rho_fluid - rho_air));
    
    float tau = 0.5f + nu/cs2;
    
    float GphiX=0.f, GphiY=0.f, GphiSQ=0.f;
    for(int i=0; i<5 ;i++){
      int idNx = (p_idX+D2Q5_cx[i]+p_simulation.getNx())%p_simulation.getNx();
      int idNy = (p_idY+D2Q5_cy[i]+p_simulation.getNy())%p_simulation.getNy();
      GphiX  += D2Q5_w[i] * D2Q5_cx[i] * p_simulation.getCell(idNx,idNy).getPhi() / cs2;
      GphiY  += D2Q5_w[i] * D2Q5_cy[i] * p_simulation.getCell(idNx,idNy).getPhi() / cs2;
      GphiSQ += D2Q5_w[i] * (p_simulation.getCell(idNx,idNy).getPhi() - phi) * 2.f/cs2;
    }
    
    float GrhoX = (rho_fluid - rho_air) * GphiX;
    float GrhoY = (rho_fluid - rho_air) * GphiY;
    
    PVector f = computeForces(p, ux, uy, Sxx, Syy, Sxy, phi, p_simulation.getForceX(p_idX, p_idY), p_simulation.getForceY(p_idX, p_idY), GphiX, GphiY, GphiSQ, GrhoX, GrhoY);
    
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
      
      int idNx = (p_idX+D2Q9_cx[j]+p_simulation.getNx())%p_simulation.getNx();
      int idNy = (p_idY+D2Q9_cy[j]+p_simulation.getNy())%p_simulation.getNy();

      CELL_TYPE typeN = p_simulation.getCell(idNx,idNy).getType();
      float pN = p_simulation.getCell(idNx,idNy).getPressure();
      float uxN = p_simulation.getCell(idNx,idNy).getVelocityX();
      float uyN = p_simulation.getCell(idNx,idNy).getVelocityY();
      float SxxN = p_simulation.getCell(idNx,idNy).getSxx();
      float SyyN = p_simulation.getCell(idNx,idNy).getSyy();
      float SxyN = p_simulation.getCell(idNx,idNy).getSxy();

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
      
      // (Guo forcing, Krueger p.233f)
      float F = D2Q9_w[i] * (1.f-1.f/(2.f*tau)) * (((D2Q9_cx[i]-ux)/cs2 + ((D2Q9_cx[i]*ux+D2Q9_cy[i]*uy)*D2Q9_cx[i])/cs4)*f.x + ((D2Q9_cy[i]-uy)/cs2 + ((D2Q9_cx[i]*ux+D2Q9_cy[i]*uy)*D2Q9_cy[i])/cs4)*f.y) ;
      _ux += fi * D2Q9_cx[i] + 0.5f * F;
      _uy += fi * D2Q9_cy[i] + 0.5f * F;
      
      _Sxx += fi * (D2Q9_cx[i]*D2Q9_cx[i]-1.f/3.f);
      _Syy += fi * (D2Q9_cy[i]*D2Q9_cy[i]-1.f/3.f);
      _Sxy += fi * (D2Q9_cx[i]*D2Q9_cy[i]);
      
      _Hxx += NeqFi * (D2Q9_cx[i]*D2Q9_cx[i]-1.f/3.f);
      _Hyy += NeqFi * (D2Q9_cy[i]*D2Q9_cy[i]-1.f/3.f);
      _Hxy += NeqFi * (D2Q9_cx[i]*D2Q9_cy[i]);
    }
    
    _H = sqrt(_Hxx*_Hxx + _Hyy*_Hyy + 2.f*_Hxy*_Hxy);
    
    for(int i=0; i<5 ;i++){
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
    }
  }
  
  public void collision(int p_idX, int p_idY, LBM p_simulation) {    
    if(type==CELL_TYPE.SOLID || type==CELL_TYPE.EQUILIBRIUM) return;   
    
    float nu = 1.f/( (1.f-_phi)/nu_air + _phi/nu_fluid );
    float rho = max(1e-3f,rho_air + _phi * (rho_fluid - rho_air));
    
    float GphiX=0.f, GphiY=0.f, GphiSQ=0.f;
    for(int i=0; i<5 ;i++){
      int idNx = (p_idX+D2Q5_cx[i]+p_simulation.getNx())%p_simulation.getNx();
      int idNy = (p_idY+D2Q5_cy[i]+p_simulation.getNy())%p_simulation.getNy();
      GphiX += D2Q5_w[i] * D2Q5_cx[i] * p_simulation.getCell(idNx,idNy).getTmpPhi() / cs2;
      GphiY += D2Q5_w[i] * D2Q5_cy[i] * p_simulation.getCell(idNx,idNy).getTmpPhi() / cs2;
      GphiSQ += D2Q5_w[i] * (p_simulation.getCell(idNx,idNy).getTmpPhi() - _phi) * 2.f/cs2;
    }
    
    float GrhoX = (rho_fluid - rho_air) * GphiX;
    float GrhoY = (rho_fluid - rho_air) * GphiY;
    
    
    PVector f = computeForces(_p, _ux, _uy, _Sxx, _Syy, _Sxy, _phi, p_simulation.getForceX(p_idX, p_idY), p_simulation.getForceY(p_idX, p_idY), GphiX, GphiY, GphiSQ, GrhoX, GrhoY);
    float fx = f.x/_p; // should be devide by rho' in text ?? and should care of divide by _p=0!
    float fy = f.y/_p; // should be devide by rho' in text ?? and should care of divide by _p=0!
    
    float _uNorm = sqrt(sq(_ux+0.5f*fx)+sq(_uy+0.5f*fy));
    if(_uNorm>cs){
      _ux *= cs/_uNorm;
      _uy *= cs/_uNorm;
    }
    
    // Smagorinsky-Lilly subgrid turbulence model, source: https://arxiv.org/pdf/comp-gas/9401004.pdf : 0.09f = 8/(PI*PI*27*cs2) => sqrt(8) / (PI*PI*sqrt(27)*sqrt(cb(Ck))*cs2) with Ck=1.5 
    float tau = 0.5f + nu/cs2 + 0.09f * _H / ((cs2+2.f*nu)*_p);
    
    p = _p; // => leeman use +1 for extinction and ddf shifting    
    ux = _ux + 0.5f*fx; // (Guo forcing, Krueger p.233f)
    uy = _uy + 0.5f*fy; // (Guo forcing, Krueger p.233f)   
    Sxx = (tau-1.f)/(2.f*tau)*(_Sxx-_Syy+_uy*_uy+fx*_ux-fy*_uy) + (tau+1.f)/(2.f*tau)*_ux*_ux + fx*_ux;
    Syy = (tau-1.f)/(2.f*tau)*(_Syy-_Sxx+_ux*_ux+fy*_uy-fx*_ux) + (tau+1.f)/(2.f*tau)*_uy*_uy + fy*_uy;
    Sxy = (1.f-1.f/tau)*_Sxy + _ux*_uy/tau + (2.f*tau-1.f)/(2.f*tau)*(fx*_uy+fy*_ux);
    phi = _phi;
  }
}
