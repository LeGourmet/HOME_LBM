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
  
  private float _p;
  private float _ux;
  private float _uy;
  private float _Sxx;
  private float _Syy;
  private float _Sxy;
  
  private float phi_1;
  private float phi_2;
  private float[] hi_1;
  private float[] hi_2; 
  
  // --------------------------------------------- DESTRUCTOR / CONSTRUCTOR ----------------------------------------------  
  public Cell(CELL_TYPE p_type, float p_p, float p_ux, float p_uy, float p_phi_1, float p_phi_2) {
    this.type = p_type;
    
    this.p = p_p;
    this.ux = p_ux;
    this.uy = p_uy;
    this.Sxx = (p_ux*p_ux - cs2);
    this.Syy = (p_uy*p_uy - cs2);
    this.Sxy = (p_ux*p_uy);  // Syx = Sxy; // symetric tensor
    
    
    this.phi_1 = constrain(p_phi_1,0.f,1.f);
    this.hi_1 = new float[9];
    for(int i=0; i<9 ;i++) this.hi_1[i] = computeHiEq(i, this.phi_1, this.ux, this.uy);
    
    this.phi_2 = constrain(p_phi_2,0.f,1.f);
    this.hi_2 = new float[9];
    for(int i=0; i<9 ;i++) this.hi_1[i] = computeHiEq(i, this.phi_2, this.ux, this.uy);
  } 
  
  // ------------------------------------------------------ GETTERS ------------------------------------------------------
  public CELL_TYPE getType() { return this.type; }
  public float getPressure() { return this.p; }
  public float getVelocityX() { return this.ux; }
  public float getVelocityY() { return this.uy; }
  public float getSxx() { return this.Sxx; }
  public float getSyy() { return this.Syy; }
  public float getSxy() { return this.Sxy; }
  public float getPhi1() { return this.phi_1; }
  public float getPhi2() { return this.phi_2; }
  public float getPhi3() { return 1.f-this.phi_1-this.phi_2; }
  
  // ------------------------------------------------------ SETTERS ------------------------------------------------------
  public void setType(CELL_TYPE p_type) { this.type = p_type; }
  public void setPressure(float p_pressure) { this.p = max(0.f,p_pressure); }
  public void setVelocityX(float p_velocity) { this.ux = p_velocity; }
  public void setVelocityY(float p_velocity) { this.uy = p_velocity; }
  
  // -------------------------------------------------- FUNCTIONS DDFs ---------------------------------------------------  
  private float computeFiEq(int p_i, float p_p, float p_ux, float p_uy) {
    return D2Q9_w[p_i] * (p_p + (D2Q9_cx[p_i]*p_ux + D2Q9_cy[p_i]*p_uy)/cs2 + 0.5f*sq(D2Q9_cx[p_i]*p_ux + D2Q9_cy[p_i]*p_uy)/cs4 - 0.5f*(p_ux*p_ux + p_uy*p_uy)/cs2 - 1.f);
  }
  
  private float computeFi(int p_i, float p_p, float p_ux, float p_uy, float p_Sxx, float p_Syy, float p_Sxy) {
    return D2Q9_w[p_i] * (p_p +
                          (D2Q9_cx[p_i]*p_ux + D2Q9_cy[p_i]*p_uy)/cs2 +
                          0.5f*( 2.f*p_Sxy*D2Q9_cx[p_i]*D2Q9_cy[p_i] + p_Sxx*(D2Q9_cx[p_i]*D2Q9_cx[p_i]-cs2) + p_Syy*(D2Q9_cy[p_i]*D2Q9_cy[p_i]-cs2))/cs4 +
                          0.5f*( (D2Q9_cx[p_i]*D2Q9_cx[p_i]*D2Q9_cy[p_i]-D2Q9_cy[p_i]*cs2) * (p_Sxx*p_uy+2.f*p_Sxy*p_ux-2.f*p_ux*p_ux*p_uy) +
                                 (D2Q9_cx[p_i]*D2Q9_cy[p_i]*D2Q9_cy[p_i]-D2Q9_cx[p_i]*cs2) * (p_Syy*p_ux+2.f*p_Sxy*p_uy-2.f*p_ux*p_uy*p_uy))/cs6 
                          - 1.f);
  }
  
  private float computeHiEq(int p_i, float p_phi, float p_ux, float p_uy) {
    return D2Q9_w[p_i] * p_phi * (1.f + (D2Q9_cx[p_i]*p_ux + D2Q9_cy[p_i]*p_uy)/cs2);
  }
  
  // -------------------------------------------------- FUNCTIONS FLOW ---------------------------------------------------  
  public void flowStreaming(int p_x, int p_y, LBM p_simulation) {
    CELL_TYPE type = p_simulation.getCell(p_x, p_y).getType();
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
      _Sxx += fi * (D2Q9_cx[i]*D2Q9_cx[i]-cs2);
      _Syy += fi * (D2Q9_cy[i]*D2Q9_cy[i]-cs2);
      _Sxy += fi * (D2Q9_cx[i]*D2Q9_cy[i]);
    }
  }
  
  public void flowCollision(int p_x, int p_y, LBM p_simulation) { 
    CELL_TYPE type = p_simulation.getCell(p_x, p_y).getType();
    if(type==CELL_TYPE.SOLID || type==CELL_TYPE.EQUILIBRIUM) return;
    
    float rho = max(1e-5f,phi_1*rho_1 + phi_2*rho_2 + (1.f-phi_1-phi_2)*rho_3);
    float nu = 1.f/(phi_1/nu_1 + phi_2/nu_2 + (1.f-phi_1-phi_2)/nu_3);
    float tau = 0.5f + nu/cs2;
    
    float Gphi1X = ( 
        D2Q9_w[1]*p_simulation.getCell(mod(p_x+1,p_simulation.getNx()),p_y).phi_1 
      - D2Q9_w[2]*p_simulation.getCell(mod(p_x-1,p_simulation.getNx()),p_y).phi_1
      + D2Q9_w[5]*p_simulation.getCell(mod(p_x+1,p_simulation.getNx()),p_y).phi_1 
      - D2Q9_w[6]*p_simulation.getCell(mod(p_x-1,p_simulation.getNx()),p_y).phi_1 
      + D2Q9_w[7]*p_simulation.getCell(mod(p_x+1,p_simulation.getNx()),p_y).phi_1 
      - D2Q9_w[8]*p_simulation.getCell(mod(p_x-1,p_simulation.getNx()),p_y).phi_1 
        )/cs2;
    
    float Gphi1Y = (
        D2Q9_w[3]*p_simulation.getCell(p_x,mod(p_y+1,p_simulation.getNy())).phi_1
      - D2Q9_w[4]*p_simulation.getCell(p_x,mod(p_y-1,p_simulation.getNy())).phi_1
      + D2Q9_w[5]*p_simulation.getCell(p_x,mod(p_y+1,p_simulation.getNy())).phi_1 
      - D2Q9_w[6]*p_simulation.getCell(p_x,mod(p_y-1,p_simulation.getNy())).phi_1 
      - D2Q9_w[7]*p_simulation.getCell(p_x,mod(p_y-1,p_simulation.getNy())).phi_1 
      + D2Q9_w[8]*p_simulation.getCell(p_x,mod(p_y+1,p_simulation.getNy())).phi_1   
        )/cs2;
    
    float Gphi2X = ( 
        D2Q9_w[1]*p_simulation.getCell(mod(p_x+1,p_simulation.getNx()),p_y).phi_2 
      - D2Q9_w[2]*p_simulation.getCell(mod(p_x-1,p_simulation.getNx()),p_y).phi_2
      + D2Q9_w[5]*p_simulation.getCell(mod(p_x+1,p_simulation.getNx()),p_y).phi_2 
      - D2Q9_w[6]*p_simulation.getCell(mod(p_x-1,p_simulation.getNx()),p_y).phi_2 
      + D2Q9_w[7]*p_simulation.getCell(mod(p_x+1,p_simulation.getNx()),p_y).phi_2 
      - D2Q9_w[8]*p_simulation.getCell(mod(p_x-1,p_simulation.getNx()),p_y).phi_2 
        )/cs2;
    
    float Gphi2Y = (
        D2Q9_w[3]*p_simulation.getCell(p_x,mod(p_y+1,p_simulation.getNy())).phi_2
      - D2Q9_w[4]*p_simulation.getCell(p_x,mod(p_y-1,p_simulation.getNy())).phi_2
      + D2Q9_w[5]*p_simulation.getCell(p_x,mod(p_y+1,p_simulation.getNy())).phi_2 
      - D2Q9_w[6]*p_simulation.getCell(p_x,mod(p_y-1,p_simulation.getNy())).phi_2 
      - D2Q9_w[7]*p_simulation.getCell(p_x,mod(p_y-1,p_simulation.getNy())).phi_2 
      + D2Q9_w[8]*p_simulation.getCell(p_x,mod(p_y+1,p_simulation.getNy())).phi_2   
        )/cs2;
    
    float GrhoX = (rho_1-rho_3) * Gphi1X + (rho_2-rho_3) * Gphi2X;
    float GrhoY = (rho_1-rho_3) * Gphi1Y + (rho_2-rho_3) * Gphi2Y;
    
    // unroll => 2*wi*(phi(x+cix,y+ciy)/cs2
    /*float GphiSQ= 2.f/cs * ( D2Q5_w[1]*(p_simulation.getCell(mod(p_x+1,p_simulation.getNx()),p_y).getPhi()-phi) + 
                             D2Q5_w[2]*(p_simulation.getCell(mod(p_x-1,p_simulation.getNx()),p_y).getPhi()-phi) + 
                             D2Q5_w[3]*(p_simulation.getCell(p_x,mod(p_y+1,p_simulation.getNy())).getPhi()-phi) + 
                             D2Q5_w[4]*(p_simulation.getCell(p_x,mod(p_y-1,p_simulation.getNy())).getPhi()-phi) );*/
    
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
    float fsx = 0.f;//(ca_fluid + ca_air) * GphiX * (24.f/interfacial_thickness * sq(missibility) * phi * (1.f-phi) * (1.f-2.f*phi) - 3.f*interfacial_thickness/2.f * GphiSQ);
    float fsy = 0.f;//(ca_fluid + ca_air) * GphiY * (24.f/interfacial_thickness * sq(missibility) * phi * (1.f-phi) * (1.f-2.f*phi) - 3.f*interfacial_thickness/2.f * GphiSQ);
        
    // external forces : Newton second law of motion
    float fx = fbx + (fpx+fvx+fsx)/rho; 
    float fy = fby + (fpy+fvy+fsy)/rho;
            
    // ajoute de la pression jusqu'au buffer overflow ! => trouver mieux
    float _uNorm = sqrt(sq(_ux+0.5f*fx)+sq(_uy+0.5f*fy));
    if(_uNorm>cs){
      _ux *= cs/_uNorm;
      _uy *= cs/_uNorm;
    }
    
    // compute moment at time t+1
    p = _p; 
    ux = _ux + 0.5f*fx; // (Guo forcing, Krueger p.233f)
    uy = _uy + 0.5f*fy; // (Guo forcing, Krueger p.233f)
    Sxx = (tau-1.f)/(2.f*tau)*(_Sxx-_Syy+_uy*_uy+fx*_ux-fy*_uy) + (tau+1.f)/(2.f*tau)*_ux*_ux + fx*_ux;
    Syy = (tau-1.f)/(2.f*tau)*(_Syy-_Sxx+_ux*_ux+fy*_uy-fx*_ux) + (tau+1.f)/(2.f*tau)*_uy*_uy + fy*_uy;
    Sxy = (1.f-1.f/tau)*_Sxy + _ux*_uy/tau + (2.f*tau-1.f)/(2.f*tau)*(fx*_uy+fy*_ux);
  } 
  
  // -------------------------------------------------- FUNCTIONS PHASE ---------------------------------------------------  
  public void phaseCollision(int p_x, int p_y, LBM p_simulation){
    if(type==CELL_TYPE.SOLID || type==CELL_TYPE.EQUILIBRIUM) return;
    
    float mo = 1.f/(phi_1/mo_1 + phi_2/mo_2 + (1.f-phi_1-phi_2)/mo_3);
    float tau = 0.5f + mo/cs2;
        
    float Gphi1X = ( 
        D2Q9_w[1]*p_simulation.getCell(mod(p_x+1,p_simulation.getNx()),p_y).phi_1 
      - D2Q9_w[2]*p_simulation.getCell(mod(p_x-1,p_simulation.getNx()),p_y).phi_1
      + D2Q9_w[5]*p_simulation.getCell(mod(p_x+1,p_simulation.getNx()),p_y).phi_1 
      - D2Q9_w[6]*p_simulation.getCell(mod(p_x-1,p_simulation.getNx()),p_y).phi_1 
      + D2Q9_w[7]*p_simulation.getCell(mod(p_x+1,p_simulation.getNx()),p_y).phi_1 
      - D2Q9_w[8]*p_simulation.getCell(mod(p_x-1,p_simulation.getNx()),p_y).phi_1 
        )/cs2;
    
    float Gphi1Y = (
        D2Q9_w[3]*p_simulation.getCell(p_x,mod(p_y+1,p_simulation.getNy())).phi_1
      - D2Q9_w[4]*p_simulation.getCell(p_x,mod(p_y-1,p_simulation.getNy())).phi_1
      + D2Q9_w[5]*p_simulation.getCell(p_x,mod(p_y+1,p_simulation.getNy())).phi_1 
      - D2Q9_w[6]*p_simulation.getCell(p_x,mod(p_y-1,p_simulation.getNy())).phi_1 
      - D2Q9_w[7]*p_simulation.getCell(p_x,mod(p_y-1,p_simulation.getNy())).phi_1 
      + D2Q9_w[8]*p_simulation.getCell(p_x,mod(p_y+1,p_simulation.getNy())).phi_1   
        )/cs2;
    
    float nGphi1 = sqrt(sq(Gphi1X)+sq(Gphi1Y));
    float nGphi1X = (nGphi1==0.f) ? 0.f : Gphi1X/nGphi1;
    float nGphi1Y = (nGphi1==0.f) ? 0.f : Gphi1Y/nGphi1;
    
    hi_1[0] += (computeHiEq(0,phi_1,ux,uy)-hi_1[0])/tau;
    hi_1[1] += (computeHiEq(1,phi_1,ux,uy)-hi_1[1])/tau +  nGphi1X            * D2Q9_w[1] * 4.f/interfacial_thickness * phi_1*(1.f-phi_1) * missibility;
    hi_1[2] += (computeHiEq(2,phi_1,ux,uy)-hi_1[2])/tau -  nGphi1X            * D2Q9_w[2] * 4.f/interfacial_thickness * phi_1*(1.f-phi_1) * missibility;
    hi_1[3] += (computeHiEq(3,phi_1,ux,uy)-hi_1[3])/tau +            nGphi1Y  * D2Q9_w[3] * 4.f/interfacial_thickness * phi_1*(1.f-phi_1) * missibility;
    hi_1[4] += (computeHiEq(4,phi_1,ux,uy)-hi_1[4])/tau -            nGphi1Y  * D2Q9_w[4] * 4.f/interfacial_thickness * phi_1*(1.f-phi_1) * missibility;
    hi_1[5] += (computeHiEq(5,phi_1,ux,uy)-hi_1[5])/tau + (nGphi1X + nGphi1Y) * D2Q9_w[5] * 4.f/interfacial_thickness * phi_1*(1.f-phi_1) * missibility;
    hi_1[6] += (computeHiEq(6,phi_1,ux,uy)-hi_1[6])/tau - (nGphi1X + nGphi1Y) * D2Q9_w[6] * 4.f/interfacial_thickness * phi_1*(1.f-phi_1) * missibility;
    hi_1[7] += (computeHiEq(7,phi_1,ux,uy)-hi_1[7])/tau + (nGphi1X - nGphi1Y) * D2Q9_w[7] * 4.f/interfacial_thickness * phi_1*(1.f-phi_1) * missibility;
    hi_1[8] += (computeHiEq(8,phi_1,ux,uy)-hi_1[8])/tau - (nGphi1X - nGphi1Y) * D2Q9_w[8] * 4.f/interfacial_thickness * phi_1*(1.f-phi_1) * missibility;
    
    float Gphi2X = ( 
        D2Q9_w[1]*p_simulation.getCell(mod(p_x+1,p_simulation.getNx()),p_y).phi_2 
      - D2Q9_w[2]*p_simulation.getCell(mod(p_x-1,p_simulation.getNx()),p_y).phi_2
      + D2Q9_w[5]*p_simulation.getCell(mod(p_x+1,p_simulation.getNx()),p_y).phi_2 
      - D2Q9_w[6]*p_simulation.getCell(mod(p_x-1,p_simulation.getNx()),p_y).phi_2 
      + D2Q9_w[7]*p_simulation.getCell(mod(p_x+1,p_simulation.getNx()),p_y).phi_2 
      - D2Q9_w[8]*p_simulation.getCell(mod(p_x-1,p_simulation.getNx()),p_y).phi_2 
        )/cs2;
    
    float Gphi2Y = (
        D2Q9_w[3]*p_simulation.getCell(p_x,mod(p_y+1,p_simulation.getNy())).phi_2
      - D2Q9_w[4]*p_simulation.getCell(p_x,mod(p_y-1,p_simulation.getNy())).phi_2
      + D2Q9_w[5]*p_simulation.getCell(p_x,mod(p_y+1,p_simulation.getNy())).phi_2 
      - D2Q9_w[6]*p_simulation.getCell(p_x,mod(p_y-1,p_simulation.getNy())).phi_2 
      - D2Q9_w[7]*p_simulation.getCell(p_x,mod(p_y-1,p_simulation.getNy())).phi_2 
      + D2Q9_w[8]*p_simulation.getCell(p_x,mod(p_y+1,p_simulation.getNy())).phi_2   
        )/cs2;
    
    float nGphi2 = sqrt(sq(Gphi2X)+sq(Gphi2Y));
    float nGphi2X = (nGphi2==0.f) ? 0.f : Gphi2X/nGphi2;
    float nGphi2Y = (nGphi2==0.f) ? 0.f : Gphi2Y/nGphi2;
    
    hi_2[0] += (computeHiEq(0,phi_2,ux,uy)-hi_2[0])/tau;
    hi_2[1] += (computeHiEq(1,phi_2,ux,uy)-hi_2[1])/tau +  nGphi2X            * D2Q9_w[1] * 4.f/interfacial_thickness * phi_2*(1.f-phi_2) * missibility;
    hi_2[2] += (computeHiEq(2,phi_2,ux,uy)-hi_2[2])/tau -  nGphi2X            * D2Q9_w[2] * 4.f/interfacial_thickness * phi_2*(1.f-phi_2) * missibility;
    hi_2[3] += (computeHiEq(3,phi_2,ux,uy)-hi_2[3])/tau +            nGphi2Y  * D2Q9_w[3] * 4.f/interfacial_thickness * phi_2*(1.f-phi_2) * missibility;
    hi_2[4] += (computeHiEq(4,phi_2,ux,uy)-hi_2[4])/tau -            nGphi2Y  * D2Q9_w[4] * 4.f/interfacial_thickness * phi_2*(1.f-phi_2) * missibility;
    hi_2[5] += (computeHiEq(5,phi_2,ux,uy)-hi_2[5])/tau + (nGphi2X + nGphi2Y) * D2Q9_w[5] * 4.f/interfacial_thickness * phi_2*(1.f-phi_2) * missibility;
    hi_2[6] += (computeHiEq(6,phi_2,ux,uy)-hi_2[6])/tau - (nGphi2X + nGphi2Y) * D2Q9_w[6] * 4.f/interfacial_thickness * phi_2*(1.f-phi_2) * missibility;
    hi_2[7] += (computeHiEq(7,phi_2,ux,uy)-hi_2[7])/tau + (nGphi2X - nGphi2Y) * D2Q9_w[7] * 4.f/interfacial_thickness * phi_2*(1.f-phi_2) * missibility;
    hi_2[8] += (computeHiEq(8,phi_2,ux,uy)-hi_2[8])/tau - (nGphi2X - nGphi2Y) * D2Q9_w[8] * 4.f/interfacial_thickness * phi_2*(1.f-phi_2) * missibility;
  }
  
  public void phaseStreaming(int p_x, int p_y, LBM p_simulation){
    if(type==CELL_TYPE.SOLID || type==CELL_TYPE.EQUILIBRIUM) return;
    
    phi_1 = 0.f;
    phi_2 = 0.f;
    
    for(int i=0; i<9 ;i++){
      int j = (i==0) ? i : ((i%2==0) ? i-1 : i+1);
      Cell cell = p_simulation.getCell(mod(p_x+D2Q9_cx[j],p_simulation.getNx()),mod(p_y+D2Q9_cy[j],p_simulation.getNy())); 
      phi_1 += (cell.type == CELL_TYPE.SOLID) ? hi_1[j] : cell.hi_1[i];
      phi_2 += (cell.type == CELL_TYPE.SOLID) ? hi_2[j] : cell.hi_2[i];
    }
  }
}
