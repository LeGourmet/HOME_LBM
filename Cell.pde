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
  
  private float _p = 0.f;
  private float _ux = 0.f;
  private float _uy = 0.f;
  private float _Sxx = 0.f;
  private float _Syy = 0.f;
  private float _Sxy = 0.f;
  
  private float phi;
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
    for(int i=0; i<5 ;i++) this.hi[i] = computeHiEq(i,this.phi,this.ux,this.uy);
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
                          0.5f*( 2.f*p_Sxy*D2Q9_cx[p_i]*D2Q9_cy[p_i] + p_Sxx*(D2Q9_cx[p_i]*D2Q9_cx[p_i]-cs2) + p_Syy*(D2Q9_cy[p_i]*D2Q9_cy[p_i]-cs2))/cs4 +
                          0.5f*( (D2Q9_cx[p_i]*D2Q9_cx[p_i]*D2Q9_cy[p_i]-D2Q9_cy[p_i]*cs2) * (p_Sxx*p_uy+2.f*p_Sxy*p_ux-2.f*p_ux*p_ux*p_uy) +
                                 (D2Q9_cx[p_i]*D2Q9_cy[p_i]*D2Q9_cy[p_i]-D2Q9_cx[p_i]*cs2) * (p_Syy*p_ux+2.f*p_Sxy*p_uy-2.f*p_ux*p_uy*p_uy))/cs6 
                          - 1.f);
  }
  
  private float computeHiEq(int p_i, float p_phi, float p_ux, float p_uy) {
    return D2Q5_w[p_i] * p_phi * (1.f + (D2Q5_cx[p_i]*p_ux + D2Q5_cy[p_i]*p_uy)/cs2);
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
      _Sxx += fi * (D2Q9_cx[i]*D2Q9_cx[i]-cs2);
      _Syy += fi * (D2Q9_cy[i]*D2Q9_cy[i]-cs2);
      _Sxy += fi * (D2Q9_cx[i]*D2Q9_cy[i]);
    }
  }
  
  public void flowCollision(int p_idX, int p_idY, LBM p_simulation) {    
    if(type==CELL_TYPE.SOLID || type==CELL_TYPE.EQUILIBRIUM) return;
    
    float rho = max(1e-5f,(1-phi)*rho_air + phi*rho_fluid);
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
    //if((nGphi*0.6f) > (4.f*phi*(1.f-phi)/interfacial_thickness))
    
    // body forces
    float fbx = p_simulation.getForceX(p_idX, p_idY);
    float fby = p_simulation.getForceY(p_idX, p_idY);
    
    // pressure force : -P * cs2 * gradient(rho) : -P * gradient(rho)
    float fpx = 0.f* -p * cs2 * GrhoX;
    float fpy = 0.f* -p * cs2 * GrhoY;
    
    // viscosity force : nu * [gradient(U) + transpose(grandient(u)] * gradient(rho)
    float fvx = ( (ux*ux-Sxx) * GrhoX + (ux*uy-Sxy) * GrhoY);
    float fvy = ( (uy*ux-Sxy) * GrhoX + (uy*uy-Syy) * GrhoY);
    
    // surface tension force : abs(ca_fluid - ca_air) < cs/100.f
    float fsx = 0.f* (ca_fluid + ca_air) * GphiX * (24.f/interfacial_thickness * (phi - 3.f*sq(phi) + 2.f*cb(phi)) + 3.f*interfacial_thickness/2.f * GphiSQ); // not missible => 1
    float fsy = 0.f* (ca_fluid + ca_air) * GphiY * (24.f/interfacial_thickness * (phi - 3.f*sq(phi) + 2.f*cb(phi)) + 3.f*interfacial_thickness/2.f * GphiSQ); // not missible => 1
    //fx += (ca_fluid - ca_air) * 4.f/interfacial_thickness * GphiSQ * GphiX; // fully missible => 0 
    //fy += (ca_fluid - ca_air) * 4.f/interfacial_thickness * GphiSQ * GphiY; // fully missible => 0 
        
    // external forces : Newton second law of motion
    float fx = fbx + (fpx+fvx+fsx)/rho; 
    float fy = fby + (fpy+fvy+fsy)/rho;
    
    // ajoute de la pression jusqu'au buffer overflow ! => trouver mieux
    float _uNorm = sqrt(sq(_ux+0.5f*fx)+sq(_uy+0.5f*fy));
    if(_uNorm>cs){
      _ux *= cs/_uNorm;
      _uy *= cs/_uNorm;
    }
    
    p = _p; 
    ux = _ux + 0.5f*fx; // (Guo forcing, Krueger p.233f)
    uy = _uy + 0.5f*fy; // (Guo forcing, Krueger p.233f)
    Sxx = (tau-1.f)/(2.f*tau)*(_Sxx-_Syy+_uy*_uy+fx*_ux-fy*_uy) + (tau+1.f)/(2.f*tau)*_ux*_ux + fx*_ux;
    Syy = (tau-1.f)/(2.f*tau)*(_Syy-_Sxx+_ux*_ux+fy*_uy-fx*_ux) + (tau+1.f)/(2.f*tau)*_uy*_uy + fy*_uy;
    Sxy = (1.f-1.f/tau)*_Sxy + _ux*_uy/tau + (2.f*tau-1.f)/(2.f*tau)*(fx*_uy+fy*_ux);
  }
  
  // -------------------------------------------------- FUNCTIONS PHASE -------------------------------------------------- 
  public void phaseCollision(int p_idX, int p_idY, LBM p_simulation){
    if(type==CELL_TYPE.SOLID || type==CELL_TYPE.EQUILIBRIUM) return;
    
    float mo = 1.f/( (1.f-phi)/mo_air + phi/mo_fluid );
    float tau = 0.5f + mo/cs2;
      
    if(false){
      // BGK => hi(x,t+1) = hi(x,t) -(fi-feq)/tau + Hi(x,t); 
      
      float GphiX = (D2Q5_w[1]*p_simulation.getCell(mod(p_idX+1,p_simulation.getNx()),p_idY).getPhi() - D2Q5_w[2]*p_simulation.getCell(mod(p_idX-1,p_simulation.getNx()),p_idY).getPhi())/cs2;
      float GphiY = (D2Q5_w[3]*p_simulation.getCell(p_idX,mod(p_idY+1,p_simulation.getNy())).getPhi() - D2Q5_w[4]*p_simulation.getCell(p_idX,mod(p_idY-1,p_simulation.getNy())).getPhi())/cs2;
      float nGphi = sqrt(sq(GphiX)+sq(GphiY));
      float nGphiX = (nGphi==0.f) ? 0.f : GphiX/nGphi;
      float nGphiY = (nGphi==0.f) ? 0.f : GphiY/nGphi;
      
      hi[0] += (computeHiEq(0,phi,ux,uy)-hi[0])/tau;
      hi[1] += (computeHiEq(1,phi,ux,uy)-hi[1])/tau + D2Q5_w[1] * 4.f/interfacial_thickness * phi*(1.f-phi) * nGphiX;
      hi[2] += (computeHiEq(2,phi,ux,uy)-hi[2])/tau - D2Q5_w[2] * 4.f/interfacial_thickness * phi*(1.f-phi) * nGphiX;
      hi[3] += (computeHiEq(3,phi,ux,uy)-hi[3])/tau + D2Q5_w[3] * 4.f/interfacial_thickness * phi*(1.f-phi) * nGphiY;
      hi[4] += (computeHiEq(4,phi,ux,uy)-hi[4])/tau - D2Q5_w[4] * 4.f/interfacial_thickness * phi*(1.f-phi) * nGphiY;
    } else {
      // MRT => hi(x,t+1) = hi(x,t) - M^(-1)*S*M*(hi(x,t)-hieq(x,t)) + M^(-1)*(1-0.5*S)*M*Hi(x,t);
      
      float sM1 = 1.f/(tau);
      float[][] M1 = {
        {1.f - 2.f * sq(ux) + 2.f * sM1 * sq(ux)                    - 2.f * sq(uy) + 2.f * sM1 * sq(uy), 
         1.f - sq(1.f - ux) - 2.f * sM1 * ( 1.f - ux) * ux - sq(ux) - 2.f * sq(uy) + 2.f * sM1 * sq(uy),
         1.f - sq(1.f - ux) - 2.f * sM1 * (-1.f - ux) * ux - sq(ux) - 2.f * sq(uy) + 2.f * sM1 * sq(uy),
         1.f - 2.f * sq(ux) + 2.f * sM1 * sq(ux) - sq(1.f - uy)     -       sq(uy) - 2.f * sM1 * (1.f - uy) * uy,
         1.f - 2.f * sq(ux) + 2.f * sM1 * sq(ux) - sq(1.f - uy)     -       sq(uy) - 2.f * sM1 * (-1.f - uy) * uy}, 
        {0.5f * ux * (1.f + ux) -         ux  * ( sM1 * (0.5f + ux)) + 0.25f * (sq(ux) - sq(uy))       + 0.25f * (sq(ux) + sq(uy)), 
         0.5f * ux * (1.f + ux) + ( 1.f - ux) * ( sM1 * (0.5f + ux)) + 0.25f * (sq(1.f - ux) - sq(uy)) + 0.25f * (sq(1.f - ux) + sq(uy)),
         0.5f * ux * (1.f + ux) + (-1.f - ux) * ( sM1 * (0.5f + ux)) + 0.25f * (sq(1.f - ux) - sq(uy)) + 0.25f * (sq(1.f - ux) + sq(uy)),
         0.5f * ux * (1.f + ux) -         ux  * ( sM1 * (0.5f + ux)) + 0.25f * (sq(ux) - sq(1.f - uy)) + 0.25f * (sq(ux) + sq(1.f - uy)),
         0.5f * ux * (1.f + ux) -         ux  * ( sM1 * (0.5f + ux)) + 0.25f * (sq(ux) - sq(1.f - uy)) + 0.25f * (sq(ux) + sq(1.f - uy))},
        {- ( sM1 * (-0.5f + ux)) *          ux + 0.5f * (-1.f + ux) * ux + 0.25f * (sq(ux) - sq(uy))       + 0.25f * (sq(ux) + sq(uy)),
           ( sM1 * (-0.5f + ux)) * ( 1.f - ux) + 0.5f * (-1.f + ux) * ux + 0.25f * (sq(1.f - ux) - sq(uy)) + 0.25f * (sq(1.f - ux) + sq(uy)),
           ( sM1 * (-0.5f + ux)) * (-1.f - ux) + 0.5f * (-1.f + ux) * ux + 0.25f * (sq(1.f - ux) - sq(uy)) + 0.25f * (sq(1.f - ux) + sq(uy)),
         - ( sM1 * (-0.5f + ux)) *          ux + 0.5f * (-1.f + ux) * ux + 0.25f * (sq(ux) - sq(1.f - uy)) + 0.25f * (sq(ux) + sq(1.f - uy)),
         - ( sM1 * (-0.5f + ux)) *          ux + 0.5f * (-1.f + ux) * ux + 0.25f * (sq(ux) - sq(1.f - uy)) + 0.25f * (sq(ux) + sq(1.f - uy))},
        {0.5f * uy * (1.f + uy) - 0.25f * (sq(ux) - sq(uy))       + 0.25f * (sq(ux) + sq(uy))       -         uy  * ( sM1 * (0.5f + uy)),
         0.5f * uy * (1.f + uy) - 0.25f * (sq(1.f - ux) - sq(uy)) + 0.25f * (sq(1.f - ux) + sq(uy)) -         uy  * ( sM1 * (0.5f + uy)),
         0.5f * uy * (1.f + uy) - 0.25f * (sq(1.f - ux) - sq(uy)) + 0.25f * (sq(1.f - ux) + sq(uy)) -         uy  * ( sM1 * (0.5f + uy)),
         0.5f * uy * (1.f + uy) - 0.25f * (sq(ux) - sq(1.f - uy)) + 0.25f * (sq(ux) + sq(1.f - uy)) +  (1.f - uy) * ( sM1 * (0.5f + uy)),
         0.5f * uy * (1.f + uy)  -0.25f * (sq(ux) - sq(1.f - uy)) + 0.25f * (sq(ux) + sq(1.f - uy)) + (-1.f - uy) * ( sM1 * (0.5f + uy))},
        {- ( sM1 * (-0.5f + uy)) * uy          + 0.5f * (-1.f + uy) * uy - 0.25f * (sq(ux) - sq(uy))       + 0.25f * (sq(ux) + sq(uy)),
         - ( sM1 * (-0.5f + uy)) * uy          + 0.5f * (-1.f + uy) * uy - 0.25f * (sq(1.f - ux) - sq(uy)) + 0.25f * (sq(1.f - ux) + sq(uy)),
         - ( sM1 * (-0.5f + uy)) * uy          + 0.5f * (-1.f + uy) * uy - 0.25f * (sq(1.f - ux) - sq(uy)) + 0.25f * (sq(1.f - ux) + sq(uy)),
           ( sM1 * (-0.5f + uy)) * ( 1.f - uy) + 0.5f * (-1.f + uy) * uy - 0.25f * (sq(ux) - sq(1.f - uy)) + 0.25f * (sq(ux) + sq(1.f - uy)),
           ( sM1 * (-0.5f + uy)) * (-1.f - uy) + 0.5f * (-1.f + uy) * uy - 0.25f * (sq(ux) - sq(1.f - uy)) + 0.25f * (sq(ux) + sq(1.f - uy))}
      };
      
    /*float[][] Mh = {
        {           1.f,               1.f,               1.f,               1.f,               1.f },
        {           -ux,            1.f-ux,           -1.f-ux,               -ux,               -ux },
        {           -uy,               -uy,               -uy,            1.f-uy,           -1.f-uy },
        { sq(ux)-sq(uy), sq(1.f-ux)-sq(uy), sq(1.f+ux)-sq(uy), sq(ux)-sq(1.f-uy), sq(ux)-sq(1.f+uy) },
        { sq(ux)+sq(uy), sq(1.f-ux)+sq(uy), sq(1.f+ux)+sq(uy), sq(ux)+sq(1.f-uy), sq(ux)+sq(1.f+uy) }
      };
            
      float[][] Sh = {
        { 1.f, 0.f, 0.f, 0.f, 0.f },
        { 0.f, sM1, 0.f, 0.f, 0.f },
        { 0.f, 0.f, sM1, 0.f, 0.f },
        { 0.f, 0.f, 0.f, 1.f, 0.f },
        { 0.f, 0.f, 0.f, 0.f, 1.f }
      };*/
        
      float[] hi_neq = {hi[0]-computeHiEq(0,phi,ux,uy), hi[1]-computeHiEq(1,phi,ux,uy), hi[2]-computeHiEq(2,phi,ux,uy), hi[3]-computeHiEq(3,phi,ux,uy), hi[4]-computeHiEq(4,phi,ux,uy)};
  
      // M1*(hi-hieq)
      float[] r1 = {0.f,0.f,0.f,0.f,0.f};
      for(int a=0; a<5 ;a++)
        for(int b=0; b<5 ;b++)
          r1[a] += M1[a][b]*hi_neq[b];
      
      float sM2 = 1.f-1.f/(2.f*tau);
      float[][] M2 = {
        {-1.f - 2.f * ux - 2.f * uy + 0.5f * (1.f - sq(ux) - sq(uy)) + (sq(ux) - sq(uy))       * (- 2.f * ux - sq(ux) - 2.f * uy - sq(uy)) -         ux  * (-2.f * sM2 * ux - sq(ux) - 2.f * uy - sq(uy)) -         uy  * (-2.f * ux - sq(ux) - 2.f * sM2 * uy - sq(uy)) + (0.5f - 2.f * ux - sq(ux) - 2.f * uy - sq(uy)) * (sq(ux) + sq(uy)),
         -1.f - 2.f * ux - 2.f * uy + 0.5f * (1.f - sq(ux) - sq(uy)) + (sq(1.f - ux) - sq(uy)) * (- 2.f * ux - sq(ux) - 2.f * uy - sq(uy)) + ( 1.f - ux) * (-2.f * sM2 * ux - sq(ux) - 2.f * uy - sq(uy)) -         uy  * (-2.f * ux - sq(ux) - 2.f * sM2 * uy - sq(uy)) + (0.5f - 2.f * ux - sq(ux) - 2.f * uy - sq(uy)) * (sq(1.f - ux) + sq(uy)),
         -1.f - 2.f * ux - 2.f * uy + 0.5f * (1.f - sq(ux) - sq(uy)) + (sq(1.f - ux) - sq(uy)) * (- 2.f * ux - sq(ux) - 2.f * uy - sq(uy)) + (-1.f - ux) * (-2.f * sM2 * ux - sq(ux) - 2.f * uy - sq(uy)) -         uy  * (-2.f * ux - sq(ux) - 2.f * sM2 * uy - sq(uy)) + (0.5f - 2.f * ux - sq(ux) - 2.f * uy - sq(uy)) * (sq(1.f - ux) + sq(uy)),
         -1.f - 2.f * ux - 2.f * uy + 0.5f * (1.f - sq(ux) - sq(uy)) + (sq(ux) - sq(1.f - uy)) * (- 2.f * ux - sq(ux) - 2.f * uy - sq(uy)) -         ux  * (-2.f * sM2 * ux - sq(ux) - 2.f * uy - sq(uy)) + ( 1.f - uy) * (-2.f * ux - sq(ux) - 2.f * sM2 * uy - sq(uy)) + (0.5f - 2.f * ux - sq(ux) - 2.f * uy - sq(uy)) * (sq(ux) + sq(1.f - uy)),
         -1.f - 2.f * ux - 2.f * uy + 0.5f * (1.f - sq(ux) - sq(uy)) + (sq(ux) - sq(1.f - uy)) * (- 2.f * ux - sq(ux) - 2.f * uy - sq(uy)) -         ux  * (-2.f * sM2 * ux - sq(ux) - 2.f * uy - sq(uy)) + (-1.f - uy) * (-2.f * ux - sq(ux) - 2.f * sM2 * uy - sq(uy)) + (0.5f - 2.f * ux - sq(ux) - 2.f * uy - sq(uy)) * (sq(ux) + sq(1.f - uy))},
        {1.f + ux + 0.25f * ux * (1.f + ux) -         ux  * (0.5f + sM2 * (0.5f + ux) + 0.5f * ux * (1.f + ux)) - (1.f + ux + 0.5f * ux * (1.f + ux)) *         uy  + (0.875f + ux + 0.5f * ux * (1.f + ux)) *       (sq(ux) - sq(uy)) + (0.875f + ux + 0.5f * ux * (1.f + ux)) * (sq(ux) + sq(uy)),
         1.f + ux + 0.25f * ux * (1.f + ux) + ( 1.f - ux) * (0.5f + sM2 * (0.5f + ux) + 0.5f * ux * (1.f + ux)) - (1.f + ux + 0.5f * ux * (1.f + ux)) *         uy  + (0.875f + ux + 0.5f * ux * (1.f + ux)) * (sq(1.f - ux) - sq(uy)) + (0.875f + ux + 0.5f * ux * (1.f + ux)) * (sq(1.f - ux) + sq(uy)),
         1.f + ux + 0.25f * ux * (1.f + ux) + (-1.f - ux) * (0.5f + sM2 * (0.5f + ux) + 0.5f * ux * (1.f + ux)) - (1.f + ux + 0.5f * ux * (1.f + ux)) *         uy  + (0.875f + ux + 0.5f * ux * (1.f + ux)) * (sq(1.f - ux) - sq(uy)) + (0.875f + ux + 0.5f * ux * (1.f + ux)) * (sq(1.f - ux) + sq(uy)),
         1.f + ux + 0.25f * ux * (1.f + ux) -         ux  * (0.5f + sM2 * (0.5f + ux) + 0.5f * ux * (1.f + ux)) + (1.f + ux + 0.5f * ux * (1.f + ux)) * ( 1.f - uy) + (0.875f + ux + 0.5f * ux * (1.f + ux)) * (sq(ux) + sq(1.f - uy)) + (0.875f + ux + 0.5f * ux * (1.f + ux)) * (sq(ux) - sq(1.f - uy)),
         1.f + ux + 0.25f * ux * (1.f + ux) -         ux  * (0.5f + sM2 * (0.5f + ux) + 0.5f * ux * (1.f + ux)) + (1.f + ux + 0.5f * ux * (1.f + ux)) * (-1.f - uy) + (0.875f + ux + 0.5f * ux * (1.f + ux)) * (sq(ux) + sq(1.f - uy)) + (0.875f + ux + 0.5f * ux * (1.f + ux)) * (sq(ux) - sq(1.f - uy))},
        {ux + 0.25f * (-1.f + ux) * ux -         ux  * (0.5f + sM2 * (-0.5f + ux) + 0.5f * (-1.f + ux) * ux) - ( ux + 0.5f * (-1.f + ux) * ux) *         uy  + (-0.125f + ux + 0.5f * (-1.f + ux) * ux) * (sq(ux) - sq(uy))       + (-0.125f + ux + 0.5f * (-1.f + ux) * ux) * (sq(ux) + sq(uy)),
         ux + 0.25f * (-1.f + ux) * ux + ( 1.f - ux) * (0.5f + sM2 * (-0.5f + ux) + 0.5f * (-1.f + ux) * ux) - ( ux + 0.5f * (-1.f + ux) * ux) *         uy  + (-0.125f + ux + 0.5f * (-1.f + ux) * ux) * (sq(1.f - ux) - sq(uy)) + (-0.125f + ux + 0.5f * (-1.f + ux) * ux) * (sq(1.f - ux) + sq(uy)),
         ux + 0.25f * (-1.f + ux) * ux + (-1.f - ux) * (0.5f + sM2 * (-0.5f + ux) + 0.5f * (-1.f + ux) * ux) - ( ux + 0.5f * (-1.f + ux) * ux) *         uy  + (-0.125f + ux + 0.5f * (-1.f + ux) * ux) * (sq(1.f - ux) - sq(uy)) + (-0.125f + ux + 0.5f * (-1.f + ux) * ux) * (sq(1.f - ux) + sq(uy)),
         ux + 0.25f * (-1.f + ux) * ux -         ux  * (0.5f + sM2 * (-0.5f + ux) + 0.5f * (-1.f + ux) * ux) + ( ux + 0.5f * (-1.f + ux) * ux) * ( 1.f - uy) + (-0.125f + ux + 0.5f * (-1.f + ux) * ux) * (sq(ux) - sq(1.f - uy)) + (-0.125f + ux + 0.5f * (-1.f + ux) * ux) * (sq(ux) + sq(1.f - uy)),
         ux + 0.25f * (-1.f + ux) * ux -         ux  * (0.5f + sM2 * (-0.5f + ux) + 0.5f * (-1.f + ux) * ux) + ( ux + 0.5f * (-1.f + ux) * ux) * (-1.f - uy) + (-0.125f + ux + 0.5f * (-1.f + ux) * ux) * (sq(ux) - sq(1.f - uy)) + (-0.125f + ux + 0.5f * (-1.f + ux) * ux) * (sq(ux) + sq(1.f - uy))},
        {0.5f + uy + 0.25f * uy * (1.f + uy) + (sq(ux) + sq(uy))       * (0.375f + uy + 0.5f * uy * (1.f + uy)) -         ux  * (0.5f + uy + 0.5f * uy * (1.f + uy)) + (sq(ux) - sq(uy))       * (0.625f + uy + 0.5f * uy * (1.f + uy)) -         uy  * ( sM2 * (0.5f + uy) + 0.5f * uy * (1.f + uy)),
         0.5f + uy + 0.25f * uy * (1.f + uy) + (sq(1.f - ux) + sq(uy)) * (0.375f + uy + 0.5f * uy * (1.f + uy)) + ( 1.f - ux) * (0.5f + uy + 0.5f * uy * (1.f + uy)) + (sq(1.f - ux) - sq(uy)) * (0.625f + uy + 0.5f * uy * (1.f + uy)) -         uy  * ( sM2 * (0.5f + uy) + 0.5f * uy * (1.f + uy)),
         0.5f + uy + 0.25f * uy * (1.f + uy) + (sq(1.f - ux) + sq(uy)) * (0.375f + uy + 0.5f * uy * (1.f + uy)) + (-1.f - ux) * (0.5f + uy + 0.5f * uy * (1.f + uy)) + (sq(1.f - ux) - sq(uy)) * (0.625f + uy + 0.5f * uy * (1.f + uy)) -         uy  * ( sM2 * (0.5f + uy) + 0.5f * uy * (1.f + uy)), 
         0.5f + uy + 0.25f * uy * (1.f + uy) + (sq(ux) + sq(1.f - uy)) * (0.375f + uy + 0.5f * uy * (1.f + uy)) -         ux  * (0.5f + uy + 0.5f * uy * (1.f + uy)) + (sq(ux) - sq(1.f - uy)) * (0.625f + uy + 0.5f * uy * (1.f + uy)) + ( 1.f - uy) * ( sM2 * (0.5f + uy) + 0.5f * uy * (1.f + uy)),
         0.5f + uy + 0.25f * uy * (1.f + uy) + (sq(ux) + sq(1.f - uy)) * (0.375f + uy + 0.5f * uy * (1.f + uy)) -         ux  * (0.5f + uy + 0.5f * uy * (1.f + uy)) + (sq(ux) - sq(1.f - uy)) * (0.625f + uy + 0.5f * uy * (1.f + uy)) + (-1.f - uy) * ( sM2 * (0.5f + uy) + 0.5f * uy * (1.f + uy))},
        {-0.5f + uy + 0.25f * (-1.f + uy) * uy -         uy  * ( sM2 * (-0.5f + uy) + 0.5f * (-1.f + uy) * uy) -         ux  * (-0.5f + uy + 0.5f * (-1.f + uy) * uy) + (-0.375f + uy + 0.5f * (-1.f + uy) * uy) * (sq(ux) - sq(uy))       + (-0.625f + uy + 0.5f * (-1.f + uy) * uy) * (sq(ux) + sq(uy)),
         -0.5f + uy + 0.25f * (-1.f + uy) * uy -         uy  * ( sM2 * (-0.5f + uy) + 0.5f * (-1.f + uy) * uy) + ( 1.f - ux) * (-0.5f + uy + 0.5f * (-1.f + uy) * uy) + (-0.375f + uy + 0.5f * (-1.f + uy) * uy) * (sq(1.f - ux) - sq(uy)) + (-0.625f + uy + 0.5f * (-1.f + uy) * uy) * (sq(1.f - ux) + sq(uy)),
         -0.5f + uy + 0.25f * (-1.f + uy) * uy -         uy  * ( sM2 * (-0.5f + uy) + 0.5f * (-1.f + uy) * uy) + (-1.f - ux) * (-0.5f + uy + 0.5f * (-1.f + uy) * uy) + (-0.375f + uy + 0.5f * (-1.f + uy) * uy) * (sq(1.f - ux) - sq(uy)) + (-0.625f + uy + 0.5f * (-1.f + uy) * uy) * (sq(1.f - ux) + sq(uy)),
         -0.5f + uy + 0.25f * (-1.f + uy) * uy + ( 1.f - uy) * ( sM2 * (-0.5f + uy) + 0.5f * (-1.f + uy) * uy) -         ux  * (-0.5f + uy + 0.5f * (-1.f + uy) * uy) + (-0.375f + uy + 0.5f * (-1.f + uy) * uy) * (sq(ux) - sq(1.f - uy)) + (-0.625f + uy + 0.5f * (-1.f + uy) * uy) * (sq(ux) + sq(1.f - uy)),
         -0.5f + uy + 0.25f * (-1.f + uy) * uy + (-1.f - uy) * ( sM2 * (-0.5f + uy) + 0.5f * (-1.f + uy) * uy) -          ux * (-0.5f + uy + 0.5f * (-1.f + uy) * uy) + (-0.375f + uy + 0.5f * (-1.f + uy) * uy) * (sq(ux) - sq(1.f - uy)) + (-0.625f + uy + 0.5f * (-1.f + uy) * uy) * (sq(ux) + sq(1.f - uy))}
      };
      
      float GphiX = (D2Q5_w[1]*p_simulation.getCell(mod(p_idX+1,p_simulation.getNx()),p_idY).getPhi() - D2Q5_w[2]*p_simulation.getCell(mod(p_idX-1,p_simulation.getNx()),p_idY).getPhi())/cs2;
      float GphiY = (D2Q5_w[3]*p_simulation.getCell(p_idX,mod(p_idY+1,p_simulation.getNy())).getPhi() - D2Q5_w[4]*p_simulation.getCell(p_idX,mod(p_idY-1,p_simulation.getNy())).getPhi())/cs2;
      float nGphi = sqrt(sq(GphiX)+sq(GphiY));
      float nGphiX = (nGphi==0.f) ? 0.f : GphiX/nGphi;
      float nGphiY = (nGphi==0.f) ? 0.f : GphiY/nGphi;
      
      float[] H = {
         0.f,
         D2Q5_w[1] * 4.f/interfacial_thickness * phi*(1.f-phi) * nGphiX,
        -D2Q5_w[2] * 4.f/interfacial_thickness * phi*(1.f-phi) * nGphiX,
         D2Q5_w[3] * 4.f/interfacial_thickness * phi*(1.f-phi) * nGphiY,
        -D2Q5_w[4] * 4.f/interfacial_thickness * phi*(1.f-phi) * nGphiY 
      };
      
    /*float[][] _Sh = {
        { 0.5f, 1.f, 1.f,  1.f,  1.f },
        {  1.f, sM2, 1.f,  1.f,  1.f },
        {  1.f, 1.f, sM2,  1.f,  1.f },
        {  1.f, 1.f, 1.f, 0.5f,  1.f },
        {  1.f, 1.f, 1.f,  1.f, 0.5f }
      };
        
      float[][] Mh_inv = {
        { -ux*ux-uy*uy+1.f, -2.f*ux, -2.f*uy,    0.f,  -1.f },
        { 0.5f*ux*(ux+1.f), ux+0.5f,     0.f,  0.25f, 0.25f },
        { 0.5f*ux*(ux-1.f), ux-0.5f,     0.f,  0.25f, 0.25f },
        { 0.5f*uy*(uy+1.f),     0.f, uy+0.5f, -0.25f, 0.25f },
        { 0.5f*uy*(uy-1.f),     0.f, uy-0.5f, -0.25f, 0.25f }
      };*/
            
      // M2*H
      float[] r2 = {0.f,0.f,0.f,0.f,0.f};
      for(int a=0; a<5 ;a++)
        for(int b=0; b<5 ;b++)
          r2[a] += M2[a][b]*H[b];
      
      // hi - M1*(hi-hieq) (fully missible) or hi - M1*(hi-hieq) + M2*(H) (not missible)
      hi[0] = hi[0] - r1[0] + H[0];// + r2[0]; 
      hi[1] = hi[1] - r1[1] + H[1];// + r2[1];
      hi[2] = hi[2] - r1[2] + H[2];// + r2[2];
      hi[3] = hi[3] - r1[3] + H[3];// + r2[3];
      hi[4] = hi[4] - r1[4] + H[4];// + r2[4];
    }
  }
  
  public void phaseStreaming(int p_idX, int p_idY, LBM p_simulation){
    if(type==CELL_TYPE.SOLID || type==CELL_TYPE.EQUILIBRIUM) return;
    
    float _phi = 0.f;
    
    for(int i=0; i<5 ;i++){
      int j = (i==0) ? i : ((i%2==0) ? i-1 : i+1);
      
      int idNx = mod(p_idX+D2Q5_cx[j],p_simulation.getNx());
      int idNy = mod(p_idY+D2Q5_cy[j],p_simulation.getNy());
      
      _phi += (p_simulation.getCell(idNx,idNy).getType()==CELL_TYPE.SOLID) ? hi[j] : p_simulation.getCell(idNx,idNy).hi[i]; 
    }
    
    //phi = constrain((phi+_phi)/2.f,0.f,1.f);
    phi = constrain(_phi,0.f,1.f);  
  }
}
