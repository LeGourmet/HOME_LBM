public class CellPhase {
  // ----------------------------------------------------- ATTRIBUTS -----------------------------------------------------
  private float phi;
  private float[] hi;
  
  // --------------------------------------------- DESTRUCTOR / CONSTRUCTOR ----------------------------------------------  
  public CellPhase(float p_phi, float p_ux, float p_uy) {
    this.phi = constrain(p_phi,0.f,1.f);
    
    this.hi = new float[5];
    for(int i=0; i<5 ;i++) this.hi[i] = computeHiEq(i,this.phi,p_ux,p_uy);
  }

  // ------------------------------------------------------ GETTERS ------------------------------------------------------
  public float getPhi() { return this.phi; }
  public float getHi(int p_i) { return this.hi[p_i]; }
  
  // ------------------------------------------------------ SETTERS ------------------------------------------------------
  public void setPhi(float p_phi) { this.phi = constrain(p_phi,0.f,1.f); }
  
  // -------------------------------------------------- FUNCTIONS DDFs ---------------------------------------------------  
  private float computeHiEq(int p_i, float p_phi, float p_ux, float p_uy) {
    return D2Q5_w[p_i] * p_phi * (1.f + (D2Q5_cx[p_i]*p_ux + D2Q5_cy[p_i]*p_uy)/cs2);
  }
  
  // -------------------------------------------------- FUNCTIONS PHASE -------------------------------------------------- 
  public void collision(int p_x, int p_y, int p_scale, CELL_TYPE p_type, float p_ux, float p_uy, LBM p_simulation){
    if(p_type==CELL_TYPE.SOLID || p_type==CELL_TYPE.EQUILIBRIUM) return;
    
    float mo = 1.f/( (1.f-phi)/mo_air + phi/mo_fluid );
    float tau = 0.5f + mo/cs2;
      
    if(false){
      // BGK => hi(x,t+1) = hi(x,t) -(fi-feq)/tau + Hi(x,t); 
       
      float GphiX = (D2Q5_w[1]*p_simulation.getMicroPhi(mod(p_x+1,p_simulation.getNx()*p_scale),p_y) - D2Q5_w[2]*p_simulation.getMicroPhi(mod(p_x-1,p_simulation.getNx()*p_scale),p_y))/cs2;
      float GphiY = (D2Q5_w[3]*p_simulation.getMicroPhi(p_x,mod(p_y+1,p_simulation.getNy()*p_scale)) - D2Q5_w[4]*p_simulation.getMicroPhi(p_x,mod(p_y-1,p_simulation.getNy()*p_scale)))/cs2;
      float nGphi = sqrt(sq(GphiX)+sq(GphiY));
      float nGphiX = (nGphi==0.f) ? 0.f : GphiX/nGphi;
      float nGphiY = (nGphi==0.f) ? 0.f : GphiY/nGphi;
      
      float ki = (phi==0.f) ? 0.f : 1.f;
      float k  = (phi==0.f || phi==1.f) ? 1.f : 2.f;
      
      hi[0] += (computeHiEq(0,phi,p_ux,p_uy)-hi[0])/tau;
      hi[1] += (computeHiEq(1,phi,p_ux,p_uy)-hi[1])/tau + D2Q5_w[1] * (1.f+ki/k) * 4.f/interfacial_thickness * phi*(1.f-phi) * missibility * nGphiX;
      hi[2] += (computeHiEq(2,phi,p_ux,p_uy)-hi[2])/tau - D2Q5_w[2] * (1.f+ki/k) * 4.f/interfacial_thickness * phi*(1.f-phi) * missibility * nGphiX;
      hi[3] += (computeHiEq(3,phi,p_ux,p_uy)-hi[3])/tau + D2Q5_w[3] * (1.f+ki/k) * 4.f/interfacial_thickness * phi*(1.f-phi) * missibility * nGphiY;
      hi[4] += (computeHiEq(4,phi,p_ux,p_uy)-hi[4])/tau - D2Q5_w[4] * (1.f+ki/k) * 4.f/interfacial_thickness * phi*(1.f-phi) * missibility * nGphiY;
    } else if(true)  {
      // MRT => hi(x,t+1) = hi(x,t) - M^(-1)*S*M*(hi(x,t) - hieq(x,t)) + M^(-1)*(1-0.5*S)*M*Hi(x,t);
      // MRT => hi(x,t+1) = hi(x,t) - M^(-1)*S*M*(hi(x,t) - hieq(x,t) - Hi(x,t)/2) + Hi(x,t);
      
      /*float[][] Mh = {
        {           1.f,               1.f,               1.f,               1.f,               1.f },
        {           -p_ux,            1.f-p_ux,           -1.f-p_ux,               -p_ux,               -p_ux },
        {           -p_uy,               -p_uy,               -p_uy,            1.f-p_uy,           -1.f-p_uy },
        { sq(p_ux)-sq(p_uy), sq(1.f-p_ux)-sq(p_uy), sq(1.f+p_ux)-sq(p_uy), sq(p_ux)-sq(1.f-p_uy), sq(p_ux)-sq(1.f+p_uy) },
        { sq(p_ux)+sq(p_uy), sq(1.f-p_ux)+sq(p_uy), sq(1.f+p_ux)+sq(p_uy), sq(p_ux)+sq(1.f-p_uy), sq(p_ux)+sq(1.f+p_uy) }
      };
            
      float[][] Sh = {
        { 1.f, 0.f, 0.f, 0.f, 0.f },
        { 0.f, sM1, 0.f, 0.f, 0.f },
        { 0.f, 0.f, sM1, 0.f, 0.f },
        { 0.f, 0.f, 0.f, 1.f, 0.f },
        { 0.f, 0.f, 0.f, 0.f, 1.f }
      };*/
      
      float sM1 = 1.f/(tau);
      float[][] M1 = {
        {1.f - 2.f * sq(p_ux) + 2.f * sM1 * sq(p_ux)                    - 2.f * sq(p_uy) + 2.f * sM1 * sq(p_uy), 
         1.f - sq(1.f - p_ux) - 2.f * sM1 * ( 1.f - p_ux) * p_ux - sq(p_ux) - 2.f * sq(p_uy) + 2.f * sM1 * sq(p_uy),
         1.f - sq(1.f - p_ux) - 2.f * sM1 * (-1.f - p_ux) * p_ux - sq(p_ux) - 2.f * sq(p_uy) + 2.f * sM1 * sq(p_uy),
         1.f - 2.f * sq(p_ux) + 2.f * sM1 * sq(p_ux) - sq(1.f - p_uy)     -       sq(p_uy) - 2.f * sM1 * (1.f - p_uy) * p_uy,
         1.f - 2.f * sq(p_ux) + 2.f * sM1 * sq(p_ux) - sq(1.f - p_uy)     -       sq(p_uy) - 2.f * sM1 * (-1.f - p_uy) * p_uy}, 
        {0.5f * p_ux * (1.f + p_ux) -         p_ux  * ( sM1 * (0.5f + p_ux)) + 0.25f * (sq(p_ux) - sq(p_uy))       + 0.25f * (sq(p_ux) + sq(p_uy)), 
         0.5f * p_ux * (1.f + p_ux) + ( 1.f - p_ux) * ( sM1 * (0.5f + p_ux)) + 0.25f * (sq(1.f - p_ux) - sq(p_uy)) + 0.25f * (sq(1.f - p_ux) + sq(p_uy)),
         0.5f * p_ux * (1.f + p_ux) + (-1.f - p_ux) * ( sM1 * (0.5f + p_ux)) + 0.25f * (sq(1.f - p_ux) - sq(p_uy)) + 0.25f * (sq(1.f - p_ux) + sq(p_uy)),
         0.5f * p_ux * (1.f + p_ux) -         p_ux  * ( sM1 * (0.5f + p_ux)) + 0.25f * (sq(p_ux) - sq(1.f - p_uy)) + 0.25f * (sq(p_ux) + sq(1.f - p_uy)),
         0.5f * p_ux * (1.f + p_ux) -         p_ux  * ( sM1 * (0.5f + p_ux)) + 0.25f * (sq(p_ux) - sq(1.f - p_uy)) + 0.25f * (sq(p_ux) + sq(1.f - p_uy))},
        {- ( sM1 * (-0.5f + p_ux)) *          p_ux + 0.5f * (-1.f + p_ux) * p_ux + 0.25f * (sq(p_ux) - sq(p_uy))       + 0.25f * (sq(p_ux) + sq(p_uy)),
           ( sM1 * (-0.5f + p_ux)) * ( 1.f - p_ux) + 0.5f * (-1.f + p_ux) * p_ux + 0.25f * (sq(1.f - p_ux) - sq(p_uy)) + 0.25f * (sq(1.f - p_ux) + sq(p_uy)),
           ( sM1 * (-0.5f + p_ux)) * (-1.f - p_ux) + 0.5f * (-1.f + p_ux) * p_ux + 0.25f * (sq(1.f - p_ux) - sq(p_uy)) + 0.25f * (sq(1.f - p_ux) + sq(p_uy)),
         - ( sM1 * (-0.5f + p_ux)) *          p_ux + 0.5f * (-1.f + p_ux) * p_ux + 0.25f * (sq(p_ux) - sq(1.f - p_uy)) + 0.25f * (sq(p_ux) + sq(1.f - p_uy)),
         - ( sM1 * (-0.5f + p_ux)) *          p_ux + 0.5f * (-1.f + p_ux) * p_ux + 0.25f * (sq(p_ux) - sq(1.f - p_uy)) + 0.25f * (sq(p_ux) + sq(1.f - p_uy))},
        {0.5f * p_uy * (1.f + p_uy) - 0.25f * (sq(p_ux) - sq(p_uy))       + 0.25f * (sq(p_ux) + sq(p_uy))       -         p_uy  * ( sM1 * (0.5f + p_uy)),
         0.5f * p_uy * (1.f + p_uy) - 0.25f * (sq(1.f - p_ux) - sq(p_uy)) + 0.25f * (sq(1.f - p_ux) + sq(p_uy)) -         p_uy  * ( sM1 * (0.5f + p_uy)),
         0.5f * p_uy * (1.f + p_uy) - 0.25f * (sq(1.f - p_ux) - sq(p_uy)) + 0.25f * (sq(1.f - p_ux) + sq(p_uy)) -         p_uy  * ( sM1 * (0.5f + p_uy)),
         0.5f * p_uy * (1.f + p_uy) - 0.25f * (sq(p_ux) - sq(1.f - p_uy)) + 0.25f * (sq(p_ux) + sq(1.f - p_uy)) +  (1.f - p_uy) * ( sM1 * (0.5f + p_uy)),
         0.5f * p_uy * (1.f + p_uy)  -0.25f * (sq(p_ux) - sq(1.f - p_uy)) + 0.25f * (sq(p_ux) + sq(1.f - p_uy)) + (-1.f - p_uy) * ( sM1 * (0.5f + p_uy))},
        {- ( sM1 * (-0.5f + p_uy)) * p_uy          + 0.5f * (-1.f + p_uy) * p_uy - 0.25f * (sq(p_ux) - sq(p_uy))       + 0.25f * (sq(p_ux) + sq(p_uy)),
         - ( sM1 * (-0.5f + p_uy)) * p_uy          + 0.5f * (-1.f + p_uy) * p_uy - 0.25f * (sq(1.f - p_ux) - sq(p_uy)) + 0.25f * (sq(1.f - p_ux) + sq(p_uy)),
         - ( sM1 * (-0.5f + p_uy)) * p_uy          + 0.5f * (-1.f + p_uy) * p_uy - 0.25f * (sq(1.f - p_ux) - sq(p_uy)) + 0.25f * (sq(1.f - p_ux) + sq(p_uy)),
           ( sM1 * (-0.5f + p_uy)) * ( 1.f - p_uy) + 0.5f * (-1.f + p_uy) * p_uy - 0.25f * (sq(p_ux) - sq(1.f - p_uy)) + 0.25f * (sq(p_ux) + sq(1.f - p_uy)),
           ( sM1 * (-0.5f + p_uy)) * (-1.f - p_uy) + 0.5f * (-1.f + p_uy) * p_uy - 0.25f * (sq(p_ux) - sq(1.f - p_uy)) + 0.25f * (sq(p_ux) + sq(1.f - p_uy))}
      };
      
  
      float GphiX = (D2Q5_w[1]*p_simulation.getMicroPhi(mod(p_x+1,p_simulation.getNx()*p_scale),p_y) - D2Q5_w[2]*p_simulation.getMicroPhi(mod(p_x-1,p_simulation.getNx()*p_scale),p_y))/cs2;
      float GphiY = (D2Q5_w[3]*p_simulation.getMicroPhi(p_x,mod(p_y+1,p_simulation.getNy()*p_scale)) - D2Q5_w[4]*p_simulation.getMicroPhi(p_x,mod(p_y-1,p_simulation.getNy()*p_scale)))/cs2;
      float nGphi = sqrt(sq(GphiX)+sq(GphiY));
      float nGphiX = (nGphi==0.f) ? 0.f : GphiX/nGphi;
      float nGphiY = (nGphi==0.f) ? 0.f : GphiY/nGphi;
      
      float[] H = {
         0.f,
         D2Q5_w[1] * 4.f/interfacial_thickness * phi*(1.f-phi) * missibility * nGphiX,
        -D2Q5_w[2] * 4.f/interfacial_thickness * phi*(1.f-phi) * missibility * nGphiX,
         D2Q5_w[3] * 4.f/interfacial_thickness * phi*(1.f-phi) * missibility * nGphiY,
        -D2Q5_w[4] * 4.f/interfacial_thickness * phi*(1.f-phi) * missibility * nGphiY 
      };
      
      float[] r1 = {0.f,0.f,0.f,0.f,0.f};
      for(int a=0; a<5 ;a++)
        for(int b=0; b<5 ;b++)
          r1[a] += M1[a][b]* (hi[b]-computeHiEq(b,phi,p_ux,p_uy)-H[b]*0.5f);
      
      // hi = hi - M1*(hi-hieq-H/2) + H
      hi[0] += - r1[0] + H[0];
      hi[1] += - r1[1] + H[1];
      hi[2] += - r1[2] + H[2];
      hi[3] += - r1[3] + H[3];
      hi[4] += - r1[4] + H[4];
    }
  }
  
  public void streaming(int p_x, int p_y, int p_scale, CELL_TYPE p_type, LBM p_simulation){
    if(p_type==CELL_TYPE.SOLID || p_type==CELL_TYPE.EQUILIBRIUM) return;
    
    phi = 0.f;
    
    for(int i=0; i<5 ;i++){
      int j = (i==0) ? i : ((i%2==0) ? i-1 : i+1);
      int x_micro = mod(p_x+D2Q5_cx[j],p_simulation.getNx()*p_scale);
      int y_micro = mod(p_y+D2Q5_cy[j],p_simulation.getNy()*p_scale);
      int x_macro = x_micro/p_scale;
      int y_macro = y_micro/p_scale;
      phi += (p_simulation.getType(x_macro,y_macro)==CELL_TYPE.SOLID) ? hi[j] : p_simulation.getHi(x_micro,y_micro,i);
    }
    
  }
  
}
