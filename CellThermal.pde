public class CellThermal {
  // ----------------------------------------------------- ATTRIBUTS -----------------------------------------------------
  private float T;
  private float[] gi;
  
  // --------------------------------------------- DESTRUCTOR / CONSTRUCTOR ----------------------------------------------  
  public CellThermal(float p_T, float p_ux, float p_uy) {
    this.T = p_T;
    
    this.gi = new float[5];
    for(int i=0; i<5 ;i++) this.gi[i] = computeGiEq(i,this.T,p_ux,p_uy);
  }

  // ------------------------------------------------------ GETTERS ------------------------------------------------------
  public float getT() { return this.T; }
  public float getGi(int p_i) { return this.gi[p_i]; }
  
  // ------------------------------------------------------ SETTERS ------------------------------------------------------
  public void setT(float p_T) { this.T = max(0.f,p_T); }
  
  // -------------------------------------------------- FUNCTIONS DDFs ---------------------------------------------------  
  private float computeGiEq(int p_i, float p_T, float p_ux, float p_uy) {
    return D2Q5_w[p_i] * p_T * (1.f + (D2Q5_cx[p_i]*p_ux + D2Q5_cy[p_i]*p_uy)/cs2);
  }
  
  // -------------------------------------------------- FUNCTIONS PHASE -------------------------------------------------- 
  public void collision(int p_xMacro, int p_yMacro, int p_xMicro, int p_yMicro, int p_scale, CELL_TYPE p_type, float p_ux, float p_uy, LBM p_simulation){ }
  
  public void streaming(int p_xMacro, int p_yMacro, int p_xMicro, int p_yMicro, int p_scale, CELL_TYPE p_type, float p_ux, float p_uy, LBM p_simulation){ }
  
}
