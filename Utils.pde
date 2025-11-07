// --------------------------------------------------- FLUIDS CONST ----------------------------------------------------
final float rho_1 = 1.f;
final float nu_1 = 0.004f;
final float mo_1 = 0.2f;
final float ca_1 = 0.26f;

final float rho_2 = 2.f;
final float nu_2 = 0.004f;
final float mo_2 = 0.2f;
final float ca_2 = 0.26f;

final float rho_3 = 3.f;
final float nu_3 = 0.004f;
final float mo_3 = 0.2f;
final float ca_3 = 0.26f;

final float missibility = 0.5f;               // 0=>fully missible ; ]0,1[=>partialy missible ; 1=>immissible
final float interfacial_thickness = 5.f;

// ----------------------------------------------------- LBM CONST -----------------------------------------------------
final float cs = 0.57735027f;
final float cs2 = 0.33333333f;
final float cs4 = 0.11111111f;
final float cs6 = 0.03703704f;
  
final float[] D2Q9_w = { 4.f/9.f, 1.f/9.f, 1.f/9.f, 1.f/9.f, 1.f/9.f, 1.f/36.f, 1.f/36.f, 1.f/36.f, 1.f/36.f };
final int[] D2Q9_cx = { 0, 1,-1, 0, 0, 1,-1, 1,-1 }; 
final int[] D2Q9_cy = { 0, 0, 0, 1,-1, 1,-1,-1, 1 };

// ------------------------------------------------------- UTILS -------------------------------------------------------
int mod(int x, int n) { return (x+n)%n; }
float cb(float x) { return x*x*x; }

boolean inSphere(int x, int y, int r, int cx, int cy){
  return ((cx-x)*(cx-x)+(cy-y)*(cy-y)) < (r*r);
}
