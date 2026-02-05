// --------------------------------------------------- FLUIDS CONST ----------------------------------------------------
final float rho_fluid = 3.f;                 // ]0, +inf[
final float nu_fluid = 0.004f;               // ]0, +inf[
final float mo_fluid = 0.2f;                 // ]0, +inf[
final float ca_fluid = 0.26f;                // [0, +inf[

final float rho_air = 1.f;                   // ]0, +inf[
final float nu_air = 0.004f;                 // ]0, +inf[
final float mo_air = 0.2f;                   // ]0, +inf[
final float ca_air = 0.26f;                  // [0, +inf[

final float missibility = 1.f;               // [0,1]      : 0=>fully missible ; ]0,1[=>partialy missible ; 1=>immissible
final float interfacial_thickness = 20.f;    // ]0, +inf[

// ----------------------------------------------------- LBM CONST -----------------------------------------------------
final float cs = 0.57735027f;
final float cs2 = 0.33333333f;
final float cs4 = 0.11111111f;
final float cs6 = 0.03703704f;
  
final float[] D2Q9_w = { 4.f/9.f, 1.f/9.f, 1.f/9.f, 1.f/9.f, 1.f/9.f, 1.f/36.f, 1.f/36.f, 1.f/36.f, 1.f/36.f };
final int[] D2Q9_cx = { 0, 1,-1, 0, 0, 1,-1, 1,-1 }; 
final int[] D2Q9_cy = { 0, 0, 0, 1,-1, 1,-1,-1, 1 };

final float[] D2Q5_w = { 1.f/3.f, 1.f/6.f, 1.f/6.f, 1.f/6.f, 1.f/6.f };
final int[] D2Q5_cx = { 0, 1,-1, 0, 0}; 
final int[] D2Q5_cy = { 0, 0, 0, 1,-1};

// ------------------------------------------------------- UTILS -------------------------------------------------------
int mod(int x, int n) { return (x+n)%n; }
float cb(float x) { return x*x*x; }

boolean inSphere(int x, int y, int r, int cx, int cy){
  return ((cx-x)*(cx-x)+(cy-y)*(cy-y)) < (r*r);
}
