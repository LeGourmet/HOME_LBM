// --------------------------------------------------- FLUIDS CONST ----------------------------------------------------
// should use 4 time finer discreate field for phi (2 per directions)
// should update hi at half time fi(t) => fi(t+1) ; hi(t) => hi(t+0.5) => hi(t+1)
// should update solid for moving solid boundary treatment
// LES is not smagorinsky lilly but WALE (nu = nu + nu')
// abs(ca_fluid - ca_air) < cs/1000.f 
// mobility : controling degree of interface splitting; smaller values imply stronger spliting as less diffusion is intoducted
// g = 1e-5f;
// surface tension water = 10-6;
// outlet condition should use 
// phi(x) = 1/2*(1-tanh(2x/interface_size) ou x est egale a sa distance signé de l'interfaces des phases 

final float rho_fluid = 3.f;                // 0.1f => idk !
final float nu_fluid = 0.001f;              // Desbrun : [0.01f, 0.0006f]
final float mo_fluid = 0.2f;                // Desbrun : 0.2f
final float ca_fluid = 0.26f;                // 0.00005f => idk !

final float rho_air = 1.f;                  // 0.00001f => idk !
final float nu_air = 0.001f;                // Desbrun : [0.01f, 0.0006f]
final float mo_air = 0.2f;                  // Desbrun : 0.2f
final float ca_air = 0.26f;                 // 0.0005f => idk !

final float interfacial_thickness = 5.f;    // Desbrun : 5.f

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

boolean inTorus(int x, int y, int r1, int r2, int cx, int cy){
  float mag = (cx-x)*(cx-x)+(cy-y)*(cy-y);
  return mag>(r1*r1) && mag<(r2*r2);
}
