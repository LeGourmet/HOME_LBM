// ------------------------------------------------------- UTILS -------------------------------------------------------
int mod(int x, int n) { return (x+n)%n; }

float cb(float a) { return a*a*a; }

float cbrt(float a) { return pow(a,1.f/3.f); }

float sign(float f) {
  if (f > 0.f) return 1.f;
  if (f < 0.f) return -1.f;
  return 0.f;
} 

boolean inSphere(int x, int y, int r, int cx, int cy){
  return ((cx-x)*(cx-x)+(cy-y)*(cy-y)) < (r*r);
}

// --------------------------------------------------- FLUIDS CONST ----------------------------------------------------
final float nu_fluid = 0.0001f;              // Desbrun : [0.01f, 0.0006f] => 0.0015f
final float st_fluid = 4e-3f;               // surface tension soulde be 10^-6surface tension soulde be 10^-6

// ----------------------------------------------------- LBM CONST -----------------------------------------------------
final float cs = 0.57735027f;
final float cs2 = 0.33333333f;
final float cs4 = 0.11111111f;
final float cs6 = 0.03703704f;
  
final float[] D2Q9_w = { 4.f/9.f, 1.f/9.f, 1.f/9.f, 1.f/9.f, 1.f/9.f, 1.f/36.f, 1.f/36.f, 1.f/36.f, 1.f/36.f };
final int[] D2Q9_cx = { 0, 1,-1, 0, 0, 1,-1, 1,-1 }; 
final int[] D2Q9_cy = { 0, 0, 0, 1,-1, 1,-1,-1, 1 };

// ----------------------------------------------------- LBM FUNCTIONS -----------------------------------------------------

float plic_cube_reduced(float V, float n1, float n2, float n3) { // optimized solution from SZ and Kawano, source: https://doi.org/10.3390/computation10020021
  float n12 = n1 + n2, n3V = n3 * V;
  if (n12 <= 2.0f * n3V) return n3V + 0.5f * n12; // case (5)
  float sqn1 = sq(n1), n26 = 6.0f * n2, v1 = sqn1 / n26; // after case (5) check n2>0 is true
  if (v1 <= n3V && n3V < v1 + 0.5f * (n2 - n1)) return 0.5f * (n1 + sqrt(sqn1 + 8.0f * n2 * (n3V - v1))); // case (2)
  float V6 = n1 * n26 * n3V;
  if (n3V < v1) return cbrt(V6); // case (1)
  float v3 = n3 < n12 ? (sq(n3) * (3.0f * n12 - n3) + sqn1 * (n1 - 3.0f * n3) + sq(n2) * (n2 - 3.0f * n3)) / (n1 * n26) : 0.5f * n12; // after case (2) check n1>0 is true
  float sqn12 = sqn1 + sq(n2), V6cbn12 = V6 - cb(n1) - cb(n2);
  boolean case34 = n3V < v3; // true: case (3), false: case (4)
  float a = case34 ? V6cbn12 : 0.5f * (V6cbn12 - cb(n3));
  float b = case34 ? sqn12 : 0.5f * (sqn12 + sq(n3));
  float c = case34 ? n12 : 0.5f;
  float t = sqrt(sq(c) - b);
  return c - 2.f * t * sin(0.33333334f * asin((cb(c) - 0.5f * a - 1.5f * b * c) / cb(t)));
}

float plic_cube(float V0, PVector n) { // unit cube - plane intersection: volume V0 in [0,1], normal vector n -> plane offset d0
  float ax = abs(n.x), ay = abs(n.y), az = abs(n.z), V = 0.5f - abs(V0 - 0.5f), l = ax + ay + az; // eliminate symmetry cases, normalize n using L1 norm
  float n1 = min(min(ax, ay), az) / l;
  float n3 = max(max(ax, ay), az) / l;
  float n2 = max(0.f, 1.f-(n1+n3)); // ensure n2>=0 (original : float n2 = fdim(1.0f, n1 + n3);) 
  float d = plic_cube_reduced(V, n1, n2, n3); // calculate PLIC with reduced symmetry
  return l * abs(0.5f - d) * sign(V0 - 0.5f); // rescale result and apply symmetry for V0>0.5 (original : return l * copysign(0.5f - d, V0 - 0.5f);) 
}

// rsqrt, calculate_normal_py
float calculate_curvature(float[] phit) {
  PVector by = new PVector(2.f*(phit[2]-phit[1])+phit[6]-phit[5]+phit[8]-phit[7], 2.f*(phit[4]-phit[3])+phit[6]-phit[5]+phit[7]-phit[8], 0.f).normalize(); // new coordinate system: bz is normal to surface, bx and by are tangent to surface
  PVector bx = by.cross(new PVector(0.0f, 0.0f, 1.0f)); // normalize() is necessary here because bz and rn are not perpendicular
  int number = 0; // number of neighboring interface points
  PVector[] p = new PVector[6]; // number of neighboring interface points is less or equal than than 8 minus 1 gas and minus 1 fluid point = 6
  float center_offset = plic_cube(phit[0], by); // calculate z-offset PLIC of center point only once
  for(int i=1; i<9; i++) { // iterate over neighbors, no loop unrolling here (50% better perfoemance without loop unrolling)
    if(phit[i]>0.f && phit[i]<1.f) { // limit neighbors to interface cells
      PVector ei = new PVector(D2Q9_cx[i], D2Q9_cy[i], 0.f); // assume neighbor normal vector is the same as center normal vector
      float offset = plic_cube(phit[i], by) - center_offset;
      p[number++] = new PVector(PVector.dot(ei, bx), PVector.dot(ei, by) + offset); // do coordinate system transformation into (x, f(x)) and apply PLIC pffsets
    }
  }
  float[] M = {0.0f,0.0f,0.0f,0.0f};
  float[] b = {0.0f,0.0f};
  for(int i=0; i<number; i++) { // f(x,y)=A*x2+H*x, x=(A,H), Q=(x2,x), M*x=b, M=Q*Q^T, b=Q*z
    float x=p[i].x, y=p[i].y, x2=x*x, x3=x2*x;
    /**/M[0]+=x2*x2;   M[1]+=x3; b[0]+=x2*y;
    /*  M[2]+=x3   ;*/ M[3]+=x2; b[1]+=x *y;
  }
  M[2] = M[1]; // use symmetry of matrix to save arithmetic operations
  
  float[] x = {0.0f,0.0f};
  int N = 2, Nsol = min(2, number);
  for(int i=0; i<Nsol ;i++) { // decompose M in M=L*U
    for(int j=i+1; j<Nsol ;j++) {
      M[N*j+i] /= M[N*i+i];
      for(int k=i+1; k<Nsol ;k++) M[N*j+k] -= M[N*j+i] * M[N*i+k];
    }
  }
  for(int i=0; i<Nsol ;i++) { // find solution of L*y=b
    x[i] = b[i];
    for(int k=0; k<i; k++) x[i] -= M[N * i + k] * x[k];
  }
  for (int i = Nsol - 1; i >= 0; i--) { // find solution of U*x=y
    for (int k = i + 1; k < Nsol; k++) x[i] -= M[N * i + k] * x[k];
    x[i] /= M[N * i + i];
  }
  
  float A=x[0], H=x[1];
  float K = 2.f*A*cb(1.f/sqrt(H*H+1.f)); // mean curvature of Monge patch (x, f(x)), note that curvature definition in 2D is different than 3D (additional factor 2)
  return constrain(K, -1.f, 1.f); // prevent extreme pressures in the case of almost degenerate matrices
}
