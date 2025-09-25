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

boolean inSphere(int x, int y, int r, int cx, int cy){
  return ((cx-x)*(cx-x)+(cy-y)*(cy-y)) < (r*r);
}

boolean inTorus(int x, int y, int r1, int r2, int cx, int cy){
  float mag = (cx-x)*(cx-x)+(cy-y)*(cy-y);
  return mag>(r1*r1) && mag<(r2*r2);
}
