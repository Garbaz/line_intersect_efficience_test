
void setup() {
  ArrayList<Integer> times = new ArrayList<Integer>();
  for (int i = 0; i < 100; i++) {
    PVector a = PVector.random2D(), b = PVector.random2D(), c = PVector.random2D(), d = PVector.random2D();
    //PVector a = PVector.random3D(), b = PVector.random3D(), c = PVector.random3D(), d = PVector.random3D(), e = PVector.random3D();
    int start_time = millis();
    for (int j = 0; j < 10000000; j++) {
      intersectLine_TOBI_IMPR2(a, b, c, d);
      //line_plane_intersect(a, b, c, d, e);
    }
    times.add(millis()-start_time);
  }
  float avg = 0;
  for (int t : times) {
    avg += t;
  }
  avg /= times.size();
  println(avg);
}

PVector intersectLine_TOBI_IMPR2(PVector l1e1, PVector l1e2, PVector l2e1, PVector l2e2) {
  final float EPS = 1e-4;

  float a1, b1;
  float delta_l1_x = l1e2.x-l1e1.x;
  float delta_l1_y = l1e1.y-l1e2.y;
  if (abs(delta_l1_x/delta_l1_y) < 1) { // No abs on d_l1_x because of the swap above
    a1 = 1;
    b1 = delta_l1_x/delta_l1_y;
  } else {
    a1 = delta_l1_y/delta_l1_x;
    b1 = 1;
  }

  float a2, b2;
  float delta_l2_x = l2e2.x-l2e1.x;
  float delta_l2_y = l2e1.y-l2e2.y;
  if (abs(delta_l2_x/delta_l2_y) < 1) { // No abs on d_l2_x because of the swap above
    a2 = 1;
    b2 = delta_l2_x/delta_l2_y;
  } else {
    a2 = delta_l2_y/delta_l2_x;
    b2 = 1;
  }

  float k1 = a1*l1e1.x+b1*l1e1.y, k2 = a2*l2e1.x+b2*l2e1.y;
  float d = a1*b2-a2*b1;

  if (abs(d) < EPS) { // Parallel or almost parallel lines
    return null;
  }

  PVector sol = new PVector((b2*k1-b1*k2)/d, (a1*k2-a2*k1)/d);

  if (min(l1e1.x, l1e2.x) <= sol.x && sol.x <= max(l1e1.x, l1e2.x) && min(l2e1.x, l2e2.x) <= sol.x && sol.x <= max(l2e1.x, l2e2.x)) { // Check that the solution actually is on the line segments
    return sol;
  } else {
    return null;
  }
}

PVector intersectLine_TOBI_IMPR(PVector l1e1, PVector l1e2, PVector l2e1, PVector l2e2) {
  final float EPS = 1e-4;

  //PVector l1e1 = line1.end1, l1e2 = line1.end2, l2e1 = line2.end1, l2e2 = line2.end2;

  if (l1e1.x > l1e2.x) { // Swapping the ends around such that l1e1.x < l1e2.x (optimization)
    float tmp_x = l1e1.x, tmp_y = l1e1.y;
    l1e1.set(l1e2.x, l1e2.y);
    l1e2.set(tmp_x, tmp_y);
  }
  if (l2e1.x > l2e2.x) { // Swapping the ends around such that l2e1.x < l2e2.x (optimization)
    float tmp_x = l2e1.x, tmp_y = l2e1.y;
    l2e1.set(l2e2.x, l2e2.y);
    l2e2.set(tmp_x, tmp_y);
  }

  float a1, b1, a2, b2;
  float d_l1_x = l1e2.x-l1e1.x;
  float d_l1_y = l1e1.y-l1e2.y;
  if (d_l1_x < abs(d_l1_y)) { // No abs on d_l1_x because of the swap above
    a1 = 1;
    b1 = d_l1_x/d_l1_y;
  } else {
    a1 = d_l1_y/d_l1_x;
    b1 = 1;
  }
  float d_l2_x = l2e2.x-l2e1.x;
  float d_l2_y = l2e1.y-l2e2.y;
  if (d_l2_x < abs(d_l2_y)) { // No abs on d_l2_x because of the swap above
    a2 = 1;
    b2 = d_l2_x/d_l2_y;
  } else {
    a2 = d_l2_y/d_l2_x;
    b2 = 1;
  }

  float k1 = a1*l1e1.x+b1*l1e1.y, k2 = a2*l2e1.x+b2*l2e1.y;
  float d = a1*b2-a2*b1;

  if (abs(d) < EPS) { // Parallel or almost parallel lines
    return null;
  }

  PVector sol = new PVector((b2*k1-b1*k2)/d, (a1*k2-a2*k1)/d);

  if (l1e1.x <= sol.x && sol.x <= l1e2.x && l2e1.x <= sol.x && sol.x <= l2e2.x) {
    return sol;
  } else {
    return null;
  }
}

PVector intersectLine_TOBI(PVector l1e1, PVector l1e2, PVector l2e1, PVector l2e2) {
  final float EPS = 1e-4;

  //PVector l1e1 = line1.end1, l1e2 = line1.end2, l2e1 = line2.end1, l2e2 = line2.end2;

  float a1, b1, a2, b2;
  if (abs(l1e1.x - l1e2.x) < abs(l1e1.y - l1e2.y)) {
    a1 = 1;
    b1 = (l1e2.x-l1e1.x)/(l1e1.y-l1e2.y);
  } else {
    a1 = (l1e1.y-l1e2.y)/(l1e2.x-l1e1.x);
    b1 = 1;
  }

  if (abs(l2e1.x - l2e2.x) < abs(l2e1.y - l2e2.y)) {
    a2 = 1;
    b2 = (l2e2.x-l2e1.x)/(l2e1.y-l2e2.y);
  } else {
    a2 = (l2e1.y-l2e2.y)/(l2e2.x-l2e1.x);
    b2 = 1;
  }

  float k1 = a1*l1e1.x+b1*l1e1.y, k2 = a2*l2e1.x+b2*l2e1.y;
  float d = a1*b2-a2*b1;

  if (abs(d) < EPS) { // Parallel or almost parallel lines
    return null;
  }

  PVector sol = new PVector((b2*k1-b1*k2)/d, (a1*k2-a2*k1)/d);

  PVector delta1 = PVector.sub(l1e2, l1e1), delta2 = PVector.sub(l2e2, l2e1);
  PVector deltasol1 = PVector.sub(sol, l1e1), deltasol2 = PVector.sub(sol, l2e1);

  if (delta1.dot(deltasol1) > 0 && deltasol1.magSq() < delta1.magSq() && delta2.dot(deltasol2) > 0 && deltasol2.magSq() < delta2.magSq()) {
    return sol;
  } else {
    return null;
  }
  //return sol;
}

PVector intersectLine_FLO(PVector lPos1, PVector lPos2, PVector rayPos, PVector rayDir) {
  PVector lDiff = PVector.sub(lPos2, lPos1);
  PVector cross;
  float lM, lC, crossX;
  //text(rayDir.x, 100, 200);
  if (lDiff.x == 0 || rayDir.x == 0) {
    lM = lDiff.x / lDiff.y;
    lC = lPos1.x - lPos1.y * lM;


    float rayM = rayDir.x / rayDir.y;
    float rayC = rayPos.x - rayPos.y * rayM;

    crossX = (lC - rayC)/(rayM-lM);
    cross = new PVector( rayM *crossX +rayC, crossX);
    if (abs((cross.y-lPos1.y) + (cross.y-lPos2.y)) <= abs(lDiff.y) && abs((cross.x-rayPos.x) + (cross.x-(rayPos.x + rayDir.x))) <= abs(rayDir.x)) {
      return cross;
    }
  } else {

    lM = lDiff.y / lDiff.x;
    lC = lPos1.y - lPos1.x * lM;


    float rayM = rayDir.y / rayDir.x;
    float rayC = rayPos.y - rayPos.x * rayM;

    crossX = (lC - rayC)/(rayM-lM);
    cross = new PVector(crossX, rayM *crossX +rayC);
    if (abs((cross.x-lPos1.x) + (cross.x-lPos2.x)) <= abs(lDiff.x) && abs((cross.y-rayPos.y) + (cross.y-(rayPos.y + rayDir.y))) <= abs(rayDir.y)) {
      return cross;
    }
  }

  return null;
}

//public PVector intersectLine_TOBI_RAY(float x0, float y0, float x1, float y1, PVector anchor) {

//  float a = cos(angle+HALF_PI), b = sin(angle+HALF_PI);
//  float a_, b_;
//  if (abs(x0 - x1) < abs(y1 - y0)) {
//    a_ = 1;
//    b_ = (x1-x0)/(y0-y1);
//  } else {
//    a_ = (y0-y1)/(x1-x0);
//    b_ = 1;
//  }

//  float k = a*anchor.x+b*anchor.y, k_ = a_*x0+b_*y0;
//  float d = a*b_-b*a_;

//  if (abs(d) < eps) { // Parallel or almost parallel lines
//    return null;
//  }

//  float x_sol = (b_*k-b*k_)/d, y_sol = (a*k_-a_*k)/d;

//  if ( abs(quadrantAngle(y_sol-anchor.y, x_sol-anchor.x)-angle) <= HALF_PI && min(x0, x1)-x_sol < eps && x_sol-max(x0, x1) < eps && min(y0, y1)-y_sol < eps && y_sol-max(y0, y1) < eps) {
//    return new PVector(x_sol, y_sol);
//  } else {
//    return null;
//  }
//}

//public PVector intersectLine(PVector p0, PVector p1) {
//  return intersectLine(p0.x, p0.y, p1.x, p1.y);
//}

PVector intersectLine_FLO_BETTER(PVector line1_Pos1, PVector line1_Pos2, PVector line2_Pos1, PVector line2_Pos2) {
  PVector line1_Diff = PVector.sub(line1_Pos2, line1_Pos1);
  PVector line2_Diff = PVector.sub(line2_Pos2, line2_Pos1);
  PVector cross;
  if (line1_Diff.x != 0 && line2_Diff.x != 0) {
    float line1M = line1_Diff.y / line1_Diff.x;
    float line1C = line1_Pos1.y - line1_Pos1.x * line1M;


    float line2M = line2_Diff.y / line2_Diff.x;
    float line2C = line2_Pos1.y - line2_Pos1.x * line2M;

    float crossX = (line1C - line2C)/(line2M-line1M);
    cross = new PVector(crossX, line2M *crossX + line2C);
  } else if (line1_Diff.x == 0 && line2_Diff.x != 0) {
    float line2M = line2_Diff.y / line2_Diff.x;
    cross = new PVector(line1_Pos1.x, line1_Pos1.x * line2M + (line2_Pos1.y - line2_Pos1.x * line2M));
  } else if (line2_Diff.x == 0 && line1_Diff.x != 0) {
    float line1M = line1_Diff.y / line1_Diff.x;
    cross = new PVector(line2_Pos1.x, line2_Pos1.x * line1M + (line1_Pos1.y - line1_Pos1.x * line1M));
  } else {
    return null;
  }

  if (isBetween(cross.x, line1_Pos1.x, line1_Pos2.x) && isBetween(cross.y, line2_Pos1.y, line2_Pos2.y) && isBetween(cross.y, line1_Pos1.y, line1_Pos2.y) && isBetween(cross.x, line2_Pos1.x, line2_Pos2.x)) {
    return cross;
  } else return null;
}
boolean isBetween(float point, float linePos1, float linePos2) {
  return abs((point-linePos1) + (point-linePos2)) <= abs(linePos2-linePos1);
}



// #####################################################################################################
// ########################################## LINE AND PLANE ###########################################
// #####################################################################################################

PVector line_plane_intersect(PVector bl, PVector l, PVector bp, PVector p1, PVector p2) {
  PVector b = PVector.sub(bl, bp);

  float u=-((b.x)*((l.z)*(p2.y)-(l.y)*(p2.z))+(l.x)*((b.y)*(p2.z)-(b.z)*(p2.y))+((b.z)*(l.y)-(b.y)*(l.z))*(p2.x))/((l.x)*((p1.z)*(p2.y)-(p1.y)*(p2.z))+(p1.x)*((l.y)*(p2.z)-(l.z)*(p2.y))+((l.z)*(p1.y)-(l.y)*(p1.z))*(p2.x)); 
  float v=((b.x)*((l.z)*(p1.y)-(l.y)*(p1.z))+(l.x)*((b.y)*(p1.z)-(b.z)*(p1.y))+((b.z)*(l.y)-(b.y)*(l.z))*(p1.x))/((l.x)*((p1.z)*(p2.y)-(p1.y)*(p2.z))+(p1.x)*((l.y)*(p2.z)-(l.z)*(p2.y))+((l.z)*(p1.y)-(l.y)*(p1.z))*(p2.x));
  float t=-((b.x)*((p1.z)*(p2.y)-(p1.y)*(p2.z))+(p1.x)*((b.y)*(p2.z)-(b.z)*(p2.y))+((b.z)*(p1.y)-(b.y)*(p1.z))*(p2.x))/((l.x)*((p1.z)*(p2.y)-(p1.y)*(p2.z))+(p1.x)*((l.y)*(p2.z)-(l.z)*(p2.y))+((l.z)*(p1.y)-(l.y)*(p1.z))*(p2.x));

  if (0 <= t && t <= 1 && 0 <= u && 0 <= v && 0 <= u+v && u+v <= 1) {
    PVector ret = l;
    ret.mult(t);
    ret.add(bl);
    return ret;
  } else {
    return null;
  }
}
