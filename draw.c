#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "ml6.h"
#include "display.h"
#include "draw.h"
#include "matrix.h"
#include "math.h"
#include "gmath.h"

/*======== void scanline_convert() ==========
  Inputs: struct matrix *points
          int i
          screen s
          zbuffer zb
  Returns:

  Fills in polygon i by drawing consecutive horizontal (or vertical) lines.

  Color should be set differently for each polygon.
  ====================*/
void scanline_convert( struct matrix *points, int i, screen s, zbuffer zb ) {

}

/*======== void add_polygon() ==========
  Inputs:   struct matrix *polygons
            double x0
            double y0
            double z0
            double x1
            double y1
            double z1
            double x2
            double y2
            double z2
  Returns:
  Adds the vertices (x0, y0, z0), (x1, y1, z1)
  and (x2, y2, z2) to the polygon matrix. They
  define a single triangle surface.
  ====================*/
void add_polygon( struct matrix *polygons,
                  double x0, double y0, double z0,
                  double x1, double y1, double z1,
                  double x2, double y2, double z2 ) {
  add_point(polygons, x0, y0, z0);
  add_point(polygons, x1, y1, z1);
  add_point(polygons, x2, y2, z2);
}

/*======== void draw_polygons() ==========
  Inputs:   struct matrix *polygons
            screen s
            color c
  Returns:
  Goes through polygons 3 points at a time, drawing
  lines connecting each points to create bounding triangles
  ====================*/
void draw_polygons( struct matrix *polygons, screen s, zbuffer zb, color c ) {
  printf("------------------------------------------\n");
  if ( polygons->lastcol < 3 ) {
    printf("Need at least 3 points to draw a polygon!\n");
    return;
  }

  int point;
  double *normal;

  // 0 to 1257
  for (point=0; point < polygons->lastcol - 2; point+=3) {

    normal = calculate_normal(polygons, point);

    if ( normal[2] > 0 ) {

      c.red = rand() % 256;
      c.green = rand() % 256;
      c.blue = rand() % 256;
      // c.red = point * 3;
      // c.green = point * 4;
      // c.blue = point * 5;

      if (point > polygons->lastcol / 4) {
        c.red = 0;
        c.blue = 128;
      }


      int y0;
      int y1;
      int y2;
      int bottom;
      int middle;
      int top;
      int temp;

      y0 = polygons->m[1][point];
      y1 = polygons->m[1][point+1];
      y2 = polygons->m[1][point+2];

      if (y0 >= y1 && y0 >= y2) {
        top = point; // top point is the one that corresponds to y0 in this case
        middle = point + 1; // assume y1 is middle, change later if wrong
        bottom = point + 2;
      }
      else if (y1 >= y0 && y1 >= y2) {
        top = point + 1;
        middle = point;
        bottom = point + 2;
      }
      else {
        top = point + 2;
        middle = point + 1;
        bottom = point;
      }

      if (polygons->m[1][bottom] > polygons->m[1][middle]) {
        temp = bottom;
        bottom = middle;
        middle = temp;
      }

      int x_bottom, x_mid, x_top, y_bottom, y_mid, y_top, z_bottom, z_mid, z_top;
      double x0, x1, y, z0, z1, dx0, dx1, dx1_1, dz0, dz1, dz1_1;

      x_bottom = polygons->m[0][bottom];
      x_mid = polygons->m[0][middle];
      x_top = polygons->m[0][top];

      y_bottom = polygons->m[1][bottom];
      y_mid = polygons->m[1][middle];
      y_top = polygons->m[1][top];

      z_bottom = polygons->m[2][bottom];
      z_mid = polygons->m[2][middle];
      z_top = polygons->m[2][top];

      printf("Point: %d\n", point);
      printf("(%d, %d, %d)  (%d, %d, %d)  (%d, %d, %d)\n",
               x_bottom, y_bottom, z_bottom, x_mid, y_mid, z_mid, x_top, y_top, z_top);

      x0 = x_bottom;
      z0 = z_bottom;

      x1 = x_bottom;
      z1 = z_bottom;

      y = y_bottom;

      dx0 = (double)(x_top - x_bottom) / (double)(y_top - y_bottom);
      dz0 = (double)(z_top - z_bottom) / (double)(y_top - y_bottom);

      if (y_bottom == y_mid) {
        x1 = x_mid;
        z1 = z_mid;

        dx1_1 = (x_top - x_mid) / (double)(y_top - y_mid);
        dx1 = dx1_1;

        dz1_1 = (z_top - z_mid) / (double)(y_top - y_mid);
        dz1 = dz1_1;
      }

      else if (y_mid == y_top) {
        dx1 = (x_mid - x_bottom) / (double)(y_mid - y_bottom);
        dx1_1 = dx1;

        dz1 = (z_mid - z_bottom) / (double)(y_mid - y_bottom);
        dz1_1 = dz1;
      }

      else {
        dx1 = (x_mid - x_bottom) / (double)(y_mid - y_bottom);
        dx1_1 = (x_top - x_mid) / (double)(y_top - y_mid);

        dz1 = (z_mid - z_bottom) / (double)(y_mid - y_bottom);
        dz1_1 = (z_top - z_mid) / (double)(y_top - y_mid);
      }
      //
      // if (point == 36) {
      //   printf("%d %lf %d\n", y_bottom, y, y_top);
      // }
      // printf("NEW POLYGON\n");
      for (y = y_bottom; y <= y_top; y++, x0 += dx0, x1 += dx1, z0 += dz0, z1 += dz1) {
          printf("(%lf, %lf, %lf) to (%lf, %lf, %lf)\n", x0, y, z0, x1, y, z1);
          draw_line(x0, y, z0, x1, y, z1, s, zb, c);
          if (y >= y_mid) {
              dx1 = dx1_1;
              dz1 = dz1_1;
            }
      }


    }
  }
}

/*======== void add_box() ==========
  Inputs:   struct matrix * edges
            double x
            double y
            double z
            double width
            double height
            double depth

  add the points for a rectagular prism whose
  upper-left-front corner is (x, y, z) with width,
  height and depth dimensions.
  ====================*/
void add_box( struct matrix *polygons,
              double x, double y, double z,
              double width, double height, double depth ) {
  double x0, y0, z0, x1, y1, z1;
  x0 = x;
  x1 = x+width;
  y0 = y;
  y1 = y-height;
  z0 = z;
  z1 = z-depth;


  //front
  add_polygon(polygons, x, y, z, x1, y1, z, x1, y, z);
  add_polygon(polygons, x, y, z, x, y1, z, x1, y1, z);
  //back
  add_polygon(polygons, x1, y, z1, x, y1, z1, x, y, z1);
  add_polygon(polygons, x1, y, z1, x1, y1, z1, x, y1, z1);

  //right side
  add_polygon(polygons, x1, y, z, x1, y1, z1, x1, y, z1);
  add_polygon(polygons, x1, y, z, x1, y1, z, x1, y1, z1);
  //left side
  add_polygon(polygons, x, y, z1, x, y1, z, x, y, z);
  add_polygon(polygons, x, y, z1, x, y1, z1, x, y1, z);

  //top
  add_polygon(polygons, x, y, z1, x1, y, z, x1, y, z1);
  add_polygon(polygons, x, y, z1, x, y, z, x1, y, z);
  //bottom
  add_polygon(polygons, x, y1, z, x1, y1, z1, x1, y1, z);
  add_polygon(polygons, x, y1, z, x, y1, z1, x1, y1, z1);
}


/*======== void add_sphere() ==========
  Inputs:   struct matrix * points
            double cx
            double cy
            double cz
            double r
            int step

  adds all the points for a sphere with center (cx, cy, cz)
  and radius r using step points per circle/semicircle.

  Since edges are drawn using 2 points, add each point twice,
  or add each point and then another point 1 pixel away.

  should call generate_sphere to create the necessary points
  ====================*/
void add_sphere( struct matrix * edges,
                 double cx, double cy, double cz,
                 double r, int step ) {

  struct matrix *points = generate_sphere(cx, cy, cz, r, step);
  int p0, p1, p2, p3, lat, longt;
  int latStop, longStop, latStart, longStart;
  latStart = 0;
  latStop = step;
  longStart = 1;
  longStop = step;

  //step++; needed for my triangles
  for ( lat = latStart; lat < latStop; lat++ ) {
    for ( longt = longStart; longt < longStop; longt++ ) {

      /*Milan's Triangles*/
      p0 = lat * (step+1) + longt;
      p1 = p0 + 1;
      p2 = (p1 + step) % (step * (step+1));
      p3 = (p0 + step) % (step * (step+1));

      add_polygon( edges, points->m[0][p0],
                   points->m[1][p0],
                   points->m[2][p0],
                   points->m[0][p1],
                   points->m[1][p1],
                   points->m[2][p1],
                   points->m[0][p2],
                   points->m[1][p2],
                   points->m[2][p2]);
      add_polygon( edges, points->m[0][p0],
                   points->m[1][p0],
                   points->m[2][p0],
                   points->m[0][p2],
                   points->m[1][p2],
                   points->m[2][p2],
                   points->m[0][p3],
                   points->m[1][p3],
                   points->m[2][p3]);

      /*My Triangles*/
      /* p0 = lat * (step) + longt; */
      /* p1 = p0+1; */
      /* p2 = (p1+step) % (step * (step-1)); */
      /* p3 = (p0+step) % (step * (step-1)); */

      /* //printf("p0: %d\tp1: %d\tp2: %d\tp3: %d\n", p0, p1, p2, p3); */
      /* if (longt < step - 2) */
      /*   add_polygon( edges, points->m[0][p0], */
      /*                points->m[1][p0], */
      /*                points->m[2][p0], */
      /*                points->m[0][p1], */
      /*                points->m[1][p1], */
      /*                points->m[2][p1], */
      /*                points->m[0][p2], */
      /*                points->m[1][p2], */
      /*                points->m[2][p2]); */
      /* if (longt > 0 ) */
      /*   add_polygon( edges, points->m[0][p0], */
      /*                points->m[1][p0], */
      /*                points->m[2][p0], */
      /*                points->m[0][p2], */
      /*                points->m[1][p2], */
      /*                points->m[2][p2], */
      /*                points->m[0][p3], */
      /*                points->m[1][p3], */
      /*                points->m[2][p3]); */
    }
  }
  free_matrix(points);
}

/*======== void generate_sphere() ==========
  Inputs:   struct matrix * points
            double cx
            double cy
            double cz
            double r
            int step
  Returns: Generates all the points along the surface
           of a sphere with center (cx, cy, cz) and
           radius r using step points per circle/semicircle.
           Returns a matrix of those points
  ====================*/
struct matrix * generate_sphere(double cx, double cy, double cz,
                                double r, int step ) {

  struct matrix *points = new_matrix(4, step * step);
  int circle, rotation, rot_start, rot_stop, circ_start, circ_stop;
  double x, y, z, rot, circ;

  rot_start = 0;
  rot_stop = step;
  circ_start = 0;
  circ_stop = step;

  for (rotation = rot_start; rotation < rot_stop; rotation++) {
    rot = (double)rotation / step;

    for(circle = circ_start; circle <= circ_stop; circle++){
      circ = (double)circle / step;

      x = r * cos(M_PI * circ) + cx;
      y = r * sin(M_PI * circ) *
        cos(2*M_PI * rot) + cy;
      z = r * sin(M_PI * circ) *
        sin(2*M_PI * rot) + cz;

      /* printf("rotation: %d\tcircle: %d\n", rotation, circle); */
      /* printf("rot: %lf\tcirc: %lf\n", rot, circ); */
      /* printf("sphere point: (%0.2f, %0.2f, %0.2f)\n\n", x, y, z); */
      add_point(points, x, y, z);
    }
  }

  return points;
}

/*======== void add_torus() ==========
  Inputs:   struct matrix * points
            double cx
            double cy
            double cz
            double r1
            double r2
            double step
  Returns:

  adds all the points required for a torus with center (cx, cy, cz),
  circle radius r1 and torus radius r2 using step points per circle.

  should call generate_torus to create the necessary points
  ====================*/
void add_torus( struct matrix * edges,
                double cx, double cy, double cz,
                double r1, double r2, int step ) {

  struct matrix *points = generate_torus(cx, cy, cz, r1, r2, step);
  int p0, p1, p2, p3, lat, longt;
  int latStop, longStop, latStart, longStart;
  latStart = 0;
  latStop = step;
  longStart = 0;
  longStop = step;

  for ( lat = latStart; lat < latStop; lat++ ) {
    for ( longt = longStart; longt < longStop; longt++ ) {
      p0 = lat * step + longt;
      if (longt == step - 1)
        p1 = p0 - longt;
      else
        p1 = p0 + 1;
      p2 = (p1 + step) % (step * step);
      p3 = (p0 + step) % (step * step);

      //printf("p0: %d\tp1: %d\tp2: %d\tp3: %d\n", p0, p1, p2, p3);
      add_polygon( edges, points->m[0][p0],
                   points->m[1][p0],
                   points->m[2][p0],
                   points->m[0][p3],
                   points->m[1][p3],
                   points->m[2][p3],
                   points->m[0][p2],
                   points->m[1][p2],
                   points->m[2][p2]);
      add_polygon( edges, points->m[0][p0],
                   points->m[1][p0],
                   points->m[2][p0],
                   points->m[0][p2],
                   points->m[1][p2],
                   points->m[2][p2],
                   points->m[0][p1],
                   points->m[1][p1],
                   points->m[2][p1]);
    }
  }
  free_matrix(points);
}

/*======== void generate_torus() ==========
  Inputs:   struct matrix * points
            double cx
            double cy
            double cz
            double r
            int step
  Returns: Generates all the points along the surface
           of a torus with center (cx, cy, cz),
           circle radius r1 and torus radius r2 using
           step points per circle.
           Returns a matrix of those points
  ====================*/
struct matrix * generate_torus( double cx, double cy, double cz,
                                double r1, double r2, int step ) {

  struct matrix *points = new_matrix(4, step * step);
  int circle, rotation, rot_start, rot_stop, circ_start, circ_stop;
  double x, y, z, rot, circ;

  rot_start = 0;
  rot_stop = step;
  circ_start = 0;
  circ_stop = step;

  for (rotation = rot_start; rotation < rot_stop; rotation++) {
    rot = (double)rotation / step;

    for(circle = circ_start; circle < circ_stop; circle++){
      circ = (double)circle / step;

      x = cos(2*M_PI * rot) *
        (r1 * cos(2*M_PI * circ) + r2) + cx;
      y = r1 * sin(2*M_PI * circ) + cy;
      z = -1*sin(2*M_PI * rot) *
        (r1 * cos(2*M_PI * circ) + r2) + cz;

      //printf("rotation: %d\tcircle: %d\n", rotation, circle);
      //printf("torus point: (%0.2f, %0.2f, %0.2f)\n", x, y, z);
      add_point(points, x, y, z);
    }
  }
  return points;
}

/*======== void add_circle() ==========
  Inputs:   struct matrix * edges
            double cx
            double cy
            double r
            double step

  Adds the circle at (cx, cy) with radius r to edges
  ====================*/
void add_circle( struct matrix *edges,
                 double cx, double cy, double cz,
                 double r, int step ) {
  double x0, y0, x1, y1, t;
  int i;

  x0 = r + cx;
  y0 = cy;
  for (i=1; i<=step; i++) {
    t = (double)i/step;
    x1 = r * cos(2 * M_PI * t) + cx;
    y1 = r * sin(2 * M_PI * t) + cy;

    add_edge(edges, x0, y0, cz, x1, y1, cz);
    x0 = x1;
    y0 = y1;
  }
}


/*======== void add_curve() ==========
Inputs:   struct matrix *edges
         double x0
         double y0
         double x1
         double y1
         double x2
         double y2
         double x3
         double y3
         double step
         int type

Adds the curve bounded by the 4 points passsed as parameters
of type specified in type (see matrix.h for curve type constants)
to the matrix edges
====================*/
void add_curve( struct matrix *edges,
                double x0, double y0,
                double x1, double y1,
                double x2, double y2,
                double x3, double y3,
                int step, int type ) {
  double t, x, y;
  int i;
  struct matrix *xcoefs;
  struct matrix *ycoefs;

  xcoefs = generate_curve_coefs(x0, x1, x2, x3, type);
  ycoefs = generate_curve_coefs(y0, y1, y2, y3, type);

  /* print_matrix(xcoefs); */
  /* printf("\n"); */
  /* print_matrix(ycoefs); */

  for (i=1; i<=step; i++) {
    t = (double)i/step;

    x = xcoefs->m[0][0] *t*t*t + xcoefs->m[1][0] *t*t+
      xcoefs->m[2][0] *t + xcoefs->m[3][0];
    y = ycoefs->m[0][0] *t*t*t + ycoefs->m[1][0] *t*t+
      ycoefs->m[2][0] *t + ycoefs->m[3][0];

    add_edge(edges, x0, y0, 0, x, y, 0);
    x0 = x;
    y0 = y;
  }

  free_matrix(xcoefs);
  free_matrix(ycoefs);
}



/*======== void add_point() ==========
Inputs:   struct matrix * points
         int x
         int y
         int z
Returns:
adds point (x, y, z) to points and increment points.lastcol
if points is full, should call grow on points
====================*/
void add_point( struct matrix * points, double x, double y, double z) {

  if ( points->lastcol == points->cols )
    grow_matrix( points, points->lastcol + 100 );

  points->m[0][ points->lastcol ] = x;
  points->m[1][ points->lastcol ] = y;
  points->m[2][ points->lastcol ] = z;
  points->m[3][ points->lastcol ] = 1;
  points->lastcol++;
} //end add_point

/*======== void add_edge() ==========
Inputs:   struct matrix * points
          int x0, int y0, int z0, int x1, int y1, int z1
Returns:
add the line connecting (x0, y0, z0) to (x1, y1, z1) to points
should use add_point
====================*/
void add_edge( struct matrix * points,
	       double x0, double y0, double z0,
	       double x1, double y1, double z1) {
  add_point( points, x0, y0, z0 );
  add_point( points, x1, y1, z1 );
}

/*======== void draw_lines() ==========
Inputs:   struct matrix * points
         screen s
         color c
Returns:
Go through points 2 at a time and call draw_line to add that line
to the screen
====================*/
void draw_lines( struct matrix * points, screen s, zbuffer zb, color c) {

 if ( points->lastcol < 2 ) {
   printf("Need at least 2 points to draw a line!\n");
   return;
 }

 int point;
 for (point=0; point < points->lastcol-1; point+=2)
   draw_line( points->m[0][point],
              points->m[1][point],
              points->m[2][point],
              points->m[0][point+1],
              points->m[1][point+1],
              points->m[2][point+1],
              s, zb, c);
}// end draw_lines

void draw_line(int x0, int y0, double z0, int x1, int y1, double z1, screen s, zbuffer zb, color c) {
    int temp;
    double temp_z;

    if (x0 > x1) {
        temp = x0;
        x0 = x1;
        x1 = temp;

        temp = y0;
        y0 = y1;
        y1 = temp;

        temp_z = z0;
        z0 = z1;
        z1 = temp_z;
    }

    int x = x0;
    int y = y0;
    double z = z0;
    double dz;

    int a = y1 - y0;
    int b = x0 - x1;

    // octants 1, 5
    if (a > 0 && a <= -b) {
        // printf("(%d, %d) to (%d, %d)\n", x0, y0, x1, y1);
        int d = 2 * a + b;

        a *= 2;
        b *= 2;

        dz = (z1 - z0) / (x1 - x0);

        while (x < x1) {
            plot(s, zb, c, x, y, z);
            if (d > 0) {
                y++;
                d += b;
            }
            x++;
            z += dz;
            d += a;
        }
    }

    // octants 2, 6
    else if (a > -b) {

        int d = a + 2 * b;

        a *= 2;
        b *= 2;

        dz = (z1 - z0) / (y1 - y0);

        while (y < y1) {
          plot(s, zb, c, x, y, z);
            if (d < 0) {
                x++;
                d += a;
            }
            y++;
            z += dz;
            d += b;
        }
    }

    // octants 3, 7
    else if (a < b) {

        int d = a + -2 * b;

        a *= 2;
        b *= 2;

        dz = (z1 - z0) / (x1 - x0);

        while (y > y1) {
          plot(s, zb, c, x, y, z);
            if (d > 0) {

                x++;
                d += a;
            }
            y--;
            z += dz;
            d -= b;
        }
    }

    // octants 4, 8
    else {

        int d = 2 * a - b;

        a *= 2;
        b *= 2;

        dz = (z1 - z0) / (x1 - x0);

        while (x < x1) {
          plot(s, zb, c, x, y, z);
            if (d < 0) {
                y--;
                d -= b;
            }
            x++;
            z += dz;
            d += a;
        }
    }
}
// void draw_line(int x0, int y0, double z0,
//                int x1, int y1, double z1,
//                screen s, zbuffer zb, color c) {
//
//
//   int x, y, d, A, B;
//   int dy_east, dy_northeast, dx_east, dx_northeast, d_east, d_northeast;
//   int loop_start, loop_end;
//
//   //swap points if going right -> left
//   int xt, yt;
//   if (x0 > x1) {
//     xt = x0;
//     yt = y0;
//     x0 = x1;
//     y0 = y1;
//     z0 = z1;
//     x1 = xt;
//     y1 = yt;
//   }
//
//   x = x0;
//   y = y0;
//   A = 2 * (y1 - y0);
//   B = -2 * (x1 - x0);
//   int wide = 0;
//   int tall = 0;
//   //octants 1 and 8
//   if ( abs(x1 - x0) >= abs(y1 - y0) ) { //octant 1/8
//     wide = 1;
//     loop_start = x;
//     loop_end = x1;
//     dx_east = dx_northeast = 1;
//     dy_east = 0;
//     d_east = A;
//     if ( A > 0 ) { //octant 1
//       d = A + B/2;
//       dy_northeast = 1;
//       d_northeast = A + B;
//     }
//     else { //octant 8
//       d = A - B/2;
//       dy_northeast = -1;
//       d_northeast = A - B;
//     }
//   }//end octant 1/8
//   else { //octant 2/7
//     tall = 1;
//     dx_east = 0;
//     dx_northeast = 1;
//     if ( A > 0 ) {     //octant 2
//       d = A/2 + B;
//       dy_east = dy_northeast = 1;
//       d_northeast = A + B;
//       d_east = B;
//       loop_start = y;
//       loop_end = y1;
//     }
//     else {     //octant 7
//       d = A/2 - B;
//       dy_east = dy_northeast = -1;
//       d_northeast = A - B;
//       d_east = -1 * B;
//       loop_start = y1;
//       loop_end = y;
//     }
//   }
//
//   while ( loop_start < loop_end ) {
//
//     plot( s, zb, c, x, y, 0);
//     if ( (wide && ((A > 0 && d > 0) ||
//                    (A < 0 && d < 0)))
//          ||
//          (tall && ((A > 0 && d < 0 ) ||
//                    (A < 0 && d > 0) ))) {
//       y+= dy_northeast;
//       d+= d_northeast;
//       x+= dx_northeast;
//     }
//     else {
//       x+= dx_east;
//       y+= dy_east;
//       d+= d_east;
//     }
//     loop_start++;
//   } //end drawing loop
//   plot( s, zb, c, x1, y1, 0 );
// } //end draw_line
