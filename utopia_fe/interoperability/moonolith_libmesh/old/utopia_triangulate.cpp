/*  triangulate.c

  Triangulate a polygon via ear-clipping, and compute the area
  of a polygon.

  by: Steven Skiena
  begun: May 7, 2002
*/

/*
Copyright 2003 by Steven S. Skiena; all rights reserved.

Permission is granted for use in non-commerical applications
provided this copyright notice remains intact and unchanged.

This program appears in my book:

"Programming Challenges: The Programming Contest Training Manual"
by Steven Skiena and Miguel Revilla, Springer-Verlag, New York 2003.

See our website www.programming-challenges.com for additional information.

This book can be ordered from Amazon.com at

http://www.amazon.com/exec/obidos/ASIN/0387001638/thealgorithmrepo/

*/

#include "utopia_triangulate.hpp"
#include <math.h>
#include <stdio.h>
#include <vector>

#define PI 3.1415926  /* ratio of circumference to diameter */
#define EPSILON 1e-10 /* a quantity small enough to be zero */

typedef struct {
    double a; /* x-coefficient */
    double b; /* y-coefficient */
    double c; /* constant term */
} line;

#define DIMENSION 2 /* dimension of points */
#define X 0         /* x-coordinate index */
#define Y 1         /* y-coordinate index */

typedef double point[DIMENSION];

#define MAXPOLY 200 /* maximum number of points in a polygon */

typedef struct {
    int n;            /* number of points in polygon */
    point p[MAXPOLY]; /* array of points in polygon */
} polygon;

typedef struct {
    point p1, p2; /* endpoints of line segment */
} segment;

typedef point triangle[3]; /* triangle datatype */

typedef struct {
    int n;             /* number of triangles in triangulation */
    int t[MAXPOLY][3]; /* indicies of vertices in triangulation */
} triangulation;

typedef struct {
    point c;  /* center of circle */
    double r; /* radius of circle */
} circle;

double signed_triangle_area(point a, point b, point c) {
    return ((a[X] * b[Y] - a[Y] * b[X] + a[Y] * c[X] - a[X] * c[Y] + b[X] * c[Y] - c[X] * b[Y]) / 2.0);
}

double triangle_area(point a, point b, point c) { return (fabs(signed_triangle_area(a, b, c))); }

bool ccw(point a, point b, point c) { return (signed_triangle_area(a, b, c) > EPSILON); }

bool cw(point a, point b, point c) { return (signed_triangle_area(a, b, c) < -EPSILON); }

void copy_point(const double *src, double *dest) {
    dest[0] = src[0];
    dest[1] = src[1];
}

bool point_in_triangle(point p, triangle t) {
    int i; /* counter */

    for (i = 0; i < 3; i++)
        if (cw(t[i], t[(i + 1) % 3], p)) return (false);

    return (true);
}

bool ear_Q(int i, int j, int k, polygon *p) {
    triangle t; /* coordinates for points i,j,k */
    int m;      /* counter */

    copy_point(p->p[i], t[0]);
    copy_point(p->p[j], t[1]);
    copy_point(p->p[k], t[2]);

    if (cw(t[0], t[1], t[2])) return (false);

    for (m = 0; m < p->n; m++) {
        if ((m != i) && (m != j) && (m != k))
            if (point_in_triangle(p->p[m], t)) return (false);
    }

    return (true);
}

void add_triangle(triangulation *t, int i, int j, int k, polygon *p) {
    int n; /* number of triangles in t */

    n = t->n;

    t->t[n][0] = i;
    t->t[n][1] = j;
    t->t[n][2] = k;

    t->n = n + 1;
}

void triangulate(polygon *p, triangulation *t) {
    int l[MAXPOLY], r[MAXPOLY]; /* left/right neighbor indices */
    int i;                      /* counter */

    for (i = 0; i < p->n; i++) { /* initialization */
        l[i] = ((i - 1) + p->n) % p->n;
        r[i] = ((i + 1) + p->n) % p->n;
    }

    t->n = 0;
    i = p->n - 1;
    while (t->n < (p->n - 2)) {
        i = r[i];
        if (ear_Q(l[i], i, r[i], p)) {
            add_triangle(t, l[i], i, r[i], p);
            l[r[i]] = l[i];
            r[l[i]] = r[i];
        }
    }
}

double area_triangulation(polygon *p) {
    triangulation t;    /* output triangulation */
    double total = 0.0; /* total area so far */
    int i;              /* counter */

    triangulate(p, &t);
    for (i = 0; i < t.n; i++) total += triangle_area(p->p[t.t[i][0]], p->p[t.t[i][1]], p->p[t.t[i][2]]);

    return (total);
}

double area(polygon *p) {
    double total = 0.0; /* total area so far */
    int i, j;           /* counters */

    for (i = 0; i < p->n; i++) {
        j = (i + 1) % p->n;
        total += (p->p[i][X] * p->p[j][Y]) - (p->p[j][X] * p->p[i][Y]);
    }

    return (total / 2.0);
}

// main(){
//   polygon p;      /* input polygon */
//   triangulation t;    /* output triangulation */
//   int i;        /* counter */
//   double area(), area_triangulate();

//   scanf("%d",&p.n);
//   for (i=0; i<p.n; i++)
//     scanf("%lf %lf",&p.p[i][X],&p.p[i][Y]);

// /*
//   print_polygon(&p);
//   triangulate(&p,&t);
//   print_triangulation(&t);
// */

//   printf(" area via triangulation = %f\n", area_triangulation(&p));
//   printf(" area slick = %f\n", area(&p));

// }

void print_triangulation(triangulation *t) {
    int i, j; /* counters */

    for (i = 0; i < t->n; i++) {
        for (j = 0; j < 3; j++) printf(" %d ", t->t[i][j]);
        /*
            for (j=0; j<3; j++)
              printf(" (%5.3f,%5.3f)",t->t[i][j][X],t->t[i][j][Y]);
        */
        printf("\n");
    }
}

void triangulate_polygon(const int n_vertices, const double *in_polygon, std::vector<int> &result) {
    polygon p;
    triangulation t;

    p.n = n_vertices;
    for (int i = 0; i < n_vertices; ++i) {
        p.p[i][X] = in_polygon[i * 2];
        p.p[i][Y] = in_polygon[i * 2 + 1];
    }

    triangulate(&p, &t);
    result.clear();
    result.reserve(t.n * 3);

    for (int i = 0; i < t.n; i++) {
        for (int j = 0; j < 3; j++) {
            result.push_back(t.t[i][j]);
        }
    }
}

void triangulate_polygon(const utopia::Polygon2 &in_polygon, std::vector<int> &result) {
    polygon p;
    triangulation t;

    int n_vertices = in_polygon.size();

    p.n = n_vertices;
    for (int i = 0; i < n_vertices; ++i) {
        p.p[i][X] = in_polygon.points[i].x;
        p.p[i][Y] = in_polygon.points[i].y;
    }

    triangulate(&p, &t);
    result.clear();
    result.reserve(t.n * 3);

    for (int i = 0; i < t.n; i++) {
        for (int j = 0; j < 3; j++) {
            result.push_back(t.t[i][j]);
        }
    }
}