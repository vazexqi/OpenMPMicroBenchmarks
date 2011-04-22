/*$Header: /home/domeika/spec2000/newver3/RCS/scanner.c,v 1.3 1997/08/31 23:18:45 domeika Exp domeika $*/
/*$Header: /home/domeika/spec2000/newver3/RCS/scanner.c,v 1.3 1997/08/31 23:18:45 domeika Exp domeika $*/
/*	Mods by Bodo Parady
 *	Changes omp_get_num_procs to omp_get_max_threads
 *		double c to cc to avoid conflict with integer c
 *       Removed scancount
 *       Removed #if 0, and removed redundant references to omp_get_thread_nu
 *  Modified by Nicholas Chen to eliminate redundant #if 0, convert paralle for to tasks
 */

#ifndef SPEC_HPG_NTOS
#include <unistd.h>
#endif
#include <string.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/timeb.h>

#include <omp.h>

#define DISPVIGILANCE 0               /*Display vigilance passing*/
#define MIN_OUTPUT_VIGILANCE 0.99   /*Minimum vigilance necessary to output the match*/
#define NUM_ALLOWABLE_RETRYS 3        /*Number of time to retry matching during scanning if no match is found*/
#define LIMIT_TRAIN_OBJ 1             /*Limit training to 'numpatterns' objects*/
#define TIMESTATS 1                   /*Display performance stats at the end*/
#define DISPSTARTUP 0                 /*Display basic runtime monitoring debug info*/
#define COLLECTDATA 1                 /*Collect Data into winner, highx, highy, highest_confidence, set_high*/
#define FLOAT_COMPARE_TOLERANCE 0.000001  /*Vigilance for floating point value compares*/

#ifndef INTS_PER_CACHELINE
#define INTS_PER_CACHELINE 4        /* Must be greater than 2*/
#endif
#ifndef DBLS_PER_CACHELINE
#define DBLS_PER_CACHELINE 4        /* Must be greater than 2*/
#endif

struct lnkNode {
  int F2Neuron;
  double Vigilance;
  int x, y;
  struct lnkNode* next;
};

struct lnkNode* head;

int numTraining;

#define TRUE 1
#define FALSE 0

#define MODE_RECOGNIZE 0
#define MODE_TRAIN1 1
#define MODE_TRAIN2 2

void alloc_td_bu();

int multiplier;

unsigned char **cimage;
double **tds;
double **bus;
double **busp;
int numthreads;
int lwidth,lheight;
int width, height,numinputs;
long i,j;

#if TIMESTATS
struct timeb* startTime;
struct timeb* parallelStartTime;
struct timeb* endTime;
#endif

int *pass_flag;
int **highx,**highy;
double **highest_confidence;
int **set_high;

int **winner, **cp;

int o;

int numf1s,numf2s,numpatterns,resonant;

#define DB1 1
#define DB2 0
#define DB3 1

FILE *fp;


double a, b, cc, d, theta, delta_t;
double rho;

typedef struct {
  double *I;
  double W;
  double X;
  double V;
  double U;
  double P;
  double Q;
  double R;
} f1_neuron;

f1_neuron **f1_layer;

typedef struct {
  double y;
  int reset;
} xyz;

xyz **Y;


void add_list_item(int neuron, double vigi, int x, int y) {
  struct lnkNode* ptr;
  struct lnkNode* temp;

  temp = (struct lnkNode*) malloc( sizeof(struct lnkNode) );
  if(temp == NULL) {
    printf("malloc problems\n");
    exit(1);
  }

  temp->x = x;
  temp->y = y;
  temp->F2Neuron = neuron;
  temp->Vigilance = vigi;


  temp->next = head;
  head = temp;
}

void sort_list_items() { /* sorts by the X values */
  struct lnkNode *nHead;
  struct lnkNode *nCurrent;
  struct lnkNode *previousgreat;
  struct lnkNode *previous;
  struct lnkNode *greatest;
  struct lnkNode *current;

  nHead = NULL;
  nCurrent = NULL;
  previousgreat = NULL;
  previous = NULL;
  greatest = head;
  current = head;

  while(head != NULL) {
    previous = NULL;
    previousgreat = NULL;
    current = head;
    greatest = head;

    while(current != NULL) {
      if(current->x > greatest->x) {
        greatest = current;
        previousgreat = previous;
      }
      previous = current;
      current = current->next;
    }

    if(previousgreat != NULL) {
      previousgreat->next = greatest->next;
    } else {
      head = greatest->next;
    }

    greatest->next = NULL;

    if(nHead == NULL) {
      nHead = greatest;
    } else {
      nCurrent->next = greatest;
    }
    nCurrent = greatest;
  }

  head = nHead;
}

double g(int i) {
  double result;

  if(i != winner[o][0]) result = 0;
  else if(Y[o][i].y > 0)
    result = d;
  else
    result = 0;
  return result;
}

void find_match(int o) {
  int i;

  winner[o][0] = 0;
  for(i = 0; i < numf2s; i++)
    if(Y[o][i].y > Y[o][winner[o][0]].y) winner[o][0] = i;
}

double simtest() {
  int j,varNumF1;
  double sum,norm;
  double temp_sum;

  varNumF1 = numf1s;

  sum = norm = 0;
  for(j = 0; j < varNumF1; j++) {
    norm += f1_layer[o][j].P * f1_layer[o][j].P;
  }
  norm = sqrt( (double) norm );
  norm *= cc;
  sum += norm;
  norm = 0;

  for(j = 0; j < varNumF1; j++) {
    temp_sum = f1_layer[o][j].U * f1_layer[o][j].U;
    norm += temp_sum;
  }
  norm = sqrt( (double) norm );
  sum += norm;


  norm = 0;

  for(j = 0; j < varNumF1; j++) {
    f1_layer[o][j].R =  (f1_layer[o][j].U + cc * f1_layer[o][j].P) / sum;
    norm += f1_layer[o][j].R * f1_layer[o][j].R;
  }



  norm = sqrt( (double) norm );
#ifdef DEBUG
  if( DB2 && (winner[o][0] == 0) ) printf("%5.3f",norm);
#endif
  return norm;
}


double simtest2(int o) {
  int j;
  double Su,Sp,numerator,denom;
  double su,sp;
  double su2,sp2;
  double sup;
  double r;
  double e = 0.0000000001;  /* 10e-10 */

  su = sp = sup = su2 = sp2 = numerator = 0.0;
  for(j = 0; j < numf1s; j++) {
    su += f1_layer[o][j].U;
    sp += f1_layer[o][j].P;
    su2 += f1_layer[o][j].U * f1_layer[o][j].U;
    sp2 += f1_layer[o][j].P * f1_layer[o][j].P;
    sup += f1_layer[o][j].U * f1_layer[o][j].P;
  }
  Su = ( (double)numf1s * su2 - su * su ) / ( (double)numf1s * ( (double)numf1s - 1.0 ) );
  Su = sqrt(Su);
  Sp = ( (double)numf1s * sp2 - sp * sp ) / ( (double)numf1s * ( (double)numf1s - 1.0 ) );
  Sp = sqrt(Sp);
  numerator = (double) numf1s * sup - su * sp;
  denom = sqrt( (double) numf1s * su2 - su * su ) * sqrt( (double) numf1s * sp2 - sp * sp );
  r = (numerator + e) / (denom + e);

  if( (numerator == 0) || (denom == 0) ) {
    fprintf(stderr,"potential div by zero");
    r = 1;
  }
  if( (numerator != 0) && (denom == 0) ) {
    fprintf(stderr,"div by zero");
    r = 1;
  }
  r *= r;      /*  trying to match inverse images */
#ifdef DEBUG
  if(DB1) printf("  simtest2(r) = %8.6f\n",r);
  if( DB2 && (winner[o][0] == 0) ) printf("%8.4f",r);
#endif
  return r;
}

void weightadj() {
  int i,j,k;
  double temp;
  double er = 0.000000001;

#ifdef DEBUG
  int bad_count;
  int good_count;
#endif


  i = winner[o][0];
#ifdef DEBUG
  fprintf(stdout,"winner %d\n",i);
#endif

  for(k = 0; k < 1; k++) {
    resonant = 0;
    for(j = 0; j < numf1s; j++) {
      temp = tds[j][i];
      tds[j][i] += g(i) * (f1_layer[o][j].P - tds[j][i]) * delta_t;
      if(fabs(temp - tds[j][i]) <= er) resonant = 1;
    }
#ifdef DEBUG
    bad_count = 0;
    good_count = 0;
#endif
    for(j = 0; j < numf1s; j++) {
      temp = bus[j][i];
      bus[j][i] += g(i) * (f1_layer[o][j].P - bus[j][i]) * delta_t;
      if( (fabs(temp - bus[j][i]) <= er) && resonant ) {
#ifdef DEBUG
        good_count++;
#endif
        resonant = 1;
      }else {
#ifdef DEBUG
        bad_count++;
#endif
        resonant = 0;
      }
    }
  }
#ifdef DEBUG
  printf("bad %d good %d\n",bad_count,good_count);
#endif
}

/* init_globs - initialize ART 2 parameters based
 *   on whether we are training a network or loading
 *   the weights of a pretrained network.
 */
void init_globs(int mode) {
  int i;

  head = NULL;
  pass_flag = (int*) malloc(sizeof(int) * numthreads);
  if(pass_flag == NULL) {
    printf("malloc error\n");
    exit(1);
  }
  winner = (int**) malloc(sizeof(int*) * numthreads);
  if(winner == NULL) {
    printf("malloc error\n");
    exit(1);
  }
  cp = (int**) malloc(sizeof(int*) * numthreads);
  if(cp == NULL) {
    printf("malloc error\n");
    exit(1);
  }

  highx = (int**) malloc(sizeof(int*) * numthreads);
  if(highx == NULL) {
    printf("malloc error\n");
    exit(1);
  }

  highy = (int**) malloc(sizeof(int*) * numthreads);
  if(highy == NULL) {
    printf("malloc error\n");
    exit(1);
  }

  highest_confidence = (double**) malloc(sizeof(double*) * numthreads);
  if(highest_confidence == NULL) {
    printf("malloc error\n");
    exit(1);
  }

  set_high = (int**) malloc(sizeof(int*) * numthreads);
  if(set_high == NULL) {
    printf("malloc error\n");
    exit(1);
  }

  for(i = 0; i < numthreads; i++) {
    #pragma omp task firstprivate(i)
    {
      winner[i] = (int*) malloc( INTS_PER_CACHELINE * sizeof(int) );
      if(winner[i] == NULL) {
        printf("malloc error\n");
        exit(1);
      }
      winner[i][0] = 0;
      cp[i] = (int*) malloc( INTS_PER_CACHELINE * sizeof(int) );
      if(cp[i] == NULL) {
        printf("malloc error\n");
        exit(1);
      }
      cp[i][0] = 0;
      highx[i] = (int*) malloc( INTS_PER_CACHELINE * sizeof(int) );
      if(highx[i] == NULL) {
        printf("malloc error\n");
        exit(1);
      }
      highx[i][0] = 0;
      highy[i] = (int*) malloc( INTS_PER_CACHELINE * sizeof(int) );
      if(highy[i] == NULL) {
        printf("malloc error\n");
        exit(1);
      }
      highy[i][0] = 0;
      highest_confidence[i] = (double*) malloc( DBLS_PER_CACHELINE * sizeof(double) );
      if(highest_confidence == NULL) {
        printf("malloc error\n");
        exit(1);
      }
      highest_confidence[i][0] = 0.0;
      set_high[i] = (int*) malloc( INTS_PER_CACHELINE * sizeof(int) );
      if(set_high == NULL) {
        printf("malloc error\n");
        exit(1);
      }
      set_high[i][0] = 0;
    }
  }



  if(mode == MODE_RECOGNIZE) {
    a       = 255;
    b       = 0.0;
    cc      = 0.11;
    d       = 0.9;
    theta   = 1 / sqrt( (double) numf1s );
    delta_t = 0.1;
    rho     = 0.70;
  }else {
    a       = 255;
    b       = 10.0;
    cc      = 0.11;
    d       = 0.9;
    theta   = 1 / sqrt( (double) numf1s );
    delta_t = 0.7;
    rho     = 0.95;
  }
} /* end of init_globs */


void init_net() {
  int i,j;

  f1_layer = (f1_neuron **)malloc( numthreads * sizeof (f1_neuron*) );
  if(f1_layer == NULL) {
    fprintf(stderr,"malloc error in init_net\n");
    exit(1);
  }
  for(i = 0; i < numthreads; i++) {
    #pragma omp task firstprivate (i)
    {
      f1_layer[i] = (f1_neuron*) malloc( numf1s * sizeof(f1_neuron) );
      if(f1_layer[i] == NULL) {
        fprintf(stderr,"malloc error in init_net\n");
        exit(1);
      }
    }
  }

  for(j = 0; j < numthreads; j++) {
    #pragma omp task firstprivate(j) private(i)
    {
      for(i = 0; i < numf1s; i++) {
        f1_layer[j][i].I = (double *)malloc( 2 * sizeof(double) );
        if(f1_layer[j][i].I == NULL) {
          fprintf(stderr,"malloc error in init_net\n");
          exit(1);
        }
        f1_layer[j][i].W = 0.0;
        f1_layer[j][i].X = 0.0;
        f1_layer[j][i].V = 0.0;
        f1_layer[j][i].U = 0.0;
        f1_layer[j][i].P = 0.0;
        f1_layer[j][i].Q = 0.0;
        f1_layer[j][i].R = 0.0;
      }
    }
  }

  Y = (xyz**)malloc( numthreads * sizeof(xyz*) );
  if(Y == NULL) {
    fprintf(stdout,"Malloc error for Y\n");
    exit(1);
  }

  for(i = 0; i < numthreads; i++) {
  #pragma omp task firstprivate(i)
    {
      Y[i] = (xyz*) malloc( numf2s * sizeof(xyz) );
      if(Y[i] == NULL) {
        fprintf(stderr,"malloc error in init_net\n");
      }
      Y[i][0].y = 0.0;
    }
  }
}


void show_pat() {
  int i;

  for(i = 0; i < numf1s; i++) {
    if( (i % 5) == 0 ) printf("\n");
    printf(" %8.5f ",f1_layer[o][i].I[cp[o][0]]);
  }
  printf("\n\n");
/*   getchar();*/
}

void reset_nodes() {
  int i,varNumF1;
  int o;

#ifdef _OPENMP
  o = omp_get_thread_num();
#else
  o = 0;
#endif

  varNumF1 = numf1s;


  for(i = 0; i < varNumF1; i++) {
    f1_layer[o][i].W = 0.0;
    f1_layer[o][i].X = 0.0;
    f1_layer[o][i].V = 0.0;
    f1_layer[o][i].U = 0.0;
    f1_layer[o][i].P = 0.0;
    f1_layer[o][i].Q = 0.0;
    f1_layer[o][i].R = 0.0;
  }
  for(i = 0; i < numf2s; i++) {
    Y[o][i].y = 0.0;
    Y[o][i].reset = 0;
  }
  winner[o][0] = 0;
  resonant = 0;
}


void reset_nodes2() {
  int i,varNumF1;
  int o;

#ifdef _OPENMP
  o = omp_get_thread_num();
#else
  o = 0;
#endif
  varNumF1 = numf1s;


  for(i = 0; i < varNumF1; i++) {
    f1_layer[o][i].W = 0.0;
    f1_layer[o][i].X = 0.0;
    f1_layer[o][i].V = 0.0;
    f1_layer[o][i].U = 0.0;
    f1_layer[o][i].P = 0.0;
    f1_layer[o][i].Q = 0.0;
    f1_layer[o][i].R = 0.0;
  }
  for(i = 0; i < numf2s; i++)
    Y[o][i].y = 0.0;
  winner[o][0] = 0;
  resonant = 0;
}

void print_weights() {
  int i,j;

  /* print td's */
  printf("============  TOP down WEIGHTS ==============\n");
  for(i = 0; i < numf1s; i++)
    for(j = 0; j < numf2s; j++)
      if( j == (numf2s - 1) ) printf(" %8.16f\n",tds[i][j]);
      else
        printf(" %8.16f ",tds[i][j]);
  /* print bu's */
  printf("============  BOTTOM up WEIGHTS ==============\n");
  for(i = 0; i < numf1s; i++)
    for(j = 0; j < numf2s; j++)
      if( j == (numf2s - 1) ) printf(" %8.16f\n",bus[i][j]);
      else
        printf(" %8.16f ",bus[i][j]);
}


void print_f12() {
  int j;
  int o;

#ifdef _OPENMP
  o = omp_get_thread_num();
#else
  o = 0;
#endif
  printf("\n\n");
  for(j = 0; j < numf2s; j += 10)
    printf(" j = %i  Y= %9.7f\n",j,Y[o][j].y);
}

void reset_unnecessary_f2s() {
  int i;
  int k;
  int o;

#ifdef _OPENMP
  o = omp_get_thread_num();
#else
  o = 0;
#endif
  for(i = numpatterns; i < numf2s; i++)
    Y[o][i].reset = 1;
}

void compute_train_match(int o, int *f1res, int *spot) {
  int varNumF1;
  int varNumF2;
  int varNumCp;
  double varNumA;
  double varNumB;
  double varNumD;
  double varNumTheta;
  double oldTnorm;

  int ti,tj,tresult;
  double tnorm,xr,qr,tsum,ttemp;

  varNumF1 = numf1s;
  varNumF2 = numf2s;
  varNumCp = cp[o][0];
  varNumA = a;
  varNumB = b;
  varNumD = d;
  varNumTheta = theta;

  /* Compute F1 layer - W values */
  tnorm = 0;

/*Q03 */
  for(ti = 0; ti < varNumF1; ti++) {
    f1_layer[o][ti].W = f1_layer[o][ti].I[varNumCp] + varNumA * (f1_layer[o][ti].U);
    tnorm += f1_layer[o][ti].W * f1_layer[o][ti].W;
  }
  tnorm =  sqrt( (double)tnorm );

  /* Compute F1 layer - V values */
  oldTnorm = tnorm;
  tnorm = 0;


  for(ti = 0; ti < varNumF1; ti++) {
    f1_layer[o][ti].X = f1_layer[o][ti].W / oldTnorm;

    if(f1_layer[o][ti].X < varNumTheta) xr = 0;
    else
      xr = f1_layer[o][ti].X;
    if(f1_layer[o][ti].Q < varNumTheta) qr = 0;
    else
      qr = f1_layer[o][ti].Q;
    f1_layer[o][ti].V = xr + varNumB * qr;
    tnorm += f1_layer[o][ti].V * f1_layer[o][ti].V;
  }

/* Compute F1 layer - U values */
  tnorm = sqrt( (double) tnorm );

/* Compute F1 layer - P values */
  oldTnorm = tnorm;
  tnorm = 0;
  tsum = 0;
  tresult = 1;


  for(ti = 0; ti < varNumF1; ti++) {
    f1_layer[o][ti].U = f1_layer[o][ti].V / oldTnorm;

    tsum = 0;
    ttemp = f1_layer[o][ti].P;
#if LIMIT_TRAIN_OBJ
    for(tj = *spot; tj < numpatterns; tj++)
#else
    for(tj = *spot; tj < varNumF2; tj++)
#endif
    {
      if( (tj == winner[o][0]) && (Y[o][tj].y > 0) ) tsum += tds[ti][tj] * varNumD;
    }

    f1_layer[o][ti].P = f1_layer[o][ti].U + tsum;

    tnorm += f1_layer[o][ti].P * f1_layer[o][ti].P;

    if(fabs(ttemp - f1_layer[o][ti].P) < FLOAT_COMPARE_TOLERANCE) tresult = 0;
  }


  *f1res = tresult;

/* Compute F1 - Q values */

  tnorm = sqrt( (double) tnorm );



  for(tj = 0; tj < varNumF1; tj++)
    f1_layer[o][tj].Q = f1_layer[o][tj].P;

/* Compute F2 - y values */

#if LIMIT_TRAIN_OBJ
  for(tj = *spot; tj < numpatterns; tj++)
#else
  for(tj = *spot; tj < varNumF2; tj++)
#endif
  {
    Y[o][tj].y = 0;
    if(!Y[o][tj].reset)
      for(ti = 0; ti < varNumF1; ti++)
        Y[o][tj].y += f1_layer[o][ti].P * bus[ti][tj];
  }

/* Find match */
  winner[o][0] = 0;
#if LIMIT_TRAIN_OBJ
  for(ti = *spot; ti < numpatterns; ti++)
#else
  for(ti = *spot; ti < varNumF2; ti++)
#endif
  {
    if(Y[o][ti].y > Y[o][winner[o][0]].y) winner[o][0] = ti;
  }
} /* end of compute_train_match()*/


void compute_values_match(int o, int *f1res, int *spot, double **busp) {
  int varNumF1;
  int varNumF2;
  int varNumCp;
  double varNumA;
  double varNumB;
  double varNumD;
  double varNumTheta;
  double oldTnorm;

  int ti,tj,tresult;
  double tnorm,xr,qr,tsum,ttemp;

  varNumF1 = numf1s;
  varNumF2 = numf2s;
  varNumCp = cp[o][0];
  varNumA = a;
  varNumB = b;
  varNumD = d;
  varNumTheta = theta;

/* Compute F1 layer - W values */
  tnorm = 0;

/*Q03*/
  for(ti = 0; ti < varNumF1; ti++) {
    f1_layer[o][ti].W = f1_layer[o][ti].I[varNumCp] + varNumA * (f1_layer[o][ti].U);
    tnorm += f1_layer[o][ti].W * f1_layer[o][ti].W;
  }
  tnorm =  sqrt( (double)tnorm );

/* Compute F1 layer - V values */
  oldTnorm = tnorm;
  tnorm = 0;


  for(ti = 0; ti < varNumF1; ti++) {
    f1_layer[o][ti].X = f1_layer[o][ti].W / oldTnorm;

    if(f1_layer[o][ti].X < varNumTheta) xr = 0;
    else
      xr = f1_layer[o][ti].X;
    if(f1_layer[o][ti].Q < varNumTheta) qr = 0;
    else
      qr = f1_layer[o][ti].Q;
    f1_layer[o][ti].V = xr + varNumB * qr;
    tnorm += f1_layer[o][ti].V * f1_layer[o][ti].V;
  }

/* Compute F1 layer - U values */
  tnorm = sqrt( (double) tnorm );

/* Compute F1 layer - P values */
  oldTnorm = tnorm;
  tnorm = 0;
  tsum = 0;
  tresult = 1;


  for(ti = 0; ti < varNumF1; ti++) {
    f1_layer[o][ti].U = f1_layer[o][ti].V / oldTnorm;

    tsum = 0;
    ttemp = f1_layer[o][ti].P;

    for(tj = *spot; tj < varNumF2; tj++) {
      if( (tj == winner[o][0]) && (Y[o][tj].y > 0) ) tsum += tds[ti][tj] * varNumD;
    }

    f1_layer[o][ti].P = f1_layer[o][ti].U + tsum;

    tnorm += f1_layer[o][ti].P * f1_layer[o][ti].P;

    if(fabs(ttemp - f1_layer[o][ti].P) < FLOAT_COMPARE_TOLERANCE) tresult = 0;
  }


  *f1res = tresult;

/* Compute F1 - Q values */

  tnorm = sqrt( (double) tnorm );



  for(tj = 0; tj < varNumF1; tj++)
    f1_layer[o][tj].Q = f1_layer[o][tj].P;

/* Compute F2 - y values */


  for(tj = *spot; tj < varNumF2; tj++) {
    Y[o][tj].y = 0;
  }
  for(ti = 0; ti < varNumF1; ti++) {
    for(tj = *spot; tj < varNumF2; tj++)
      if(!Y[o][tj].reset) Y[o][tj].y += f1_layer[o][ti].P * busp[ti][tj];
  }

/* Find match */
  winner[o][0] = 0;

  for(ti = *spot; ti < varNumF2; ti++) {
    if(Y[o][ti].y > Y[o][winner[o][0]].y) winner[o][0] = ti;
  }
} /* end of compute_values_,match() */

void train_match(int spot) {
  int j,matched,f1res,mt;
  char matchtest;
  double match_confidence;
  int o;

#ifdef _OPENMP
  o = omp_get_thread_num();
#else
  o = 0;
#endif
  numTraining = numTraining + 1;
  f1res = 0;
  reset_nodes();
  cp[o][0] = spot;
  matched = 0;

  reset_unnecessary_f2s();
  while(!matched) {
    f1res = 0;
    for(j = 0; j < 9 && !f1res ; j++) {
      compute_train_match(o, &f1res, &spot);
    }
#ifdef DEBUG
    if(DB1) print_f12();
    if(DB1) printf("\n num iterations for p to stabalize = %i \n",j);
#endif
    match_confidence = simtest();
#ifdef DEBUG
    fprintf(stdout,"rho %e\n",match_confidence);
#endif
    if( (match_confidence) > rho ) {
#ifdef DEBUG
      if( DB2 && (winner[o][0] == 0) ) printf("#%i",winner[o][0]);
#endif
      weightadj();
      matched = 1;
    }else {Y[o][winner[o][0]].y = 0;
           Y[o][winner[o][0]].reset = 1;
#ifdef DEBUG
           if(DB1) printf("#%iN",winner[o][0]);
#endif
           matchtest = 0;
           for(mt = spot; mt < numf2s; mt++)
             if(Y[o][mt].reset == 0) matchtest = 1;
           if(matchtest) find_match(o);
           else
             matched = 1;}
  } /* end while */
} /* end of train_match() */

int match(int o, int xcoor, int ycoor, double *mcp, double **busp) {
  int j,matched,f1res,mt;
  long c;
  char matchtest;
  double match_confidence;
  int spot;
  int ret;
  int mprint;

  c = 0;
  ret = 0;
  f1res = 0;
  spot = 0;


  cp[o][0] = 0;
  reset_nodes();

  matched = 0;
  mprint = 1;
  while(!matched) {
    if(c++ > NUM_ALLOWABLE_RETRYS) break;
    reset_nodes2();
    f1res = 0;
    for(j = 0; j < 9 && !f1res ; j++) {
      compute_values_match(o, &f1res, &spot, busp);
    }
#ifdef DEBUG
    if(DB1) print_f12();
    if(DB1) printf("\n num iterations for p to stabilize = %i \n",j);
#endif
    match_confidence = simtest2(o);
    if(mprint == 1) {
      *mcp = match_confidence;
      mprint = 0;
    }

    if( (match_confidence) > rho ) {
      /* If the winner[o][0] is not the default F2 neuron (the highest one)
       * we have a match.
       */

      if(winner[o][0] != numf2s - 1) {
        pass_flag[o] = 1;
        ret = 1;
#if DISPVIGILANCE
        fprintf(stdout,"\nF2 neuron %d passes vigilance with a value of %0.4f",winner[o][0],match_confidence);
#if 0
        print_f12();
#endif
#endif
        if(match_confidence >= MIN_OUTPUT_VIGILANCE) {
          add_list_item(winner[o][0], match_confidence,xcoor,ycoor);
#if DISPSTARTUP
          printf("\nNeuron# %d -->Added winner: x %d, y %d, Vigi %f",winner[o][0],xcoor,ycoor,match_confidence);
#endif
        }
        if(match_confidence > highest_confidence[o][winner[o][0]]) {
          highest_confidence[o][winner[o][0]] = match_confidence;
          set_high[o][winner[o][0]] = TRUE;
        }
      }

      matched = 1;
    }else {Y[o][winner[o][0]].y = 0;
           Y[o][winner[o][0]].reset = 1;
#ifdef DEBUG
           if(DB1) printf("#%i No",winner[o][0]);
#endif
           matchtest = 0;
           for(mt = 0; mt < numf2s; mt++)
             if(Y[o][mt].reset == 0) matchtest = 1;
           if(matchtest) find_match(o);
           else
             matched = 1;}
  } /* end while */

  return ret;
}


/*
 *  loadimage - load image to scan
 *  This was rewritten because Windows NT seems to have
 *  problems with sequential calls to getc, scanf, fread,
 *  and fread.  The bug is flaky.  It appears that I'd
 *  only get one big read and all of the rest of the reads
 *  would contain bogus data, yet no error condition is
 *  generated by the read.  Solution: one big read of the
 *  whole image and then a copy to where I really want it,
 *  cimage.
 */
void loadimage(char *input_file) {
  int i,j,r,c;
  int fd;
  int cimgwidth, cimgheight;
  char buffer[64];
  char *superbuffer;

  if( ( fd = open(input_file,O_RDONLY) ) == -1 ) {
    fprintf(stderr,"Error opening %s\n",input_file);
    exit(1);
  }
#ifdef DEBUG
  printf("made it to loadimage\n");
#endif

  /* Strip Format descriptor */
  read(fd,buffer,8);
  /* Read width */
  read(fd,buffer,4);
  for(i = 0; i < 4; i++)
    if(buffer[i] != ' ') width = width * 10 + buffer[i] - '0';

  /* Read height */
  read(fd,buffer,4);
  for(i = 0; i < 4; i++)
    if(buffer[i] != ' ') height = height * 10 + buffer[i] - '0';

  /* calculate cimage size based on multiplier*/
  cimgwidth = width * multiplier;
  cimgheight = height * multiplier;

#ifdef DEBUG
  fprintf(stderr,"width %d, height %d\n",width,height);
#endif

  superbuffer = (char *)malloc( width * height * sizeof(char) );
  if(superbuffer == NULL) {
    fprintf(stderr,"Problems with malloc in loadimage()\n");
    exit(1);
  }
  cimage = (unsigned char **) malloc(sizeof(unsigned char *) * cimgheight);
  if(cimage == NULL) {
    fprintf(stderr,"Problems with malloc in loadimage()\n");
    exit(1);
  }

  for(i = 0; i < cimgheight; i++) {
    cimage[i] = (unsigned char *) malloc( cimgwidth * sizeof(unsigned char) );
    if(cimage[i] == NULL) {
      fprintf(stderr,"Problems with malloc in loadimage()\n");
      exit(1);
    }
  }
#if DISPSTARTUP
  printf("\nRight before reading the image into the big one");
#endif
  read(fd,superbuffer,width * height);
  for(i = 0; i < cimgheight; i += height) {
    for(j = 0; j < cimgwidth; j += width) {
      for(r = 0; r < height; r++) {
        for(c = 0; c < width; c++) {
          cimage[i + r][j + c] = superbuffer[r * width + c];
        }
      }
    }
  }


#ifdef DEBUG
  printf(" \nupper left 10X10 corner of image\n");
  for(i = 0; i < 20; i++) {
    for(j = 0; j < 20; j++)
      printf("%4d",cimage[i][j]);
    printf("\n");
  }
#endif
  width = cimgwidth;
  height = cimgheight;
#if DISPSTARTUP
  printf("\nwidth:  %d  height:  %d  ",width,height);
  printf("\nended loadimage(scanfile) RUNNING BENCHMARK\n");
#endif
}

/* load_weights - load neural net weights
 * This seems to function properly which is odd because
 * loadimage does not.  The only difference really is that
 * loadimage is trying to load bytes of any value, whereas
 * load_weights is looking at a file that only contains
 * ascii characters and reading them as ints or doubles.
 */
void load_weights(char *weightfile) {
  double a;
  long i,j;

  FILE *inp;

  if( ( inp = fopen(weightfile,"r") ) == NULL ) {
    fprintf(stderr,"Unable to open %s\n",weightfile);
    exit(1);
  }

  printf("made it to load_weights\n");
  fscanf(inp,"%d %d",&lwidth,&lheight);
  numf1s = numinputs = lwidth * lheight;
  numf2s = numpatterns + 1;

  alloc_td_bu();

  j = 0;
  for(i = 0; i < numf1s; i++) {
    fscanf(inp,"%le",&a);
    bus[i][j] = tds[i][j] = a;
  }
} /* end of load_weights */

/* alloc_td_bu - memory alloc of top down and bottom up
 *   connections
 */
void alloc_td_bu() {
  bus = (double **)malloc( numf1s * sizeof(double *) );
  tds = (double **)malloc( numf1s * sizeof(double *) );
  if( (bus == NULL) || (tds == NULL) ) {
    fprintf(stderr,"Malloc problem in load_weights\n");
    exit(1);
  }
  for(i = 0; i < numf1s; i++) {
    bus[i] = (double *)malloc( numf2s * sizeof(double) );
    tds[i] = (double *)malloc( numf2s * sizeof(double) );
    if( (bus[i] == NULL) || (tds[i] == NULL) ) {
      fprintf(stderr,"Malloc problem in load_weights, i=%d\n",i);
      exit(1);
    }
  }
} /* end of alloc_td_bu */

/* init_td - initialize top down weights
 *   start signifies which F2 neuron to initialize for.  Enables
 *   training on more than one image.
 */
void init_td(int start) {
  int i,j;

  for(i = 0; i < numf1s; i++)
    for(j = start; j < numf2s; j++)
      tds[i][j] = 0.0;
} /* end of init_td */

/* init_bu - initialize bottom up weights
 */
void init_bu(int start) {
  int i,j;

  for(i = 0; i < numf1s; i++)
    for(j = start; j < numf2s; j++)
      bus[i][j] = 1 / (1.0 - d) / sqrt( (double)numf1s );
} /* end of init_bu */

/* load_train - load a training file into
 *   location f1_layer[o][].I[spot]
 */
void load_train(char *trainfile,int mode, int objects) {
  int i;
  int fd;
  char buffer[64];
  char *superbuffer;
  unsigned char t;
  int spot;

  if(mode == MODE_TRAIN1) {
    spot = 0;
  }else {
    spot = 1;
  }

  if( ( fd = open(trainfile,O_RDONLY) ) == -1 ) {
    fprintf(stderr,"Error opening %s\n",trainfile);
    exit(1);
  }
#ifdef DEBUG
  printf("made it to load_train. opening %s\n",trainfile);
#endif

  lwidth = 0;
  lheight = 0;

  /* Strip Format descriptor */
  read(fd,buffer,8);
  /* Read width */
  read(fd,buffer,4);
  for(i = 0; i < 4; i++)
    if(buffer[i] != ' ') lwidth = lwidth * 10 + buffer[i] - '0';
  /* Read height */
  read(fd,buffer,4);
  for(i = 0; i < 4; i++)
    if(buffer[i] != ' ') lheight = lheight * 10 + buffer[i] - '0';

#ifdef DEBUG
  fprintf(stderr,"width %d, height %d\n",lwidth,lheight);
#endif

  /* The first time through we set up the network
   * based on what is read from the file.
   * The second time through (if we have more than
   * one training file, we make sure the parameters
   * match what was read the first time, e.g. the
   * f1 layer is the same size.
   */
  if(mode == MODE_TRAIN1) {
    numf1s = numinputs = lwidth * lheight;
    numf2s = objects + 1;
    init_globs(MODE_TRAIN1);
    init_net();
  }else {
    if( (lwidth * lheight) != numf1s ) {
      fprintf(stderr,"Dimensions of first image do not match");
      fprintf(stderr," dimensions of second.\n");
      exit(1);
    }
  }

  superbuffer = (char *)malloc( lwidth * lheight * sizeof(char) );
  if(superbuffer == NULL) {
    fprintf(stderr,"Problems with malloc in loadimage()\n");
    exit(1);
  }

  read(fd,superbuffer,lwidth * lheight);
  for(i = 0; i < lheight * lwidth; i++) {
    t = superbuffer[i];
    f1_layer[o][i].I[spot] = (double) t;
  }

  free(superbuffer);
} /* end of load_train */

/* This routine is used to simulate training of other objects.
 * Training on multiple objects would consume the entire execution
 * time of the benchmark.  Instead we train on a few then we simulate
 * others by copying the interconnections for the objects we are trained
 * on and then adding noise to the connections.  This simulates training on
 * other objects.  Need to blur the objects enough to overcome ART's
 * noise filtering.
 */
void sim_other_objects(int low, int high, int stop) {
  int i,j,varNumF1,varHigh;
  int noise1;
  double noise2;

#ifdef DEBUG
  printf("sim other low %d high %d stop %d\n",low,high,stop);
  printf("sim other numf2s %d numpat %d\n",numf2s,numpatterns);
#endif
  if(high <= low) {
    return;
  }
  srand(10);

  varNumF1 = numf1s;
  varHigh = high;

  for(i = low; i < varHigh; i++) {
    for(j = 0; j < varNumF1; j++) {
      if(i % low) {
        tds[j][i] = tds[j][0];
        tds[j][i] = bus[j][0];
      } else {
        tds[j][i] = tds[j][1];
        tds[j][i] = bus[j][1];
      }
    }
  }

  for(i = low; i < varHigh; i++) {
    for(j = 0; j < varNumF1; j++) {
      noise1 = rand() & 0xffff;
      noise2 = (double)noise1 / (double)0xffff;
      tds[j][i] -= noise2;
      bus[j][i] -= noise2;
    }
  }
} /* sim_other_objects */

void setup_base_pattern(int spot) {
  int i,j;

  for(i = 0; i < numf1s; i++) {
    for(j = spot; j < numf2s; j++) {
      tds[i][j] = bus[i][j] = 1.0 / sqrt( (double)numf1s ) / (1 - d);
    }
  }
}

void scan_recognize(int startx, int starty, int endx, int endy, int stride) {
  int gStartY,gEndY,gStartX,gEndX,gLheight,gLwidth,gStride,gPassFlag,avgy;
  int i,j,m,n, r;
  long k;
  int ij, ijmx, jnum, inum;
  long elapsedTimeSec;
  int elapsedTimeMilli;
  int modval;
  int addsubval_out, addsubval_in;
  struct _timeb* startTime;
  struct _timeb* endTime;
  double* mat_con;

  gStartY = starty;
  gStartX = startx;
  gEndY = endy;
  gEndX = endx;
  gLheight = lheight;
  gLwidth = lwidth;
  gStride = stride;


  if( ( starty > (height - lheight + 1) ) || ( startx > (width - lwidth + 1) ) ) {
    fprintf(stderr,"Startx %d or Starty %d is out of range\n", startx, starty);
    exit(1);
  }
  if( ( endy > (height - lheight + 1) ) || ( endx > (width - lwidth + 1) ) ) {
    fprintf(stderr,"endx %d or endy %d is out of range\n", endx, endy);
    exit(1);
  }
#ifdef DEBUG
  if(DB3) {
    fprintf(stdout,"made it to scan_recognize\n");
    fprintf(stdout,"w= %d h = %d lw = %d lh = %d\n",width,height,lwidth,lheight);
  }
#endif

  avgy = abs( (gEndY - gStartY) / numthreads );
  modval = avgy % 2;
  if(modval != 0) {
    addsubval_out = avgy + 1;
    addsubval_in = avgy - 1;
  } else {
    addsubval_out = avgy;
    addsubval_in = avgy;
  }
#ifdef DEBUG
  printf("\naddsubval_out: %d",addsubval_out);
  printf("\naddsubval_in: %d", addsubval_in);
#endif


#if DISPSTARTUP
  printf("\nStartX = %d, EndX = %d, StartY = %d, EndY = %d",gStartX,gEndX,gStartY,gEndY);
  printf("\nThis is avgy=%d",avgy);
#endif

#if TIMESTATS
  ftime(parallelStartTime);
#endif

/*  for (r=gStartY; r < gEndY; r+= addsubval_out)
   {   */
  inum  = (gEndX - gStartX) / gStride;
  jnum  = (gEndY - gStartY) / gStride;
  ijmx = (gEndX - gStartX) * (gEndY - gStartY) / (gStride * gStride);
  mat_con = (double *) malloc( ijmx * sizeof (double) );
  if( (mat_con == NULL) ) {
    fprintf(stderr,"Malloc problem with mat_con\n");
    exit(1);
  }


#pragma omp parallel private (busp, o, i, j)
  {
/*   create a local copy of the bus array   */


    busp = (double **)malloc( numf1s * sizeof(double *) );
    if( (busp == NULL) ) {
      fprintf(stderr,"Malloc problem in private load_weights\n");
      exit(1);
    }
    busp[0] = (double *)malloc( numf1s * numf2s * sizeof(double) );
    if( (busp[0] == NULL) ) {
      fprintf(stderr,"Malloc problem in private load_weights 0\n");
      exit(1);
    }

    for(i = 1; i < numf1s; i++) {
      busp[i] = busp[i - 1] + numf2s;
    }

    for(i = 0; i < numf1s; i++) {
      for(j = 0; j < numf2s; j++) {
        busp[i][j] = bus [i][j];
      }
    }
#pragma omp single
    {
      for(ij = 0; ij < ijmx; ij++) {
        #pragma omp task private(k, m, n, gPassFlag) firstprivate(ij)
        {
          o = omp_get_thread_num();
          j = ( (ij / inum) * gStride ) + gStartY;
          i = ( (ij % inum) * gStride ) + gStartX;
          k = 0;
          for(m = j; m < (gLheight + j); m++)
            for(n = i; n < (gLwidth + i); n++)
              f1_layer[o][k++].I[0] = cimage[m][n];


          gPassFlag = 0;
          gPassFlag = match(o,i,j, &mat_con[ij], busp);

          if(gPassFlag == 1) {
                                #ifdef DEBUG
            printf(" at X= %d Y = %d\n",i,j);
                                #endif
            if(set_high[o][0] == TRUE) {
              highx[o][0] = i;
              highy[o][0] = j;
              set_high[o][0] = FALSE;
            }
            if(set_high[o][1] == TRUE) {
              highx[o][1] = i;
              highy[o][1] = j;
              set_high[o][1] = FALSE;
            }
          }

        #ifdef DEBUG
          else if(DB3)
            printf("0.00#%dx%da%2.1fb%2.1f\n",i,j,a,b);
        #endif
        }
      }
    }
  }
  for(ij = 0; ij < ijmx; ij++) {
    j = ( (ij / inum) * gStride ) + gStartY;
    i = ( (ij % inum) * gStride ) + gStartX;
    printf("x = %i, y = %i, match confidence = %9.7f \n", i, j, mat_con[ij]);
  }
}


int main(int argc, char *argv[]) {
  int k;
  int startx, starty;
  int endx, endy;
  int stride;
  int objects;
  int arg_index;
  int max_threads;
  char *scanfile = NULL;
  char *weightfile = NULL;
  char *trainfile1 = NULL;
  char *trainfile2 = NULL;
  struct lnkNode* ptr;

#if TIMESTATS
  long SerialSec,SerialMilli,ParallelSec,ParallelMilli,TotalSec,TotalMilli;
  double percentSerial,percentParallel;


  startTime = (struct timeb*) malloc( sizeof(struct timeb) );
  parallelStartTime = (struct timeb*) malloc( sizeof(struct timeb) );
  endTime = (struct timeb*) malloc( sizeof(struct timeb) );


  if(startTime == NULL || parallelStartTime == NULL || endTime == NULL) {
    printf("malloc error\n");
    exit(1);
  }
#endif

#ifdef _OPENMP
  max_threads = omp_get_max_threads();
#if DISPSTARTUP
  printf("\n*** RUNNING OPENMP VERSION OF ART BENCHMARK ***\n");
  printf("num procs= %d",max_threads);
#endif

#ifdef NO_OPENMP
  numthreads = 1;
#else
  numthreads = max_threads;
#endif

#if DISPSTARTUP
  printf("\nnum theads= %d\n",numthreads);
#endif
  omp_set_num_threads(numthreads);
#else
  numthreads = 1;
  printf("\n");
#endif
  numTraining = 0;

#if TIMESTATS
  ftime(startTime);
#endif


  if(argc < 2) {
    goto Usage;
  }
  if(argc == 2) {
    if(strcmp(argv[1],"-v") == 0) goto Version;
    else if(strcmp(argv[1],"-h") == 0)
      goto Usage;
  }

  stride = 0;
  startx = 0;
  starty = 0;
  endy = 0;
  endx = 0;
  objects = 0;
  multiplier = 1;
  /* Read command line options */
  arg_index = 1;
#if DISPSTARTUP
  printf("before reading arguments...\n");
#endif
  while(arg_index < argc - 1) {
    if(strcmp(argv[arg_index],"-scanfile") == 0) {
      scanfile = argv[arg_index + 1];
    }else if(strcmp(argv[arg_index],"-weightfile") == 0) {
      weightfile = argv[arg_index + 1];
    }else if(strcmp(argv[arg_index],"-trainfile1") == 0) {
      trainfile1 = argv[arg_index + 1];
    }else if(strcmp(argv[arg_index],"-trainfile2") == 0) {
      trainfile2 = argv[arg_index + 1];
    }else if(strcmp(argv[arg_index],"-startx") == 0) {
      startx = atoi(argv[arg_index + 1]);
    }else if(strcmp(argv[arg_index],"-starty") == 0) {
      starty = atoi(argv[arg_index + 1]);
    }else if(strcmp(argv[arg_index],"-endx") == 0) {
      endx = atoi(argv[arg_index + 1]);
    }else if(strcmp(argv[arg_index],"-endy") == 0) {
      endy = atoi(argv[arg_index + 1]);
    }else if(strcmp(argv[arg_index],"-stride") == 0) {
      stride = atoi(argv[arg_index + 1]);
    }else if(strcmp(argv[arg_index],"-objects") == 0) {
      objects = atoi(argv[arg_index + 1]);
    }else if(strcmp(argv[arg_index],"-tile") == 0) {
      multiplier = atoi(argv[arg_index + 1]);
    } else{
      fprintf(stderr,"ERROR: Unknown option -> %s\n",argv[arg_index]);
      goto Usage;
    }
    arg_index += 2; /* this works as long as options are duals!!! */
  }

  /* Some basic error checking. */
#if DISPSTARTUP
  printf("read arguments, checking...\n");
#endif

  if(scanfile == NULL) {
    fprintf(stderr,"ERROR: Must specify input files\n");
    goto Usage;
  }
  if( (weightfile == NULL) && (trainfile1 == NULL) ) {
    fprintf(stderr,"ERROR: Must specify weightfile or trainfile1\n");
    goto Usage;
  }
  if( (weightfile != NULL) && (trainfile1 != NULL) ) {
    fprintf(stderr,"ERROR: Cannot specify weightfile and trainfile1\n");
    goto Usage;
  }

#ifdef DEBUG
  fprintf(stdout,
          "scanfile = %s\n weightfile = %s\n startx = %d\n starty = %d\n stride = %d\n",
          scanfile,
          weightfile,
          startx,
          starty,
          stride);
#endif
#if DISPSTARTUP
  printf("calling loadimage...\n");
#endif

  printf("Running ART-II Object Recognition...\n");
  loadimage(scanfile);

  /* Differentiate between loading pretrained network (load_weights)
   * and training.  Currently, loading a pretrained network
   * supports only 1 object.  If we train, we can learn to
   * recognize two objects.
   */
  if(weightfile != NULL) {
    numpatterns = 1;
    if(objects == 0) {
      objects = numpatterns;
    }
    load_weights(weightfile);
    init_globs(MODE_RECOGNIZE);
    init_net();
  }else {
    if(trainfile2 != NULL) {
      numpatterns = 2;
      if(objects < numpatterns) {
        objects = numpatterns;
      }
      load_train(trainfile1,MODE_TRAIN1,objects);
      alloc_td_bu();
      init_td(0);
      init_bu(0);
      resonant = k = 0;
      while(!resonant) {
#ifdef DEBUG
        fprintf(stdout,"k=%d\n",k);
#endif
        train_match(0);
        k++;
      }
      load_train(trainfile2,MODE_TRAIN2,objects);
      init_globs(MODE_TRAIN2);
      init_td(1);
      init_bu(1);
      resonant = k = 0;
      while(!resonant) {
#ifdef DEBUG
        fprintf(stdout,"k=%d\n",k);
#endif
        train_match(1);
        k++;
      }
      init_globs(MODE_RECOGNIZE);
      init_td(objects);
      init_bu(objects);
      sim_other_objects(numpatterns,objects,numf2s);
      setup_base_pattern(objects);
    }else {
      numpatterns = 1;
      if(objects < numpatterns) {
        objects = numpatterns;
      }
      load_train(trainfile1,MODE_TRAIN1,objects);

      alloc_td_bu();
      init_td(0);
      init_bu(0);
      resonant = k = 0;
      while(!resonant) {
#ifdef DEBUG
        fprintf(stdout,"k=%d\n",k);
#endif
        train_match(0);
        k++;
      }
      init_globs(MODE_RECOGNIZE);
      init_td(1);
      init_bu(1);
      setup_base_pattern(1);
    }
  }

  /* Set endx and endy if user never specified */
  if(endy == 0) {
    endy = height - lheight;
  }
  if(endx == 0) {
    endx = width - lwidth;
  }
  for(i = 0; i < numthreads; i++) {
    highest_confidence[i][0] = 0.0;
    highest_confidence[i][1] = 0.0;
    highx[i][0] = 0;
    highx[i][1] = 0;
    highy[i][0] = 0;
    highy[i][1] = 0;
    set_high[i][0] = FALSE;
    set_high[i][1] = FALSE;
  }

  printf("\nWorking with %d F2 objects...\n",numf2s);
  scan_recognize(startx, starty, endx, endy, stride);
#if TIMESTATS
  ftime(endTime);
#endif


#if _DEBUG
  printf("\n\nResults:\n");
  for(i = 0; i < numthreads; i++) {
    fprintf(stdout,"THREAD# %d ->Highest vigilance for 1 = %0.4f for object at X = %d, Y = %d\n",i,
            highest_confidence[i][0], highx[i][0], highy[i][0]);
    if(numpatterns == 2) {
      fprintf(stdout,
              "THREAD# %d ->Highest vigilance for 2 = %0.4f for object at X = %d, Y = %d\n",
              i,
              highest_confidence[i][1],
              highx[i][1],
              highy[i][1]);
    }
  }
#endif
#if DISPSTARTUP
  printf("\nnumTraining: %d\n", numTraining);
#endif

  sort_list_items();

  ptr = head;
  while(ptr != NULL) {
    printf("\nResult: Matched image %d with confidence of %f at coordinates x = %d, y = %d",
           ptr->F2Neuron,
           ptr->Vigilance,
           ptr->x,
           ptr->y);
    ptr = ptr->next;
  }

#if TIMESTATS
  if(parallelStartTime->millitm >= startTime->millitm) {
    SerialSec = parallelStartTime->time - startTime->time;
    SerialMilli = parallelStartTime->millitm - startTime->millitm;
  }else{
    SerialSec = parallelStartTime->time - startTime->time - 1;
    SerialMilli = (1000 + parallelStartTime->millitm) - startTime->millitm;
  }
  if(endTime->millitm >= parallelStartTime->millitm) {
    ParallelSec = endTime->time - parallelStartTime->time;
    ParallelMilli = endTime->millitm - parallelStartTime->millitm;
  }else{
    ParallelSec = endTime->time - parallelStartTime->time - 1;
    ParallelMilli = (1000 + endTime->millitm) - parallelStartTime->millitm;
  }
  if(endTime->millitm >= startTime->millitm) {
    TotalSec = endTime->time - startTime->time;
    TotalMilli = endTime->millitm - startTime->millitm;
  }else{
    TotalSec = endTime->time - startTime->time - 1;
    TotalMilli = (1000 + endTime->millitm) - startTime->millitm;
  }

  printf("\n\nART-II Benchmark Runtime Performance Statistics:\n");
  printf("\nSerial Runtime:\t\t %ld.%ld sec",SerialSec,SerialMilli);
  printf("\nParallel Runtime:\t %ld.%ld sec",ParallelSec,ParallelMilli);
  printf("\nTotal Runtime:\t\t %ld.%ld sec",TotalSec,TotalMilli);

  percentSerial =
    ( (double)SerialSec + 0.001 * (double)SerialMilli ) / ( (double)TotalSec + 0.001 * (double)TotalMilli ) * 100;
  percentParallel =
    ( (double)ParallelSec + 0.001 * (double)ParallelMilli ) / ( (double)TotalSec + 0.001 * (double)TotalMilli ) * 100;
  printf("\n");
  printf("\nPercent Serial:\t\t %f %%",percentSerial);
  printf("\nPercent Parallel:\t %f %%",percentParallel);
  printf("\n");

#endif


  printf("\nDone\n");


  return 0;
Usage:
  fprintf(
    stderr,
    "Usage: scanner [-startx <num>] [-starty <num>] [-endx <num>] [-endy <num>] [-stride <num>] -scanfile <filename> -trainfile1 <filename> [-trainfile2 <filename>]\n");
  exit(1);
Version:
  fprintf(stderr,"Version 1.00 \n");
  exit(1);
}
