/* PageRank 
   The genetic.dat dataset comes from:
   http://www.cs.toronto.edu/~tsap/experiments/datasets/
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* allocate one object of given type */
#define NEW(type) ((type*)calloc((size_t)1,(size_t)sizeof(type)))

/* allocate num objects of given type */
#define NEW_A(num,type) ((type*)calloc((size_t)(num),(size_t)sizeof(type)))

typedef unsigned int u_int;

/* vector */
typedef struct
{
  u_int dim;
  double *e;
} VEC;

/* row of a sparse matrix */
typedef struct
{
  u_int  nnz;  /* # of non-zero (nz) value on this row */
  u_int  *col; /* column identifier for each nz value */
  double *val; /* value for each nz value */
} SROW;

/* sparse matrix */
typedef struct
{
  u_int m, n;
  SROW *row;
} SMAT;

/* v_new -- gets a VEC of dimension 'dim'
   Precondition: size >= 0
   Postcondition: initialized to zero */
VEC *v_new( u_int size )
{
  VEC *v;
  
  if( (v = NEW(VEC)) == (VEC *)NULL )
  {
    fprintf( stderr, "v_new memory error" );
    exit( -1 );
  }
  
  v->dim = size;
  if( (v->e = NEW_A(size,double)) == (double *)NULL )
  {
    free( v );
    fprintf( stderr, "v_new memory error" );
    exit( -1 );
  }
  
  return v;
}

/* v_free -- returns VEC & associated memory back to memory heap */
int v_free( VEC *v )
{
  if( v == (VEC *)NULL )
    return -1;
  
  if( v->e == (double *)NULL ) 
  {
    free( v );
  }
  else
  {
    free( v->e );
    free( v );
  }
  
  return 0;
}

/* sm_new -- gets an mxn sparse matrix by dynamic memory allocation 
   Precondition: m>=0 && n>=0
   Postcondition: each row is empty*/
SMAT *sm_new( u_int m, u_int n )
{
  SMAT *M;
  u_int i;
  
  if( (M = NEW(SMAT)) == (SMAT *)NULL )
  {
    fprintf( stderr, "sm_new memory error" );
    exit( -1 );
  }
  
  M->m = m ; M->n = n;

  if( (M->row = NEW_A(m,SROW)) == (SROW *)NULL )
  {
    free( M );
    fprintf( stderr, "sm_new memory error" );
    exit( -1 );
  }
  
  for( i = 0 ; i < m ; i++ )
  {
    (M->row[i]).nnz = 0;
    (M->row[i]).col = (u_int *) NULL;
    (M->row[i]).val = (double *) NULL;
  }

  return M;
}

/* sm_free -- returns SMAT & associated memory back to memory heap */
int sm_free( SMAT *M )
{
  u_int i;
  SROW *ri;

  if( M == (SMAT *)NULL ) return -1;
  
  if( M->row == (SROW *)NULL ) 
  {
    free( M );
  }
  else
  {
    for( i = 0 ; i < M->m ; i++ )
    {
      ri = &(M->row[i]);
      if( ri->nnz > 0 )
      {
        free( ri->col );
        free( ri->val );
      }
    }
    free( M->row );
    free( M );
  }
  
  return 0;
}

/* sm_input -- file input of sparse matrix 
   Precondition: will only work with a binary matrix. */
SMAT *sm_input( FILE *fp )
{
  SMAT *M;
  u_int *col; /* temp array to store the nz col of the current row */
  u_int r;    /* index of current row */
  int c;      /* index of current column */
  SROW *ri;   /* pointer to the current row in M */
  u_int m,n,i,j,k;
  
  /* get dimension */
  if( fscanf( fp, " SparseMatrix: %u by %u", &m, &n ) < 2 )
  {
    fprintf( stderr, "sm_input error reading dimensions" );
    exit( -1 );
  }
  
  if( (col = NEW_A(n,u_int)) == (u_int *)NULL )
  {
    fprintf( stderr, "sm_input memory error" );
    exit( -1 );
  }

  M = sm_new( m, n );
  
  /* get entries */
  for( i=0 ; i<m ; i++ )
  {
    if( fscanf( fp, " row %u:", &r ) < 1 )
    {
      fprintf( stderr, "sm_input error reading line %u", i );
      exit( -1 );
    }
    ri = &(M->row[i]);
    j = 0;
    for( ; ; )
    {
      if( fscanf( fp, "%d", &c ) < 1 )
      {
        fprintf( stderr, "sm_input error reading line %u col x", i );
        exit( -1 );
      }
      if( c < 0 ) break;
      col[j] = c;
      j++;
    } /* j is the number of nz value in row i */

    if( j!= 0 )
    {
      if( (ri->col = NEW_A(j,u_int)) == (u_int *)NULL )
      {
        fprintf( stderr, "sm_input memory error" );
        exit( -1 );
      }
      if( (ri->val = NEW_A(j,double)) == (double *)NULL )
      {
        fprintf( stderr, "sm_input memory error" );
        exit( -1 );
      }
    }

    ri->nnz = j;

    for( k = 0 ; k < j ; k++ )
    {
      ri->col[k] = col[k]; 
      ri->val[k] = 1.0;
    }
  }
  
  free( col );

  return M;
}

/* sm_output -- file output of sparse matrix 
   Postcondition: the result is not a valid entry for sm_input,
     since it also works for a non-binary matrix. */
void sm_output( FILE *fp, SMAT *M )
{
  u_int i,j;
  SROW *ri;

  fprintf( fp, "SparseMatrix: %d by %d\n", M->m, M->n );

  for( i = 0 ; i < M->m ; i++ )
  {
    fprintf( fp, "row %u: ", i ); 
    ri = &( M->row[i] );
    for( j = 0 ; j < ri->nnz ; j++ )
    {
      fprintf( fp, "%u:%1.5g ", ri->col[j], ri->val[j] );
    }
    fprintf( fp, "-1\n" );
  }
}

/* v_output -- file output of vector */
void v_output( FILE *fp, VEC *v )
{
  u_int i;

  fprintf( fp, "Vector: %d\n", v->dim );
  for( i = 0 ; i < v->dim ; i++ ) fprintf( fp, "%1.5g ", v->e[i] );
  putc( '\n', fp );
}

/* Calcul le produit d'un vecteur et d'une matrice creuse */
VEC* vec_prod_smat(VEC* base, SMAT* H)
{
  VEC* next = v_new(H->n); /* Vecteur résultat */
  SROW* crow = NULL; /* Pointeur de ligne */

  /* Initialise le vecteur résultat */
  for (int i=0; i < next->dim; i++) {next->e[i] = 0;}

  /* Parcours chaque ligne de la matrice */
  for (int i=0; i < H->m; i++)
  {
    /* Récupère la ligne courante */
    crow = &(H->row[i]);

    /* Ajoute chaque valeur non-nulles */
    for (int j=0; j < crow->nnz; j++)
    {
      next->e[crow->col[j]] += base->e[i] * crow->val[j];
    }
  }

  v_free(base);
  return next;
}

/* Affiche le contenu d'un vecteur */
void vec_print(VEC* v)
{
  if (v->dim < 1) 
  {
    printf("[]\n");
    return;
  }

  printf("[%f", v->e[0]);
  for (int i=1; i < v->dim; i++)
  {
    printf(", %f", v->e[i]);
  }
  printf("]\n");
}

int main()
{
  FILE *fp;
  SMAT *SM;

  SROW* crow = NULL; /* Pointeur de ligne */

  fp = fopen( "exemple.dat", "r" );
  SM = sm_input( fp );
  fclose( fp );

  /* Corrige les valeurs de SM */
  for (int i=0; i < SM->m; i++)
  {
    /* Récupère la ligne courante */
    crow = &(SM->row[i]);

    /* Combinaison linéaire */
    for (int j=0; j < crow->nnz; j++)
    {
      crow->val[j] /= crow->nnz;
    }
  } 

  /* Calcul de r0 */
  VEC* r0 = v_new(SM->n);
  for (int i=0; i < r0->dim; i++)
  {r0->e[i] = (double) 1/r0->dim;}

  VEC* r1 = r0;

  /* Itére sur la matrice */
  for (int i=0; i < 5; i++)
  {
    /* Calcul de Rk+1 */
    r1 = vec_prod_smat(r1, SM);
  }

  /* Affiche les résultats */
  vec_print(r1);

  // sm_output( stdout, SM );

  sm_free( SM );
  v_free(r1);

  return 0;
}
