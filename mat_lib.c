#include "mkl_cblas.h"
#include "mkl_lapack.h"
#include "mkl_lapacke.h"
#include "mkl_service.h"
#include "mkl_vml_functions.h"
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <mkl.h>
#include <omp.h>
#include <search.h>



#define limit 5200 //Limit of len to run serial or parallel




//-----------------------------Dynamic_Arrays(double)----------------------------------------------------

typedef struct
{
    size_t len; 
    double *data; 
}array;

void array_print(const array *input)
{
    printf("\n");
    for(int i=0; i< input->len; i++)
    {
        printf("[%d]:%lf ",i,input->data[i]);
    }
}

void array_init(size_t len,array *main_array)
{
    main_array->len = len;
    main_array->data=(double *)calloc(len, sizeof(*main_array->data));
}


//-----------------------------Dynamic_Arrays(int)----------------------------------------------------
typedef struct
{
    size_t len; 
    int *data; 
}iarray;

void iarray_print(const iarray *input)
{
    printf("\n");
    for(int i=0; i< input->len; i++)
    {
        printf("[%d]:%d ",i,input->data[i]);
    }
}

void iarray_init(size_t len,iarray *main_array)
{
    main_array->len = len;
    main_array->data=(int *)calloc(len, sizeof(*main_array->data));
}


//-----------------------------HASH-TABLE(MAP)----------------------------------------------------
//NOT RESPONSIBLE FOR INITIALIZATION OF INPUTS
//-----------------------------MATRIX-INDEX-DATA-STRUCTURE----------------------------------------------------

typedef struct 
{
    int x,y;
}mat_index;

static const mat_index origin={0,0};

//-----------------------------COORD_PAIR-DATA-STRUCTURE----------------------------------------------------

typedef struct 
{
    int key,x,y;
}coord_pair;

int compar_coord(const void *l, const void *r) //Function to compare keys
{
    const coord_pair *lm = l;
    const coord_pair *lr = r;
    return lm->key - lr->key;
}

void coord_init_list(coord_pair *input,const array *list_x,const array *list_y, void **root) //Just for syntax template (Not usefull)
{
    for(int i=0; i<list_x->len; i++) //Create a tree
    {
        input->key=i;
        input->x=list_x->data[i];
        input->y=list_y->data[i];
        if(tsearch((void *)input, root, compar_coord) == NULL){printf("Error at insertetion of nodes");} /* insert */
        input++;
    }
}

void coord_internal_print(const void *ptr, VISIT order, int level) //Used only by twalk
{
    const coord_pair *p = *(const coord_pair **) ptr;

    if (order == postorder || order == leaf)  {
        printf("key = %d,  value (x,y)= %d,%d\n", p->key, p->x,p->y);
    }
}


void coord_print_tree(void **root)
{
    twalk(root, coord_internal_print);
}

int coord_update(int key, double x,double y, void **root) {
    coord_pair *dummy =malloc(sizeof(coord_pair)); //Unfreed memory
    dummy->key=key;

    void *found = tfind(dummy, root, compar_coord);
    if (!found)
    {
        printf("Error at updating values");
        return -1;
    } 
    
    coord_pair *entry = *(coord_pair **)found;
    entry->x = x;
    entry->y=y;
    return 0;
}

coord_pair * coord_pair_insert(int key, int x,int y, void **root) { 
    coord_pair *p = malloc(sizeof(coord_pair));
    
    p->key = key;
    p->x= x;
    p->y= y;

    void *ret = tsearch(p, root, compar_coord);
    if (!ret) {
        free(p);
        return NULL;
    } 

   // If a duplicate was found, free your new node and return the existing one
    if (*(coord_pair **)ret != p) {
        free(p);
    } 
    return *(coord_pair **)ret;
}



int coord_pair_search(int key, double *out_x,double *out_y, void **root) {
    if (!root || !out_x || !out_y) 
    {
        printf("Error at input(coord_pair)");
        return NAN;
    }

    coord_pair dummy;
    dummy.key=key;
    void *found = tfind(&dummy, root, compar_coord);
    if (!found)
    {
        printf("Not found in this tree(coord_pair)");
        *out_x = NAN;
        *out_y = NAN;
        return NAN;
    } 
    

    const coord_pair *entry = *(coord_pair **)found;
    *out_x = entry->x;
    *out_y = entry->y;
    return 1;
}


static int counter = 0; //Works on any kind of tree(map)

void count_action(const void *node, VISIT which, int depth) {
    if (which == postorder || which == leaf) {
        counter++;
    }
}

int count(void **root) {
    if (!root || !*root) {
        printf("Empty tree (count)\n");
        return 0;
    }
    counter = 0; 
    twalk(*root, count_action); 
    return counter;
}


//-----------------------------PAIR-DATA-STRUCTURE----------------------------------------------------

typedef struct
{
   int key;
   double val; 
}pair;


int compar(const void *l, const void *r) //Function to compare keys
{
    const pair *lm = l;
    const pair *lr = r;
    return lm->key - lr->key;
}

void pair_init_list(pair *input,const array *list, void **root) //Just for syntax template (Not usefull)
{
    for(int i=0; i<list->len; i++) //Create a tree
    {
        input->key=i;
        input->val=list->data[i];
        if(tsearch((void *)input, root, compar) == NULL){printf("Error at insertetion of nodes");} /* insert */
        input++;
    }
}

void pair_internal_print(const void *ptr, VISIT order, int level) //Used only by twalk
{
    const pair *p = *(const pair **) ptr;

    if (order == postorder || order == leaf)  {
        printf("key = %d,  value = %lf\n", p->key, p->val);
    }
}


void pair_print_tree(void **root)
{
    twalk(root, pair_internal_print);
}


int pair_update( int key, double new_value, void **root) {
    pair dummy;
    dummy.key=key;

    void *found = tfind(&dummy, root, compar);
    if (!found)
    {
        printf("Error at updating values");
        return -1;
    }
    
    pair *entry = *(pair **)found;
    entry->val = new_value;
    return 0;
}

pair * pair_insert(int key, double val, void **root) { 
    pair *p = malloc(sizeof(pair));
    
    p->key = key;
    p->val = val;

    void *ret = tsearch(p, root, compar);
    if (!ret) {
        free(p);
        return NULL;
    }

    // If a duplicate was found, free your new node and return the existing one
    if (*(pair **)ret != p) {
        free(p);
    }
    return *(pair **)ret;
}

double pair_search(int key, void **root) 
{
    if (!root ) 
     {
        printf("Error at input");
        return NAN;
    }

    pair dummy;
    dummy.key=key;
    void *found = tfind(&dummy, root, compar);
    if (!found) 
    {
        printf("Not found in this tree(pair)");
        return NAN;
    }

    pair *entry = *(pair **)found;
    return entry->val;
}

//-----------------------------SORTING-FUCNTION(QUICK-SORT)----------------------------------------------------


void *global_root; //Global root for tree because compare func does not support extra arguments

int sorting_func(const void *key1, const void *key2)
{
    if (pair_search(*(int *)key1, &global_root) < pair_search(*(int *)key2, &global_root))
        return -1;  
    else if (pair_search(*(int *)key1, &global_root) > pair_search(*(int *)key2, &global_root))
        return 1; 
    return 0;
}



//-----------------------------Polynomials----------------------------------------------------
//!RULE :ALWAYS THE FUNCTION IS RESPONSIBLE FOR POLYNOMIAL INITIALIZATION


typedef struct
{
    size_t degree; 
    double *data; //[0]=data[0]*x^0+data[1]*x^1+...+data[d]*x^{d} x does not matter since it will be a Matrix
}poly;

void init_poly(size_t degree,poly *main_poly)
{
    main_poly->degree = degree;
    main_poly->data=(double *)calloc(degree+1, sizeof(*main_poly->data));
    
}
#define term(m,i) ((m).data[i]) //omit .data every time
#define termr(m,i) ((m)->data[i]) //omit -> every time (For refrences only)

void print_poly(const poly *main_poly)
{
    for (u_int i=0; i<main_poly->degree+1; i++)
        printf("%lf x^%d+",termr(main_poly, i),i);      
    printf("\n"); 
}

void read_file_poly(const char* filename,poly *main_poly)
{
    FILE* file = fopen(filename, "r");
    if (!file) 
    {
        fprintf(stderr, "Failed to open file");
        exit(EXIT_FAILURE);
    }
    u_int degree;
    if(fscanf(file, "%u ",&degree) != 1 )
    {
        fprintf(stderr, "Wrong input dimensions");
        exit(EXIT_FAILURE);
    } 
    init_poly(degree,main_poly);
    for (u_int i=0; i<=degree; i++) {
            if(fscanf(file, "%lf",&termr(main_poly,i)) != 1 )
            {
                fprintf(stderr, "Error at reading input polynomial");
                exit(EXIT_FAILURE);
            }
    }
    fclose(file);
}

void free_poly(poly *main_poly)
{
    free(main_poly->data);
}

//-------------------------MATRIX--------------------------------------------------------

//CONSIDER ONLY SQUARE MATRICES!!!!
//!RULE :ALWAYS THE FUNCTION IS RESPONSIBLE FOR MATRIX INITIALIZATION

typedef struct 
{
    int rows;
    int cols;
    size_t len;
    double *data;
    int leading_dim;
}mat;

void init_mat(u_int rows,u_int cols,mat *main_mat)
{
    main_mat->rows=rows;
    main_mat->cols=cols;
    main_mat->len=(size_t)cols * (size_t)rows;
    main_mat->data=(double *)calloc(main_mat->len, sizeof(*main_mat->data)); //All elements are initialized to 0
}

#define cell(m,i,j) ((m).data[(j)*(m).rows + (i)]) //Access 1d array like 2d column_major
#define cellr(m,i,j) ((m)->data[(j)*(m)->rows + (i)]) //Access 1d array like 2d (For refrences only)

void freemat(mat *main_mat)
{
    if(main_mat)
    {
        free(main_mat->data);
        main_mat->data=NULL;
    }
}


void mat_mult_gemm(const mat *A,const mat *B,mat *C) //Using mkl dgemm is already optimized see page 24 Developer Refrence MLK
{
    if(A->cols != B->rows)
    {
        fprintf(stderr, "Wrong dimensions");
        exit(EXIT_FAILURE);
    }

    init_mat(A->rows, B->cols, C);
    const u_int m = A->rows;
    const u_int k = A->cols;
    const u_int n = B->cols;
    const double alpha = 1;
    const double beta = 0;
    // It is not so fast
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha,
                A->data, k, B->data, n, beta, C->data, n);
}

void mat_mult_options(const mat *A,const mat *B,mat *C,const int m,const int n,const int k,mat_index fromA,mat_index fromB,mat_index fromC,double alpha_in,double beta_in) //More advanced matrix multiplication that support submatrices and -= or += to final result 
{
 
    //DOES INTITIALIE THE RESULT
    const double alpha = alpha_in; //For addition 1 , for subtraction -1 , scallar to the product Α * Β
    const double beta = beta_in; //Scallar for the C matrix 

    if (fromA.x + m > A->rows || fromA.y + k > A->cols ) 
    {
        printf("A submatrix out of bounds");
        return;
    }
    if (fromB.x + k > B->rows || fromB.y + n > B->cols ) 
    {
        printf("B submatrix out of bounds");
        return;
    }
     if (fromC.x + m > C->rows || fromC.y + n > C->cols ) 
    {
        printf("C submatrix out of bounds");
        return;
    }  

/*     const double *Ap_start = A->data + fromA.y + A->cols * fromA.x;
    const double *Bp_start = B->data + fromB.x * B->cols + fromB.y;
    double *Cp_start = C->data + fromC.x * C->cols + fromC.y;  */
    


    const double *Ap_start = A->data + fromA.y * A->rows + fromA.x;
    const double *Bp_start = B->data + fromB.y * B->rows + fromB.x;
    double *Cp_start       = C->data + fromC.y * C->rows + fromC.x; 

    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha,
                Ap_start, k, Bp_start, n, beta, Cp_start, n);
}







void mat_mult_naive(const mat *A,const mat *B,mat *C) 
{
    if(A->cols != B->rows)
    {
        fprintf(stderr, "Wrong dimensions");
        exit(EXIT_FAILURE);
    }
    #pragma omp parallel for collapse(2)
    for (int i=0; i<A->rows; i++) {
        for(int j=0; j<A->cols; j++ )
        {
            cellr(C, i, j) = 0;
            for (int k = 0; k < A->rows; k++)
            {
                cellr(C, i, j) += cellr(A, i, k) * cellr(B, k, j);
            }
        }
        }
}





void read_file_mat(const char* filename,mat *main_mat)
{
    FILE* file = fopen(filename, "r");
    if (!file) 
    {
        fprintf(stderr, "Failed to open file");
        exit(EXIT_FAILURE);
    }
    u_int rows,cols;
    if(fscanf(file, "%u %u",&rows,&cols) != 2 )
    {
        fprintf(stderr, "Wrong input dimensions");
        exit(EXIT_FAILURE);
    } 
    init_mat(rows, cols,main_mat);
    for (u_int i=0; i<rows; i++) {
        for (u_int j=0; j<cols; j++) {
            if(fscanf(file, "%lf",&cellr(main_mat,i,j)) != 1 )
            {
                fprintf(stderr, "Error at reading input mstrix");
                exit(EXIT_FAILURE);
            }
        }
    }
    fclose(file);
}

void print_mat(const mat *main_mat)
{
    for (u_int i=0; i<main_mat->rows; i++) {
        for (u_int j=0; j<main_mat->cols; j++) {
            
                printf("%8.6lf ",cellr(main_mat, i, j));
        }
        printf("\n");
    }
    printf("\n");
}

void mat_copy(const mat *A,mat *B) //B=A
{
    init_mat(A->rows, A->cols, B);
    cblas_dcopy(A->len, A->data, 1, B->data, 1); //Might not be optimal for large matrices
}

void mat_copy_options(const mat *A,mat *B,const mat_index from) 
{
    //Copy the whole matrix A into B starting from from(x,y)(On the B matrix)

    if (from.x + A->rows > B->rows || from.y + A->cols > B->cols )
    {
         printf("Cannot assign matrix, out of bounds lenA:%zu -> lenB:%zu + offset:%d",A->len,B->len,(from.x + from.y * B->rows));
        return;
    }
    //Calculate offset +(from.x + from.y * A->rows)
    int i=LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'A', A->rows, A->cols, A->data, A->rows, B->data+(from.x + from.y * B->rows), B->cols);
    if(i!=0)
    {
        printf("Error at copying:%d",i);
    }
}

/* void mat_copy_options_smaller(const mat *A,mat *B,const mat_index from) 
{
    //Copy the matrix A starting from(x,y) into B 
    //Calculate offset +(from.x + from.y * A->rows)

    if (A->len-(from.x + from.y * A->rows) > B->len )
    {
        printf("Cannot assign matrix, out of bounds (smaller) lenA:%zu -> lenB:%zu",A->len,B->len);
        return;
    }

    int i=LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'A', A->rows, A->cols, A->data+(from.x + from.y * A->rows), A->rows, B->data, B->cols);
    if(i!=0)
    {
        printf("Error at copying:%d(smaller)",i);
    }
} */


void mat_copy_options_smaller(const mat *A, mat *B, const mat_index from)
{
    // Copy submatrix from A starting at (from.x, from.y) into B

    /* if (from.x < 0 || from.y < 0 || from.x + B->rows > A->rows || from.y + B->cols > A->cols)
    {
        fprintf(stderr, "Submatrix copy out of bounds: from=(%d,%d), B=(%d,%d), A=(%d,%d)\n",
                from.x, from.y, B->rows, B->cols, A->rows, A->cols);
        exit(EXIT_FAILURE);
    } */

    const double *A_start = A->data + from.y * A->rows + from.x;

    int info = LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'A',
                              B->rows, B->cols,
                              A_start, A->rows,
                              B->data, B->rows);

    if (info != 0)
    {
        printf("Error during LAPACKE_dlacpy: %d (smaller)\n", info);
    }
}
 

void mat_scalar_mult(const mat *main_mat,const double scalar,mat *result_mat)
{
    init_mat(main_mat->rows, main_mat->cols, result_mat);
    mat_copy(main_mat, result_mat);
    cblas_dscal(result_mat->len, scalar, result_mat->data, 1);
}

void mat_add(const mat *A,const mat *B,mat *result_mat)
{
    init_mat(A->rows, A->cols, result_mat);
    vdAdd(A->len, A->data, B->data, result_mat->data);
}

void mat_add_openmp(const mat *A,const mat *B,mat *result_mat)
{
    init_mat(A->rows, A->cols, result_mat);
    #pragma omp parallel for if(A->rows>limit)
    for (int i=0; i<A->rows; i++) {
        for(int j=0; j<A->cols; j++ )
        {
            cellr(result_mat, i, j) = cellr(A, i, j) + cellr(B, i, j);
        }
    }
}

void mat_identity(u_int rows,u_int cols, mat *A)
{
    init_mat(rows,cols,A);
    int i=LAPACKE_dlaset(LAPACK_COL_MAJOR, 'A', A->rows, A->cols, 0, 1, A->data, A->cols);
    if(i!=0)
    {
        printf("Error at identity matrix initialization");
    }
}

//-----------------------------------------Polynomial evaluation------------------------------------------------


void eval_poly_mat_naive(const mat *A,poly *pol,mat *result_mat)
{
    init_mat(A->rows, A->cols, result_mat);
    mat I; mat_identity(A->rows, A->cols, &I); //Intentity matrix;
    mat temp,temp2,temp3; mat_copy(A, &temp); mat_copy(A, &temp2);init_mat(A->rows, A->cols, &temp3);
    cblas_dscal(I.len, termr(pol, 0), I.data, 1); //c0*I=I
    cblas_dscal(temp.len, termr(pol, 1), temp.data, 1); //c1*A=A
    vdAdd(A->len, I.data, temp.data, result_mat->data); //c0*I+c1*A
    //mat_copy(A, &temp); //Clear temp
    for(int i=2; i<=(pol->degree); i++)
    {
        mat_mult_gemm(A, &temp, &temp2); //temp=P  
        
        cblas_dcopy(temp2.len, temp2.data, 1, temp.data, 1); //P=P * A
        cblas_dscal(temp2.len, termr(pol, i), temp2.data, 1); //P=c_i * P
        vdAdd(result_mat->len, temp2.data, result_mat->data, temp3.data); // result_mat=Q
        mat_copy(&temp3, result_mat); //Q=Q+P
        //print_mat(&temp3);
        memset(temp3.data, 0, temp3.len * sizeof(double));
        free(temp2.data);
    }
}




void eval_poly_mat_paterson_stockmeyer(const mat *A,poly *pol,mat *qA) //We dont care about dimension everything is square and share the same dimensions
{
    if (pol->degree < 3)
    {
        eval_poly_mat_naive(A, pol, qA);
        return;
    }

    init_mat(A->rows, A->cols, qA); //Every function is responsilbe for initializing the results
    mat I; mat_identity(A->rows, A->cols, &I); //Identity matrix;
    mat Q; init_mat(A->rows, A->cols, &Q); //Intermediate results
    

    size_t d = pol->degree ;
    size_t p = (size_t)sqrt(d);
    size_t s = (d + 1)/p;

    mat* powersA=(mat *)malloc((p+1) * sizeof(*A)); //might need p-1 allocations
    init_mat(I.rows, I.cols, &powersA[0]);
    cblas_dcopy(I.len, I.data, 1, powersA[0].data, 1);
    //mat_copy_options(&I, &powersA[0], origin);
    mat_copy(A, &powersA[1]); //A^1
    mat temp;init_mat(A->rows, A->cols, &temp);
    
    //powersA[1].len=I.len; powersA[1].cols=I.cols; powersA[1].rows=I.rows; 
    //Evaluate and cache I,A,A^2,...,A^{p-1}
    for(int i=2; i<=p; i++) //Cannot be parallized since we have dependency between A_i and A_{i-1}
    {
        //Initialization phase
        init_mat(I.rows, I.cols, &powersA[i]);
        cblas_dcopy(powersA[i-1].len, powersA[i-1].data, 1, powersA[i].data, 1); 

        cblas_dgemm_64(CblasColMajor, CblasNoTrans, CblasNoTrans, (const int)A->cols, (const int)A->cols, (const int)A->cols, 1,
                A->data, (const int)A->cols, powersA[i].data, (const int)A->cols, 0, temp.data, (const int)A->cols);

        cblas_dcopy(powersA[i].len, temp.data, 1, powersA[i].data, 1); //Copy temp back 
        memset(temp.data, 0, temp.len * sizeof(double)); //temp reset          
    }
   
    const int incx=1,incy=1,len=A->len; 
    const double a=1;
    double *tempfor = malloc(Q.len * sizeof(double));
    
    for(int si=s-1; si>=0; si--)
    {
        #pragma omp parallel for //thread safety from mkl
        for(size_t i=0; i<p; i++)
        {
            int pi= si * p + i; 
            cblas_daxpy(len, termr(pol, pi), powersA[i].data, incx, tempfor, incy); //tempfor=c_{q*q} * I , c_{q*(p+1)} * A , ... ,
            cblas_daxpy(len, a, tempfor, incx, Q.data, incy); //Q_q+=tempfor=(c_{q*q} * I + c_{q*(p+1)} * A + ... +)                   
        }
        
        cblas_dgemm_64(CblasColMajor, CblasNoTrans, CblasNoTrans, (const int)A->cols, (const int)A->cols, (const int)A->cols, 1,
                qA->data, (const int)A->cols, powersA[p].data, (const int)A->cols, 0, temp.data, (const int)A->cols); //temp=qA * A^p
        
        vdAdd(temp.len, temp.data, Q.data, qA->data);//qA=temp+Q
        memset(temp.data, 0, temp.len * sizeof(double)); //temp reset 
        memset(Q.data, 0, Q.len * sizeof(double)); //Q reset
    } 

    //Free memory section
    free(tempfor);
    for (int i=0; i<=p; i++) {
        //print_mat(&powersA[i]); 
        free(powersA[i].data);  
    }  
    free(powersA);
    free(temp.data);
    free(I.data);
    free(Q.data);
    
}


void eval_poly_mat_paterson_stockmeyer_options(const mat *A,const poly *pol,mat *qA,const mat_index from,const int n)
{

    init_mat(n, n, qA); //Every function is responsilbe for initializing the results
    mat I; mat_identity(n, n, &I); //Identity matrix;
    mat Q; init_mat(n, n, &Q); //Intermediate results
    

    size_t d = pol->degree ;
    size_t p = (size_t)sqrt(d);
    size_t s = (d + 1)/p;

    mat* powersA=(mat *)malloc((p+1) * sizeof(*A)); //might need p-1 allocations
    init_mat(I.rows, I.cols, &powersA[0]);
    cblas_dcopy(I.len, I.data, 1, powersA[0].data, 1);

    init_mat(I.rows, I.cols, &powersA[1]);
    //mat_copy_options(A, &powersA[1], from);
    mat_copy_options_smaller(A, &powersA[1], from);
    //print_mat(&powersA[1]);
    mat temp;init_mat(n, n, &temp);
    
    //powersA[1].len=I.len; powersA[1].cols=I.cols; powersA[1].rows=I.rows; 
    //Evaluate and cache I,A,A^2,...,A^{p-1}
    for(int i=2; i<=p; i++) //Cannot be esasily parallized since we have dependency between A_i and A_{i-1}
    {
        //Initialization phase
        init_mat(I.rows, I.cols, &powersA[i]);
        cblas_dcopy(powersA[i-1].len, powersA[i-1].data, 1, powersA[i].data, 1); 
        

        mat_mult_options(&powersA[1], &powersA[i], &temp,  n, n, n, from, origin, origin, 1, 0);
        cblas_dcopy(powersA[i].len, temp.data, 1, powersA[i].data, 1); //Copy temp back 
        memset(temp.data, 0, temp.len * sizeof(double)); //temp reset          
    }
   
    const int incx=1,incy=1,len=Q.len; 
    const double a=1;
    double *tempfor = calloc(Q.len,sizeof(double));
    for(int si=s-1; si>=0; si--)
    {
        #pragma omp parallel for //MKL gurantees thread safety
        for(size_t i=0; i<p; i++)
        {
            int pi= si * p + i;
            
            cblas_daxpy(len, termr(pol, pi), powersA[i].data, incx, tempfor, incy); //tempfor=c_{q*q} * I , c_{q*(p+1)} * A , ... ,
            cblas_daxpy(len, a, tempfor, incx, Q.data, incy); //Q_q+=tempfor=(c_{q*q} * I + c_{q*(p+1)} * A + ... +)                  
           
        }
        mat_mult_options(qA, &powersA[p], &temp, qA->rows, Q.cols, temp.rows, origin, origin, origin, 1, 0);
        
        vdAdd(temp.len, temp.data, Q.data, qA->data);//qA=temp+Q
        memset(temp.data, 0, temp.len * sizeof(double)); //temp reset 
        memset(Q.data, 0, Q.len * sizeof(double)); //Q reset
    } 
    //Free memory section
     free(tempfor); 
    for (int i=0; i<=p; i++) {
        //print_mat(&powersA[i]); 
        free(powersA[i].data);  
    }  
    free(powersA);
    free(temp.data);
    free(I.data);
    free(Q.data);
}


//----------------------------------------------DECOMPOSIOTIONS-TRANSFORMATIONS------------------------------------------

void schur_decomposition(const mat *A, mat *schur, mat *schur_co, array *eigenvalues_real, array *eigenvalues_complex) //No need for complex eigenvalues
{
    //must initialize the output 
    mat_copy(A, schur);
    init_mat(A->rows, A->cols, schur_co);
    array_init(A->rows, eigenvalues_complex);
    array_init(A->rows, eigenvalues_real);
    const int n = A->rows;
    const int leading_dim = n > 1 ? n : 1;
    int info=-2; //used see if execution is successful
    int lwork=-1;
    double work_size=0;
    int sdim = 0;
    //Calculate the len for workspace
    dgees("V", "N", NULL, &n, schur->data, &leading_dim, &sdim, eigenvalues_real->data, eigenvalues_complex->data, schur_co->data, &leading_dim, &work_size, &lwork, NULL, &info);
    lwork = (int)work_size;
    array work; array_init(lwork, &work);
    //Do the actual schur decomposition
    dgees("V", "N", NULL, &n, schur->data, &leading_dim, &sdim, eigenvalues_real->data, eigenvalues_complex->data, schur_co->data, &leading_dim, work.data, &lwork, NULL, &info);
    //printf("\nInfo:%lld\n",info);
    //array_print(&work);
    if(info != 0){printf("Problem at schur decomposiotion");} //For more see page 1062 MLK developer refernce
}


void mat_transpose(mat *A)
{
    #pragma omp parallel for collapse(2)
    for (int r=0; r< A->rows; r++)
    {
        for (int c=r+1; c<A->cols; c++) //Swap only the upper triangular with the lower triangular
        {
            double temp =cellr(A, r, c); //In place swap does not work 
            cellr(A, r, c)=cellr(A, c, r);
            cellr(A, c, r)=temp;
        }
    }
}

void mat_inverse(mat *A)
{
   //array ipivc;array_init(A->rows, &ipivc);
   int *ipiv = malloc(A->len * sizeof(long long));
   const int lwork= A->len;
   array work;array_init(A->rows, &work);
   int infoLU=-2,infoinv=-2;
   const int n= A->rows;
   const int m= A->cols;

   //First we need LU decomposition
   dgetrf(&n, &m, A->data, &n, ipiv, &infoLU);

   //Then calculate the inverse
   dgetri(&n, A->data, &n, ipiv, work.data, &lwork, &infoinv);
   if (infoinv != 0 || infoLU != 0) {printf("Problem at inverse function");}
}

//----------------------------------------------EIGENVALUES-CLUSTERS------------------------------------------
//ONLY REAL EIGENVALUES
//FUNCTIONS ARE NOT RESPOSIBLE FOR INITIALIZATIONS

void cluster_values(const double *values,const int n,const double delta,int *clusters)
{
    
    int p=1;
    
    for (int i = 0; i < n; i++)
	{
        double lamda_i= values[i];
        if (clusters[i] == 0 || clusters[i] >= p) // q == p 
		{
			clusters[i] = p;
			p++;
		}
        int qi=clusters[i];
        for (int j = i + 1; j < n; j++)
		{
            double lamda_j= values[j];
            if(clusters[j] != clusters[i])
            {
                double distance=0;
                double diff= lamda_j - lamda_i;
                distance=fabs(diff); // |lambda_j - lambda_i|
                if (distance <= delta) 
				{
                    qi = clusters[i];
                    if (clusters[j] == 0 || clusters[j] >= p) 
                    {
                        clusters[j] = qi;
                    }
                    else 
                    {
                        int qj = clusters[j];

						int min_qiqj = (qi < qj) ? qi : qj;
						int max_qiqj = (qi > qj) ? qi : qj;

                        #pragma omp parallel for 
                        for(int m=0; m<n; m++)
                        {
                            if (clusters[m] == max_qiqj)
								clusters[m] = min_qiqj; // Move the elements of Smax(qi,qj ) to Smin(qi,qj )
							if (clusters[m] > max_qiqj)
								clusters[m]--; // Reduce by 1 the indices of sets Sq for q > max(qi, qj ).
                        }
                        p--;
                    }
                }
            }
        }
    }     
}

iarray cluster_permutation(const int *clusters,const int n) //Need the size of th returned array
{
    void *root_index=0,*root_average=0;//Roots for tree construct that acts like map
     //coord_pair *cluster_index=malloc(sizeof(coord_pair) * n);
    //pair *cluster_average=malloc(sizeof(pair) * n);
    coord_pair *last_inserted_index = NULL; //Substitute for .end
    const pair *last_inserted_average = NULL; //Substitute for .end

    for (int i=0; i<n; i++) 
    {
      coord_pair *ptr_index = malloc(sizeof(coord_pair)); //Unfreed memory 
      ptr_index->key= clusters[i];
      if(tfind(ptr_index, &root_index, compar_coord) == NULL) //We need NULL return
      {
         last_inserted_index = coord_pair_insert(clusters[i], 0, 0, &root_index);
         //printf("key:%d val x:%d val y:%d \n",last_inserted_index->key,last_inserted_index->x,last_inserted_index->y);
      }
      //coord_print_tree(root_index);
      double x,y;
      coord_pair_search(clusters[i], &x, &y, &root_index);
      x++; //first
      y+= (i+1); //second
      coord_update(clusters[i],  x,  y, &root_index);

      pair_insert(clusters[i], y/x, &root_average); //Add a new pair average map
    }
    //coord_print_tree(root_index);
    //pair_print_tree(root_average);
    int nclusters=count(&root_index); //Enumerate how many elements you have
    //*nclusters_return = nclusters; //For the return type
    iarray ordered_clusters;
    iarray_init(nclusters, &ordered_clusters);
    //int *ordered_clusters= calloc(nclusters, sizeof(int));
    #pragma omp parallel for
    for (int i=0; i<nclusters; i++) 
    {
        ordered_clusters.data[i]=i+1;
    }
    global_root=root_average;
    qsort(ordered_clusters.data, nclusters, sizeof(int), sorting_func);
    return ordered_clusters;
}

void bubblesort_eigenvalues(const mat *schur,const mat *schur_co,int n,int *clusters,const iarray *permutation) //Modifies the clusters
{
    iarray permutation_indexes; iarray_init(permutation->len, &permutation_indexes);
    for(unsigned int i=0;i< permutation->len; i++)
    {
        permutation_indexes.data[permutation->data[i]-1]=i;
    }
    int leading_dimension = n > 1 ? n : 1;

    // Bubble sort
    int temp_n=n;
    int swaps_count = 0, total_swaps =0;
    while (temp_n > 0)
    {
        int new_n=0;
        for (int i=1; i < temp_n; i++) 
        {
            int cluster_i_1 = clusters[i-1];
            int cluster_i = clusters[i];
            if(permutation_indexes.data[cluster_i_1 - 1] > permutation_indexes.data[cluster_i-1])
            {
                // Replace in T

                int ifst = i - 1 + 1;
				int ilst = i + 1;
				int orig_ifst = ifst, orig_ilst = ilst;
				int ifst_cluster = clusters[ifst - 1];
				int ilst_cluster = clusters[ilst - 1];

                int result = 0;
                array work; array_init(n, &work);
                dtrexc("V", &n, schur->data, &leading_dimension, schur_co->data, &leading_dimension, &ifst, &ilst, work.data, &result);

                if(result != 0)
                {
                    fprintf(stderr,"Error at bubblesort with dtrexc");
                    return ;
                }
                swaps_count++;
				total_swaps += ilst - ifst;

                clusters[ifst - 1] = ilst_cluster;
				clusters[ilst - 1] = ifst_cluster;
				if (ifst == orig_ifst - 1) // first block (now in ilst) was 2*2, need to set the cluster for bottom value
					clusters[ilst + 1 - 1] = ifst_cluster;
				if (ilst != orig_ilst - 1) // second block (now in ifst) was 2*2, need to set the cluster for bottom value
					clusters[ifst + 1 - 1] = ilst_cluster;
				new_n = i;
            }
        }
        temp_n = new_n;
    }
    //return clusters;
}



iarray cluster_schur_eigenvalues(const mat *schur,const mat *schur_co,array *eigenvalues,const double tolerance) //Returns the clusters
{
    int n=schur->rows;
    if (eigenvalues->len != n) //Dont khow why, schur decomposition seems to always return the right amount eigenvalues
	{
		array_init(n, eigenvalues);
		for (int i = 0; i < n; i++)
        {
			eigenvalues->data[i]= cellr(schur, i, i);
        }
	}
    //int *clusters= malloc(sizeof(int) * n);
    iarray clusters; iarray_init(n, &clusters);
    cluster_values(eigenvalues->data, n, tolerance, clusters.data);
    //int nclusters; //Total number for permutation array
    //int *permutation= cluster_permutation(clusters, n, &nclusters); //test the function
    iarray permutation=cluster_permutation(clusters.data, n); 
    int leading_dimension = n > 1 ? n : 1;

    bubblesort_eigenvalues(schur, schur_co, n, clusters.data, &permutation); 
    return clusters;
}






void f(const mat *A,const poly *ax,const int n,const mat_index from,mat *output)
{
    mat qA; eval_poly_mat_paterson_stockmeyer_options(A, ax, &qA, from, n);
    mat_copy_options(&qA, output, from);
}

typedef struct 
{
    mat *matrix;
    mat_index from;
    int rows;
    int cols;
    int dim;
    double *data; //Must a special way to intilizae dat using the mat_index from
    int len_array;
}matblock;


void matblock_init(const mat *matrix_in,const mat_index from,const int rows,const int cols,matblock *out)
{
    out->matrix=malloc(sizeof(mat));
    out->matrix->rows=matrix_in->rows;
    out->matrix->cols=matrix_in->cols;
    out->matrix->len= (size_t)out->matrix->cols * (size_t)out->matrix->cols;
    out->matrix->data=(double *)calloc(out->matrix->len, sizeof(double));
    cblas_dcopy(matrix_in->len, matrix_in->data, 1, out->matrix->data, 1);
    
    out->from =from;
    out->rows = rows;
    out->cols=cols;
    out->dim=matrix_in->leading_dim;
    out->data=matrix_in->data+(from.x + from.y * matrix_in->rows);
    //out ->len_array =0;//Used in case of array of blocks;
    
}



void handle_block(const mat *T,const poly *pol, const int cluster_start, const int cluster_end,int n,void f(const mat*,const poly*,const int,mat_index,mat*),mat *F,matblock *blocks)
{
    
    int block_size= cluster_end - cluster_start;
    int leading_dim=n; //Leading dimension for the block data
    mat block_data; init_mat(block_size, block_size, &block_data);
    block_data.leading_dim=n;
    block_data.data=T->data+(cluster_start + cluster_start * leading_dim);

    mat_index from= {cluster_start,cluster_start};
    
    
    matblock *new=realloc(blocks, (blocks[0].len_array +1) * sizeof(matblock));
    blocks[0].len_array++;
    if(new == NULL)
    {
        fprintf(stderr, "Error at realloc");
    }
    blocks=new;

    matblock_init(&block_data, from, block_size, block_size, &blocks[blocks[0].len_array]); 
    f(T,pol,cluster_end-cluster_start,from,F);
    
    //freemat(&block_data);
}




//RETURNS AN ARRAY OF MATRICES
matblock* calculate_diagonal_blocks(const mat *T,const poly *pol,const iarray *clusters,void f(const mat*,const poly*,const int,mat_index,mat*),mat *F)
{
    int n= T->rows;
    int cluster_start=1;
    matblock *blocks=malloc(1*sizeof(matblock)); //Not intialized YET
    blocks[0].len_array=1; //Used in handle block
    for (unsigned int i=0; i < clusters->len; i++)
    {
        if( clusters->data[cluster_start] != clusters->data[i])
        {
            int cluster_end=i;
            handle_block(T, pol, cluster_start, cluster_end, n, f, F, blocks);
            cluster_start =cluster_end;
        }
    }
    handle_block(T, pol, cluster_start, clusters->len, n, f, F, blocks);
    return blocks;
}


void sylvester_solver(const matblock *blocks,mat *F,const mat *T)
{
    int nblocks= blocks->len_array;
    for (int j=1; j< nblocks; j++)
    {
        for(int i=j-1; i>=0; i--)
        {
            mat *Tii=blocks[i].matrix; //Does not copy right 
            mat *Tjj=blocks[j].matrix; //Does not copy right
            int i_row = blocks[i].from.x;
            int j_col = blocks[j].from.y;

            mat_index ii = blocks[i].from;
            mat_index ij = {i_row,j_col};
            mat_index jj = blocks[j].from;

            int block_m = Tii->rows;
            int block_n = Tjj->cols;

            //Multiplications
            //same
            mat_mult_options(F, T, F, block_m, block_n, block_m, ii, ij, ij, 1, 1);
            mat_mult_options(T, F, F, block_m, block_n, block_n, ij, jj, ij, 1, -1);

            for(int k=i+1; k<j; k++)
            {
                int k_row = blocks[k].from.x;
                int k_col = blocks[k].from.y;

                int ik_m = blocks[i].rows;
                int ik_n = blocks[k].rows;

                int kj_n = blocks[j].rows;

                mat_index ik = {i_row,k_col};
                mat_index kj = {k_row,j_col};

                mat_mult_options(F, T, F, ik_m, kj_n, ik_n, ik, kj, ij, 1, 1);
                mat_mult_options(T, F, F, ik_m, kj_n, ik_n, ik, kj, ij, 1, -1);
            }
            mat *Fij; Fij->rows=block_m; Fij->cols=block_n; Fij->data=F->data+(ij.x+ij.y*T->rows);
            dtrsyl("N", "N", (const int *)-1, &Tii->rows, &Tjj->cols, Tii->data,&Tii->rows , Tjj->data,&Tjj->rows , Fij->data, &T->rows , (double *)1, NULL);            
        }
    }

}





void  calculate_block_parlett_recurrence(const mat *T,const poly *pol,iarray *clusters, void f(const mat*,const poly*,const int,mat_index,mat*),mat *out)
{
    int n= T->rows;
    if(n != T->cols)
    {
        fprintf(stderr, "Error this algorithm work only on square matrices");
    }
    init_mat(n, n, out);
    const matblock *blocks=calculate_diagonal_blocks(T, pol, clusters, f, out);
    sylvester_solver(blocks, out, T);
}



void calculate_polynomial(const mat *raw_input,const poly *pol,const double eigenvalues_cluster_delta,mat *result)
{
    mat output,schur,schur_co; init_mat(raw_input->rows, raw_input->cols, &output); //Initialize output
    array eigenvalues,eigenvalues_complex;
    schur_decomposition(raw_input, &schur, &schur_co, &eigenvalues, &eigenvalues_complex);

    int n=schur.rows;    
    iarray clusters= cluster_schur_eigenvalues(&schur, &schur_co, &eigenvalues, eigenvalues_cluster_delta);
    
    const double *co=pol->data;
    mat Ft; calculate_block_parlett_recurrence(&schur, pol, &clusters, f, &Ft);
    mat temp; init_mat(n, n, &temp); init_mat(n, n, result);
    mat_mult_gemm(&schur_co, &Ft, &temp);
    mat_mult_gemm(&schur_co, &temp, result);
}




int main(int argc, char* argv[])
{
     if (argc < 2) {
       fprintf( stderr,"Usage: %s filename",argv[0]);
        return 1;
    }
    // Input section
    mat A,C,C2,C3,I,schur,schur_co;
    array eigenvalues_real;
    read_file_mat(argv[1],&A);
    //read_file_mat("Amat2048",&A);
    //read_file_mat("A",&A);  
    //ead_file_mat("Amat6",&A);  

    poly ax;
    read_file_poly("polynomial.in",&ax);
  
   /*  print_poly(&ax);
    printf("\n"); */
    mkl_set_dynamic;
    mkl_verbose(1); //Used to print the time mlk functions
    FILE* debug = fopen("Debug.txt", "w");fclose(debug); //Erase the contents of file because verbose just appends
    mkl_verbose_output_file("Debug.txt");
   
    
    unsigned long long clock_start,clock_end;
    double start=dsecnd(); //Start time
    mkl_get_cpu_clocks(&clock_start); //Start cpu clock time
    
    //eval_poly_mat_paterson_stockmeyer(&A, &ax, &C2);
    //eval_poly_mat_naive(&A, &ax, &C2);
    //print_mat(&C2); 
    //eval_poly_mat_naive(&A, &ax, &C3);
    //print_mat(&C3);
    mat result;
    calculate_polynomial(&A, &ax, 0.005,&result);


    double end=dsecnd(); //Stop time
    mkl_get_cpu_clocks(&clock_end); //Stop cpu clock time 
    //print_mat(&C);
    printf("\n");
    

    printf("\n Time:%lf Cpu:%lf Ghz Cpu clock:%llu degree:%zu matrix:%d \n",end-start,mkl_get_clocks_frequency(),clock_end-clock_start,ax.degree,A.cols);

    // Free memory
    freemat(&A);
    freemat(&C);
    freemat(&C2);
    freemat(&I);
    freemat(&schur);
    freemat(&schur_co);
    freemat(&C3);
    free_poly(&ax);

    return 0;
}