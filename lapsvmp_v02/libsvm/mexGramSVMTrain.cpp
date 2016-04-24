/*
Modified and updated for the latest libsvm version by Stefano Melacci (MAY 2009), mela@dii.unisi.it
*/
extern "C" {
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "svm.h"
#include "mex.h"
    
#define Malloc(type,n) (type *)malloc((n)*sizeof(type))
    
/* Matlab Input Arguments */
#define SAMPLES prhs[0]
#define LABELS  prhs[1]
#define PARAMETERS prhs[2]

/* Matlab Output Arguments */
#define ALPHAY plhs[0]
#define SVs    plhs[1]
#define BIAS   plhs[2]
#define NSV        plhs[3]
#define NLABEL     plhs[4]
   
static void parse_command_line(double params[]);
static void read_problem(double examples[], double labels[],size_t dim, size_t num);
void svm_copy_model(double[],double[],double[],double[],double[],svm_model*);
       
struct svm_parameter param;		// set by parse_command_line
struct svm_problem prob;		// set by read_problem
struct svm_model *model;
struct svm_node *x_space;
    
 
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[] )
{
    double *parameters;
    double *labels;
    double *samples;
    size_t dimension;
    size_t numberofexamples;
    int l, n;
    
    parameters = mxGetPr(PARAMETERS);
    parse_command_line(parameters);
    
    dimension = mxGetN(SAMPLES);
    numberofexamples = mxGetM(SAMPLES);
    
    samples = mxGetPr(SAMPLES);
    labels  = mxGetPr(LABELS);
    
    read_problem(samples,labels,dimension,numberofexamples);
    
    model = svm_train(&prob,&param);
    
    l=model->l;
    n=svm_get_nr_class(model);
    
    ALPHAY = mxCreateDoubleMatrix(n-1, l, mxREAL);
    SVs    = mxCreateDoubleMatrix( 1,l, mxREAL);
    BIAS   = mxCreateDoubleMatrix( 1,n*(n-1)/2, mxREAL);
    NSV    =  mxCreateDoubleMatrix( 1,n, mxREAL);
    NLABEL = mxCreateDoubleMatrix( 1,n, mxREAL);
    
    svm_copy_model(mxGetPr(ALPHAY),mxGetPr(SVs),mxGetPr(BIAS),mxGetPr(NSV),mxGetPr(NLABEL),model);
    
    svm_destroy_model(model);
    free(prob.y);
    free(prob.x);
    free(x_space);
}
}


static void parse_command_line(double params[])
{
    // default values
    param.svm_type = (int) params[7];
    param.kernel_type = PRECOMPUTED;
    param.degree = params[1];
    param.gamma = params[2];	// 1/k
    param.coef0 = params[3];
    param.nu = params[8];
    param.cache_size = params[5];
    param.C = params[4];
    param.eps = params[6];
    param.p = params[9];
    param.shrinking = (int) params[10];
    param.nr_weight = 0;
    param.weight_label = NULL;
    param.weight = NULL;
}


// read in a problem (in svmlight format)
static void read_problem(double examples[],double labels[],size_t dim,size_t num)
{
    size_t elements;
    int i,j,p,t;
    
    prob.l = num;
    prob.y = Malloc(double,prob.l);
    prob.x = Malloc(struct svm_node *,prob.l);
    elements=dim*num;
    
    x_space = Malloc(struct svm_node,prob.l + elements);

    j=0;
    p=0;
    
    for(i=0;i<prob.l;i++)
    {
        double label;
        prob.x[i] = &x_space[j];
        prob.y[i] = labels[i];
        
        x_space[j].index = 1;
        x_space[j].value = i+1;
        j++;
        for(t=0;t<dim;t++)
        {
            x_space[j].index=t+2;
            x_space[j].value=examples[p];
            j++;
            p++;
        }  
    }
}


void svm_copy_model(double alphay[],double svs[],double bias[],double nsv[],double nlabel[],svm_model* model)
{
    int i,t;
    int n,l;
    n=model->nr_class;
    l=model->l;
    
    for(i=0;i<n;i++)
    {
        nsv[i]=model->nSV[i];
        nlabel[i]=model->label[i];
    }
    
    for(i=0;i<n*(n-1)/2;i++)
    {
        bias[i]=model->rho[i];
    }
    
    for(i=0;i<l;i++)
    {
        svs[i]=(model->SV[i])->value; // I removed the "- 1" and fixed matlab files;
    }
    
    t=0;
    for(int i=0;i<l;i++)
    {
        for(int j=0;j<n-1;j++)
        {
            alphay[t]=model->sv_coef[j][i];
            t++;
        }
    }
}
