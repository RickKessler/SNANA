#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <Python.h>

// ========= GLOBAL DECLARATIONS =================

// -----------------------------------------------
//        skewNormal Global Declarations
// -----------------------------------------------

//#define  IGNORE_SNANA_PYTHON   
#define  SNANA_PYTHON   

#define  INITDONE_SKEWNORMAL  8888   // to check that init was called.

int INITFLAG_SKEWNORMAL = 0;

/* xxxx mark delete; moved elsewhere xxx
// function prototypesx
void  init_skewNormal(int seed); // one-time init
double skewNormalRan(int seed, double loc, double scale, double skew) ;
xxxxxxx */

// =======================================================
//          Start Code
// =======================================================

/*
void  Initialize ()
{
    Py_Initialize();
}

void Finalize ()
{
    Py_Finalize();
}
*/

void init_skewNormal(int seed) {


  //E. Jennings. dummy call to RNG to init and set random seed in python in sync with SNANA
  double loc = -0.1;
  double scale = 0.3;
  double skew = 5.0;
  skewNormalRan(seed,loc,scale,skew);

  INITFLAG_SKEWNORMAL = INITDONE_SKEWNORMAL  ;

  return ;

}  //end init_skewNormal



double skewNormalRan(int seed, double loc, double scale, double skew) {

  // Created Sep 2016 E. Jennings 
  // Return random draw from skewNormal defined by loc,scale,skew parameters

  //  char fnam[] =  "skewNormalRan" ;
  
  // ---------------BEGIN -------------

#ifndef SNANA_PYTHON
  sprintf(c1err,"Cannot use skewNormal option because");
  sprintf(c2err,"Env SNANA_PYTHON_DIR is not defined. See install guide");
  errmsg(SEV_FATAL, 0, fnam, c1err, c2err);     
#endif

/*
  PyObject *pName, *pModule, *pDict, *pFunc, *pValue, *pArgs, *rtn_rvs, *obj_rtn;
  Py_ssize_t j=0;
  char str1[1024] = "py_func";
  char str2[1024] = "return_random";
  float rvs;
  int num_in=5; //number of params passed to python func
  int N = 1; //num of rvs to draw
  double flt_rtn[N];
  double  params_of_pdf[num_in];

  params_of_pdf[0] = loc; 
  params_of_pdf[1] = scale;
  params_of_pdf[2] = skew;
  params_of_pdf[3] = (float)N;
  if( INITFLAG_SKEWNORMAL ==0){
	params_of_pdf[4] = (float)RNG_SEED;
  }
  else{ params_of_pdf[4] = (float)0.0;}

  if (INITFLAG_SKEWNORMAL == 0)
        {
                call_num++;
                Initialize();
                atexit(Finalize);
        }

  rtn_rvs = PyTuple_New(N);

  PyObject *sys = PyImport_ImportModule("sys");
  PyObject *path = PyObject_GetAttrString(sys, "path");
  PyList_Append(path, PyString_FromString("."));

  pName = PyString_FromString(str1);
  pModule = PyImport_Import(pName);
  pDict = PyModule_GetDict(pModule);

  pFunc = PyDict_GetItemString(pDict, str2);

  pArgs = PyTuple_New(num_in);

  for (i = 0; i < num_in; i++)
  {
	pValue = PyFloat_FromDouble(params_of_pdf[i]);
	if (!pValue)
	{
		PyErr_Print();
		return 1;
	}
	PyTuple_SetItem(pArgs, i, pValue);
  }

  if (PyCallable_Check(pFunc))
  {
	rtn_rvs=PyObject_CallObject(pFunc, pArgs);
  }
	else
  {
	PyErr_Print();
  }

  for (j=0;j<N;j++)
  {
	obj_rtn = PyTuple_GetItem(rtn_rvs,j);
	flt_rtn[j] = PyFloat_AsDouble(obj_rtn);
  }

  if (rtn_rvs != NULL)
  {
	for (i=0;i<N;i++)
	{
		printf("Random variable received from python : %lf\n", flt_rtn[i]);
	}
  }

  Py_DECREF(pModule);
  Py_DECREF(pDict);
  Py_DECREF(pName);
  Py_DECREF(pArgs);
  Py_DECREF(pValue);
  Py_DECREF(rtn_rvs);

  return (flt_rtn[0]);
*/
  return -0.1;

} // end skewNormalRan


 
