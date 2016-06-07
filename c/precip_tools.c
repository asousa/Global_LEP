#include <Python.h>
#include <stdio.h>
#include <math.h>
#include "time.h"
#include <unistd.h>
# include <stdlib.h>
#include <stdbool.h>

#include <numpy/arrayobject.h>

// You can use it in code like so:

// DPRINTF("key == %d\n", key);


// // WORKING HELLO WORLD
// static PyObject* helloworld(PyObject* self)
// {
//     return Py_BuildValue("s", "Hello, Python extensions!!");
// }

// static char helloworld_docs[] =
//     "helloworld( ): Any message you want to put here!!\n";

// static PyMethodDef helloworld_funcs[] = {
//     {"helloworld", (PyCFunction)helloworld, 
//      METH_NOARGS, helloworld_docs},
//     {NULL}
// };

// void inithelloworld(void)
// {
//     Py_InitModule3("helloworld", helloworld_funcs,
//                    "Extension module example!");
// }


#define D2R  (3.14159/180.0)
#define H_IONO 100000.0
#define R_E 6370000.0
#define I0_sim -100000.0


// Prototypes ----------------------------------------------
static double scale_factor_single(double inp_lat, double inp_lon, double out_lat, double out_lon, double I0);


// // Totally working version (but it's slow!)  ----------
// -------------------------------------------------------
// static PyObject* example (PyObject *self, PyObject *args) {
//     PyObject *arg1=NULL, *arg2=NULL, *arg3=NULL;
//     PyArrayObject *arr1=NULL, *arr2=NULL, *oarr=NULL;
//     int nd;
//     int i,j;
//     double inp_lat, inp_lon, I0;
//     double * lat_ptr;
//     double * lon_ptr;
//     double * optr;
//     double tmp;
//     if (!PyArg_ParseTuple(args, "ddOOOd", &inp_lat, &inp_lon, &arg1, &arg2, &arg3, &I0)) {
//       return NULL;
//     }

//     // In lats
//     arr1 = (PyArrayObject*)PyArray_FROM_OTF(arg1, NPY_DOUBLE, NPY_IN_ARRAY);
//     if (arr1 == NULL) return NULL;

//     // In lons
//     arr2 = (PyArrayObject*)PyArray_FROM_OTF(arg2, NPY_DOUBLE, NPY_IN_ARRAY);
//     if (arr2 == NULL) return NULL;

//     // output space
//     oarr = (PyArrayObject*)PyArray_FROM_OTF(arg3, NPY_DOUBLE, NPY_INOUT_ARRAY);
//     if (oarr == NULL) return NULL;

//     /*vv* code that makes use of arguments *vv*/

//     nd = PyArray_NDIM(arr1);   //number of dimensions
//     npy_intp *lat_shape = PyArray_DIMS(arr1);  // npy_intp array of length nd showing length in each dim.
//     npy_intp *lon_shape = PyArray_DIMS(arr2);  // npy_intp array of length nd showing length in each dim.
    
//     lat_ptr = (double *)PyArray_DATA(arr1);
//     lon_ptr = (double *)PyArray_DATA(arr2);

//     // for (i=0; i<nd; ++i) {
//     //     printf("shape: %ld\n",lat_shape[i]);
//     // }
    
//     for (i=0; i < lat_shape[0]; i++) {
//       printf("%d: ",i);
//       for (j=0; j < lon_shape[0]; j++) {


//         // printf("(%g, %g)",lat_ptr[i], lon_ptr[j]);
//         optr = (double *) PyArray_GETPTR2(oarr, i, j);

//        // printf("sent: %g %g %g %g %g\n",inp_lat, inp_lon, lat_ptr[i], lon_ptr[j], I0);

//         tmp = scale_factor_single(inp_lat, inp_lon, lat_ptr[i], lon_ptr[j], I0);
//         *optr = tmp;
//         // *optr = lon_ptr[j]*lat_ptr[i];
//         // printf("%f ",tmp);

//       }
//       // printf("\n");
//     }

//     /*^^* code that makes use of arguments *^^*/

//     // Tidy up
//     Py_DECREF(arr1);
//     Py_DECREF(arr2);
//     Py_DECREF(oarr);
//     Py_INCREF(Py_None);
//     return Py_None;
// }
// -------------------------------------------------------------


/*  wrapped cosine function */
static PyObject* cos_func_np(PyObject* self, PyObject* args)
{

    PyArrayObject *lat_in;
    pyArrayObject *lon_in;
    PyObject      *out_array;
    NpyIter *in_iter;
    NpyIter *out_iter;
    NpyIter_IterNextFunc *in_iternext;
    NpyIter_IterNextFunc *out_iternext;

    /*  parse single numpy array argument */
    if (!PyArg_ParseTuple(args, "OO",  &lat_in, &lon_in))
        return NULL;

    /*  construct the output array, like the input array */
    out_array = PyArray_NewLikeArray(in_array, NPY_ANYORDER, NULL, 0);
    if (out_array == NULL)
        return NULL;

    /*  create the iterators */
    in_iter = NpyIter_New(in_array, NPY_ITER_READONLY, NPY_KEEPORDER,
                             NPY_NO_CASTING, NULL);
    if (in_iter == NULL)
        goto fail;

    out_iter = NpyIter_New((PyArrayObject *)out_array, NPY_ITER_READWRITE,
                          NPY_KEEPORDER, NPY_NO_CASTING, NULL);
    if (out_iter == NULL) {
        NpyIter_Deallocate(in_iter);
        goto fail;
    }

    in_iternext = NpyIter_GetIterNext(in_iter, NULL);
    out_iternext = NpyIter_GetIterNext(out_iter, NULL);
    if (in_iternext == NULL || out_iternext == NULL) {
        NpyIter_Deallocate(in_iter);
        NpyIter_Deallocate(out_iter);
        goto fail;
    }
    double ** in_dataptr  = (double **) NpyIter_GetDataPtrArray(in_iter);
    double ** out_dataptr = (double **) NpyIter_GetDataPtrArray(out_iter);

    /*  iterate over the arrays */
    do {
        **out_dataptr = cos(**in_dataptr);
    } while(in_iternext(in_iter) && out_iternext(out_iter));





    /*  clean up and return the result */
    NpyIter_Deallocate(in_iter);
    NpyIter_Deallocate(out_iter);
    Py_INCREF(out_array);
    return out_array;

    /*  in case bad things happen */
    fail:
        Py_XDECREF(out_array);
        return NULL;
}





static PyObject *scale_factor(PyObject *self, PyObject *args) {
  // Calculate the scaling factor at a point due to flash.
  // (single weight)
    double inp_lat, inp_lon, out_lat, out_lon;
    double I0;
    double ratio;
    // double dlong_sim = 0.7;
    // // double dlat, dlong, clat1, clat2, slat1, slat2;
    // double dist_iono_0, dist_iono_1;
    // double dist_lat_0;
    // double dist_long_0;
    // double R_0, R_1;
    // double xi_0, xi_1;
    // double a, b;
    // double new_weight, old_weight, ratio;
    // double dlat, dlong;

    // gather inputs from Python
   if (!PyArg_ParseTuple(args, "ddddd", &inp_lat, &inp_lon, &out_lat, &out_lon, &I0)) {
      return NULL;
   }
   
   ratio = scale_factor_single(inp_lat, inp_lon, out_lat, out_lon, I0);

  //  // printf("received: %g %g %g %g %g\n",inp_lat, inp_lon, out_lat, out_lon, I0);


  //  dlat = fabs(D2R*(out_lat - inp_lat));
  //  dlong= fabs(D2R*(out_lon - inp_lon));


  // // # ------------- More-realistic attempt, using same scaling factor as latitude.
  // // # Ratio of (2.1) --> equation (5.5) in Jacob's thesis
  // dist_lat_0 = (R_E + H_IONO/2.0) * dlat;
  // dist_long_0= (R_E + H_IONO/2.0) * dlong_sim * D2R;
  // dist_iono_0  = hypot(dist_lat_0, dist_long_0);

  // R_0 =  hypot(dist_iono_0, H_IONO);
  // xi_0 = atan2(dist_iono_0, H_IONO);

  // // # Compute current (latitude and longitude dependent) weighting:
  // // # (Use Haversine formula since we're moving further out than reasonable for Cartesian)
  // // # ----> Does not wrap around the north / south poles! What gives?
  // a = pow(sin(D2R*(out_lat - inp_lat)/2.0),2);
  // b = cos(D2R*inp_lat)*cos(D2R*out_lat)*pow(sin(D2R*(out_lon - inp_lon)/2.0),2);
  // dist_iono_1 = 2.0*R_E*asin(sqrt(a + b));
  // R_1 = hypot(dist_iono_1, H_IONO);
  // xi_1= atan2(dist_iono_1, H_IONO);


  // new_weight = sin(xi_1)/R_1;
  // old_weight = sin(xi_0);

  // ratio = new_weight/old_weight;
  // ratio = fabs(I0/I0_sim)*ratio;

    // return Py_BuildValue("f", scalefactor);
  return Py_BuildValue("f", ratio);
}



static double scale_factor_single(double inp_lat, double inp_lon, double out_lat, double out_lon, double I0) {
  // double inp_lat, inp_lon, out_lat, out_lon;
  // double I0;
  double dlong_sim = 0.7;
  // double dlat, dlong, clat1, clat2, slat1, slat2;
  double dist_iono_0, dist_iono_1;
  double dist_lat_0;
  double dist_long_0;
  double R_0, R_1;
  double xi_0, xi_1;
  double a, b;
  double new_weight, old_weight, ratio;
  double dlat, dlong;


   // printf("received: %g %g %g %g %g\n",inp_lat, inp_lon, out_lat, out_lon, I0);

  dlat = fabs(D2R*(out_lat - inp_lat));
  dlong= fabs(D2R*(out_lon - inp_lon));


  // # ------------- More-realistic attempt, using same scaling factor as latitude.
  // # Ratio of (2.1) --> equation (5.5) in Jacob's thesis
  dist_lat_0 = (R_E + H_IONO/2.0) * dlat;
  dist_long_0= (R_E + H_IONO/2.0) * dlong_sim * D2R;
  dist_iono_0  = hypot(dist_lat_0, dist_long_0);

  R_0 =  hypot(dist_iono_0, H_IONO);
  xi_0 = atan2(dist_iono_0, H_IONO);

  // # Compute current (latitude and longitude dependent) weighting:
  // # (Use Haversine formula since we're moving further out than reasonable for Cartesian)
  // # ----> Does not wrap around the north / south poles! What gives?
  a = pow(sin(D2R*(out_lat - inp_lat)/2.0),2);
  b = cos(D2R*inp_lat)*cos(D2R*out_lat)*pow(sin(D2R*(out_lon - inp_lon)/2.0),2);
  dist_iono_1 = 2.0*R_E*asin(sqrt(a + b));
  R_1 = hypot(dist_iono_1, H_IONO);
  xi_1= atan2(dist_iono_1, H_IONO);


  new_weight = sin(xi_1)/R_1;
  old_weight = sin(xi_0);

  ratio = new_weight/old_weight;
  ratio = fabs(I0/I0_sim)*ratio;

  return ratio;
}



// Interconnection stuff ------------------------------------------------
static char func_docs[] = "func( ): Any message you want to put here!!\n";

static PyObject* say_hi(PyObject* self)
{
    return Py_BuildValue("s", "Hi from precip_tools!");
}


static PyMethodDef module_funcs[] = {
    {"say_hi", (PyCFunction)say_hi, METH_NOARGS, func_docs},
    { "scale_factor", scale_factor, METH_VARARGS, NULL },
    // { "precip_single_flash", precip_single_flash, METH_VARARGS, NULL },
    { "example", example, METH_VARARGS, NULL },
    { NULL, NULL, 0, NULL }
};

void initprecip_tools(void)
{
    Py_InitModule3("precip_tools", module_funcs,
                   "Extension module example!");
    import_array();
}
// End iterconnection stuff --------------------------------------------



