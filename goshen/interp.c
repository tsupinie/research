
#include <stdio.h>
#include <math.h>

#define PI 3.1415926536

#ifndef DEG2RAD
#define DEG2RAD(x) (PI / 180. * x)
#endif

/*

PyArrayObject documentation:
If data a pointer to the object (type PyArrayObject*), then
    data->nd            is the number of dimensions (type int)
    data->dimensions    is an array of dimension lengths (type long int*, length data->nd)
    data->strides       is an array of strides (type long int*, length data->nd)
    data->data          is the actual data (type char* - cast to whatever type the data should be, for example float*)

To loop over the data using 3-dimenisonal indexes:

    data_values = (float*)(data->data);

    for (int kdz = 0; kdz < data->dimensions[0]; kdz++) {
        for (int jdy = 0; jdy < data->dimensions[1]; jdy++) {
            for (int idx = 0; idx < data->dimensions[2]; idx++) {
                index = kdz * data->dimensions[1] * data->dimensions[2] + jdy * data->dimensions[2] + idx;

                // do whatever with data_values[index] 
            }
        }
    }
 */

int _find_lb_index(float* x_data, int len_data, float x_target) {
    int idx_lbound = -1;

//  if (x_data[0] > 50000.)
//      printf("\n\nIn _find_lb_index ...\n");

    for (int idx = 0; idx < len_data - 1; idx++) {
        if ((x_data[idx] <= x_target && x_target < x_data[idx + 1]) || (x_data[idx] >= x_target && x_target > x_data[idx + 1]))
            idx_lbound = idx;

//      if (x_data[0] > 50000.)
//          printf("x_data[%d] = %f, idx_lbound = %d\n", idx, x_data[idx], idx_lbound);
    }

    return idx_lbound;
}

float _interp2d(float* x_data, float* y_data, float** u_data, int len_x_data, int len_y_data, float x_target, float y_target) {
   
    int idx_lbound = _find_lb_index(x_data, len_x_data, x_target);
    int idx_ubound;
    int jdy_lbound = _find_lb_index(y_data, len_y_data, y_target);
    int jdy_ubound;

//  printf("idx_lbound = %d, jdy_lbound = %d\n", idx_lbound, jdy_lbound);

    float numer, denom;
    int x_neg, y_neg;

    if (idx_lbound < 0) {
#ifdef INTERP_DEBUG
        printf("_interp2d(): x_target (%f) not in range of x_data (%f to %f).\n", x_target, x_data[0], x_data[len_x_data - 1]);
#endif
        return nan("");
    }

    if (jdy_lbound < 0) {
#ifdef INTERP_DEBUG
        printf("_interp2d(): y_target (%f) not in range of y_data (%f to %f).\n", y_target, y_data[0], y_data[len_y_data - 1]);
#endif
        return nan("");
    }

    idx_ubound = idx_lbound + 1;
    jdy_ubound = jdy_lbound + 1;

    numer = 0;
    for (int idx = 0; idx <= 1; idx++) {
        for (int jdy = 0; jdy <= 1; jdy++) {
            x_neg = (idx == 0 ? 1 : -1);
            y_neg = (jdy == 0 ? 1 : -1);
            numer += u_data[jdy_ubound + jdy][idx_lbound + idx] * x_neg * (x_data[idx_lbound + (1 - idx)] - x_target) * y_neg * (y_data[jdy_lbound + (1 - jdy)] - y_target);
        }
    }

    denom = (x_data[idx_ubound] - x_data[idx_lbound]) * (y_data[jdy_ubound] - y_data[jdy_lbound]);

    return numer / denom;
}

float _interp1d(float* x_data, float* u_data, int len_data, float x_target, bool wrap) {
    int idx_lbound = _find_lb_index(x_data, len_data, x_target);
    int idx_ubound = idx_lbound + 1;

    if (idx_lbound < 0) {
#ifdef INTERP_DEBUG
        printf("_interp1d(): x_target (%f) not in range of x_data (%f to %f).\n", x_target, x_data[0], x_data[len_data - 1]);
#endif
        if (!wrap) {
            return nan("");
        }
        else {
            if (x_data[0] < x_data[len_data - 1]) {
                if (x_target < x_data[0]) {
                    return u_data[0];
                }
                else if (x_target > x_data[len_data - 1]) {
                    return u_data[len_data - 1];
                }
            }
            else {
                if (x_target > x_data[0]) {
                    return u_data[0];
                }
                else if (x_target < x_data[len_data - 1]) {
                    return u_data[len_data - 1];
                }
            }
        }
    }

    return u_data[idx_lbound] + ((x_target - x_data[idx_lbound]) * u_data[idx_ubound] - (x_target - x_data[idx_lbound]) * u_data[idx_lbound]) / (x_data[idx_ubound] - x_data[idx_lbound]);
}

void interppt(PyArrayObject* data, PyArrayObject* z_axis, float height, PyArrayObject* y_axis, float y_distance, PyArrayObject* x_axis, float x_distance, float* interp_data, bool wrap) {
    float* data_values = NULL;
    float* z_values = NULL;
    float* y_values = NULL;
    float* x_values = NULL;

    float** tmp_data_slice = NULL;
    float** tmp_z_slice = NULL;
    float* tmp_data_column = NULL;
    float* tmp_z_column = NULL;

    int slice_index;

    data_values = (float*)data->data;
    z_values = (float*)z_axis->data;
    y_values = (float*)y_axis->data;
    x_values = (float*)x_axis->data;

    tmp_data_slice = _malloc_float_2d(data->dimensions[1], data->dimensions[2]);
    tmp_z_slice = _malloc_float_2d(data->dimensions[1], data->dimensions[2]);
    tmp_data_column = (float*)malloc(data->dimensions[0] * sizeof(float));
    tmp_z_column = (float*)malloc(data->dimensions[0] * sizeof(float));

/*
    printf("data->nd = %d\n", data->nd);

    for (int dim = 0; dim < data->nd; dim++) {
        printf("data->dimensions[%d] = %d\n", dim, data->dimensions[dim]);
    }
*/
    int num_vertical_levels = 0;

    for (int kdz = 0; kdz < data->dimensions[0]; kdz++) {
        for (int jdy = 0; jdy < data->dimensions[1]; jdy++) {
            for (int idx = 0; idx < data->dimensions[2]; idx++) {
                slice_index = _unravel_index(data->dimensions, kdz, jdy, idx);
                tmp_data_slice[jdy][idx] = data_values[slice_index];
                tmp_z_slice[jdy][idx] = z_values[slice_index];
            }
        }

        tmp_data_column[kdz] = _interp2d(x_values, y_values, tmp_data_slice, data->dimensions[2], data->dimensions[1], x_distance, y_distance);
        tmp_z_column[kdz] = _interp2d(x_values, y_values, tmp_z_slice, data->dimensions[2], data->dimensions[1], x_distance, y_distance);

//      printf("tmp_data_column[%d] = %f, tmp_z_column[%d] = %f\n", kdz, tmp_data_column[kdz], kdz, tmp_z_column[kdz]);

        num_vertical_levels++;

        if (kdz > 0 && ((tmp_z_column[kdz - 1] <= height && tmp_z_column[kdz] > height) || (tmp_z_column[kdz - 1] >= height && tmp_z_column[kdz] < height)))
            break;
    }

    *interp_data = _interp1d(tmp_z_column, tmp_data_column, num_vertical_levels, height, wrap);

//  printf("interp_data = %f\n", *interp_data);

    free(tmp_data_column);
    free(tmp_z_column);
    _free_float_2d(tmp_data_slice, data->dimensions[1], data->dimensions[2]);
    _free_float_2d(tmp_z_slice, data->dimensions[1], data->dimensions[2]);

    return;
}

void interppts(PyArrayObject* data, PyArrayObject* z_axis, PyArrayObject* heights, PyArrayObject* y_axis, PyArrayObject* y_distances, PyArrayObject* x_axis, PyArrayObject* x_distances, PyArrayObject* interp_data, bool wrap) {
    float* height_values = (float*)heights->data;
    float* y_dist_values = (float*)y_distances->data;
    float* x_dist_values = (float*)x_distances->data;
    float* interp_values = (float*)interp_data->data;

    for (int idx = 0; idx < heights->dimensions[0]; idx++) {
        interppt(data, z_axis, height_values[idx], y_axis, y_dist_values[idx], x_axis, x_dist_values[idx], &(interp_values[idx]), wrap);
    }

    return;
}

void interpsounding(PyArrayObject* data, PyArrayObject* y_axis, float y_distance, PyArrayObject* x_axis, float x_distance, PyArrayObject* interp_data) {
    float* data_values = (float*)data->data;
    float* x_values = (float*)x_axis->data;
    float* y_values = (float*)y_axis->data;
    float* interp_values = (float*)interp_data->data;

    float** tmp_data_slice = _malloc_float_2d(data->dimensions[1], data->dimensions[2]);

    int slice_index;

#ifdef SND_INTERP_DEBUG
    printf("data->nd = %d\n", data->nd);

    for (int dim = 0; dim < data->nd; dim++) {
        printf("data->dimensions[%d] = %d\n", dim, data->dimensions[dim]);
    }

    printf("x_distance = %f, y_distance = %f\n", x_distance, y_distance);
#endif

    for (int kdz = 0; kdz < data->dimensions[0]; kdz++) {
        for (int jdy = 0; jdy < data->dimensions[1]; jdy++) {
            for (int idx = 0; idx < data->dimensions[2]; idx++) {
                slice_index = _unravel_index(data->dimensions, kdz, jdy, idx);
                tmp_data_slice[jdy][idx] = data_values[slice_index];
            }
        }

        interp_values[kdz] = _interp2d(x_values, y_values, tmp_data_slice, data->dimensions[2], data->dimensions[1], x_distance, y_distance);
    }

    _free_float_2d(tmp_data_slice, data->dimensions[1], data->dimensions[2]);

    return;
}

void interpz(PyArrayObject* data, PyArrayObject* z_axis, float height, PyArrayObject* interp_data, bool wrap) {
    float* data_values = NULL;
    float* z_values = NULL;
    float* interp_values = NULL;

    float* temp_data_array = NULL;
    float* temp_z_array = NULL;
    int interp_index, data_index;

    if (data->nd != 3) {
        printf("Number of dimensions in the input data should be 3 (is %d).\n", data->nd);
        return;
    }

    if (z_axis->nd != 3) {
        printf("Number of dimensions in the z-axis should be 3 (is %d).\n", z_axis->nd);
        return;
    }

    if (interp_data->nd != 2) {
        printf("Number of dimensions in the output array should be 2 (is %d).\n", interp_data->nd);
        return;
    }

    for (int idx = 0; idx < data->nd; idx++) {
        if (data->dimensions[idx] != z_axis->dimensions[idx]) {
            printf("Dimensions of the input data and the z-axis are different.\n");
            return;
        }

        if (idx > 0 && data->dimensions[idx] != interp_data->dimensions[idx - 1]) {
            printf("Dimensions of the output array don't match the input data.\n");
            return;
        }
    }

    data_values = (float*)(data->data);
    z_values = (float*)(z_axis->data);
    interp_values = (float*)(interp_data->data);

    temp_data_array = (float*)malloc(data->dimensions[0] * sizeof(float));
    temp_z_array = (float*)malloc(z_axis->dimensions[0] * sizeof(float));

    for (int jdy = 0; jdy < data->dimensions[1]; jdy++) {
        for (int idx = 0; idx < data->dimensions[2]; idx++) {
            interp_index = _unravel_index(interp_data->dimensions, jdy, idx);

            for (int kdz = 0; kdz < data->dimensions[0]; kdz++) {
                data_index = _unravel_index(data->dimensions, kdz, jdy, idx);
                temp_data_array[kdz] = data_values[data_index];
                temp_z_array[kdz] = z_values[data_index];
            }

            interp_values[interp_index] = _interp1d(temp_z_array, temp_data_array, data->dimensions[0], height, wrap);
        }
    }
   
    free(temp_data_array);
    free(temp_z_array);

    temp_data_array = NULL;
    temp_z_array = NULL;

//  interp_data->data = (char*)interp_values;

    return;
}

float _cone_height(float distance, float base_height, float elev_angle) {
    const double earth_radius = 6371000.;
    const double eff_factor = 4. / 3.;
    const double eff_earth_radius = earth_radius * eff_factor;

    elev_angle = DEG2RAD(elev_angle);
   
    double z = (cos(elev_angle) / cos(elev_angle + distance / eff_earth_radius) - 1) * eff_earth_radius + base_height;
    return float(z);
}

void interpcone(PyArrayObject* data, PyArrayObject* z_axis, PyArrayObject* y_axis, PyArrayObject* x_axis, float elev_angle, float base_z, float base_y, float base_x, PyArrayObject* interp_data, PyArrayObject* cone_height, bool wrap) {
    float* z_values = (float*)z_axis->data;
    float* y_values = (float*)y_axis->data;
    float* x_values = (float*)x_axis->data;

    float* interp_values = (float*)interp_data->data;
    float* data_values = (float*)data->data;
    float* cone_height_values = (float*)cone_height->data;

    float* temp_data_array = (float*)malloc(data->dimensions[0] * sizeof(float));
    float* temp_z_array = (float*)malloc(data->dimensions[0] * sizeof(float));

    int interp_index = 0;
    int data_index = 0;
    int num_vertical_levels;

    float interp_height = 0;
    float gp_distance;

    for (int jdy = 0; jdy < data->dimensions[1]; jdy++) {
        for (int idx = 0; idx < data->dimensions[2]; idx++) {
            interp_index = _unravel_index(interp_data->dimensions, jdy, idx);

            gp_distance = sqrt(pow(y_values[jdy] - base_y, 2) + pow(x_values[idx] - base_x, 2));
            interp_height = _cone_height(gp_distance, base_z, elev_angle);

            cone_height_values[interp_index] = interp_height;

            num_vertical_levels = 0;

            for (int kdz = 0; kdz < data->dimensions[0]; kdz++) {
                data_index = _unravel_index(data->dimensions, kdz, jdy, idx);
                temp_data_array[kdz] = data_values[data_index];
                temp_z_array[kdz] = z_values[data_index];

                num_vertical_levels++;

                if (kdz > 0 && temp_z_array[kdz - 1] <= interp_height && temp_z_array[kdz] > interp_height)
                    break;
            }

            interp_values[interp_index] = _interp1d(temp_z_array, temp_data_array, num_vertical_levels, interp_height, wrap);
        }
    }

    free(temp_data_array);
    free(temp_z_array);

    temp_data_array = NULL;
    temp_z_array = NULL;

    return;
}

void interpheights(PyArrayObject* data, PyArrayObject* z_axis, PyArrayObject* heights, PyArrayObject* interp_data, bool wrap) {
    float* z_values = (float*)z_axis->data;

    float* interp_values = (float*)interp_data->data;
    float* data_values = (float*)data->data;
    float* height_values = (float*)heights->data;

    float* temp_data_array = (float*)malloc(data->dimensions[0] * sizeof(float));
    float* temp_z_array = (float*)malloc(data->dimensions[0] * sizeof(float));

    int interp_index = 0;
    int data_index = 0;
    int num_vertical_levels;

    for (int jdy = 0; jdy < data->dimensions[1]; jdy++) {
        for (int idx = 0; idx < data->dimensions[2]; idx++) {
            interp_index = _unravel_index(interp_data->dimensions, jdy, idx);

            num_vertical_levels = 0;

            for (int kdz = 0; kdz < data->dimensions[0]; kdz++) {
                data_index = _unravel_index(data->dimensions, kdz, jdy, idx);
                temp_data_array[kdz] = data_values[data_index];
                temp_z_array[kdz] = z_values[data_index];

                num_vertical_levels++;

                if (kdz > 0 && temp_z_array[kdz - 1] <= height_values[interp_index] && temp_z_array[kdz] > height_values[interp_index])
                    break;
            }

            interp_values[interp_index] = _interp1d(temp_z_array, temp_data_array, num_vertical_levels, height_values[interp_index], wrap);
        }
    }

    free(temp_data_array);
    free(temp_z_array);

    temp_data_array = NULL;
    temp_z_array = NULL;

}
