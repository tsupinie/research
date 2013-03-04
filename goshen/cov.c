
#include <stdio.h>
#include <math.h>

float _mean(float* array, int array_len) {
    float array_sum = 0;

    for (int idx = 0; idx < array_len; idx++) {
        array_sum += array[idx];
    }

    return array_sum / array_len;
}

float _covariance(float* array1, float* array2, int array_len, float* means) {
    float cov_sum = 0;
    float array1_mean = 0;
    float array2_mean = 0;

    if (means == NULL) {
        array1_mean = _mean(array1, array_len);
        array2_mean = _mean(array2, array_len);
    }
    else {
        array1_mean = means[0];
        array2_mean = means[1];
    }

#ifdef COV_DEBUG
    printf("array1_mean = %f, array2_mean = %f\n", array1_mean, array2_mean);
#endif

    for (int idx = 0; idx < array_len; idx++) {
        cov_sum += (array1[idx] - array1_mean) * (array2[idx] - array2_mean);
    }

#ifdef COV_DEBUG
    printf("covariance = %f\n", cov_sum / (array_len - 1));
#endif

    return cov_sum / (array_len - 1);
}

float _correlation(float* array1, float* array2, int array_len, float* means) {
    float var1 = _covariance(array1, array1, array_len, means);
    float var2 = _covariance(array2, array2, array_len, means);

    return _covariance(array1, array2, array_len, means) / (sqrt(var1) * sqrt(var2));
}

void ens_covariance(PyArrayObject* ens_data, PyArrayObject* ens_obs, PyArrayObject* cov_data, bool normalize) {
    int grid_index;
    int ens_index;

    float* ens_data_values = (float*)ens_data->data;
    float* ens_obs_values = (float*)ens_obs->data;
    float* cov_data_values = (float*)cov_data->data;

    float* tmp_ens_data = (float*)malloc(ens_data->dimensions[0] * sizeof(float));
    float* tmp_means = (float*)malloc(2 * sizeof(float));

    tmp_means[1] = _mean(ens_obs_values, ens_data->dimensions[0]);

#ifdef COV_DEBUG_STUFF
    printf("ens_data->nd = %d\n", ens_data->nd);

    for (int idx = 0; idx < ens_data->nd; idx++) {
        printf("ens_data->dimensions[%d] = %d\n", idx, ens_data->dimensions[idx]);
    }
#endif

    for (int kdz = 0; kdz < ens_data->dimensions[1]; kdz++) {
        for (int jdy = 0; jdy < ens_data->dimensions[2]; jdy++) {
            for (int idx = 0; idx < ens_data->dimensions[3]; idx++) {
                grid_index = _unravel_index(cov_data->dimensions, kdz, jdy, idx);

#ifdef COV_DEBUG_STUFF
                printf("grid_index = %d\n", grid_index);
#endif

                for (int lde = 0; lde < ens_data->dimensions[0]; lde++) {
                    ens_index = _unravel_index(ens_data->dimensions, lde, kdz, jdy, idx);

                    tmp_ens_data[lde] = ens_data_values[ens_index];
#ifdef COV_DEBUG
                    printf("tmp_ens_data[%d] = %f\n", lde, tmp_ens_data[lde]);
#endif
                }

                tmp_means[0] = _mean(tmp_ens_data, ens_data->dimensions[0]);
                if (normalize) {
                    cov_data_values[grid_index] = _correlation(tmp_ens_data, ens_obs_values, ens_data->dimensions[0], tmp_means);
                }
                else {
                    cov_data_values[grid_index] = _covariance(tmp_ens_data, ens_obs_values, ens_data->dimensions[0], tmp_means);
                }
            }
        }
    }

    free(tmp_ens_data);
    tmp_ens_data = NULL;

    free(tmp_means);
    tmp_means = NULL;
}
