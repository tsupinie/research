
float** _malloc_float_2d(int x_length, int y_length) {
    float** array = (float**)malloc(x_length * sizeof(float*));
    for (int idx = 0; idx < x_length; idx++) {
        array[idx] = (float*)malloc(y_length * sizeof(float));
    }

    return array;
}

void _free_float_2d(float** array, int x_length, int y_length) {
    for (int idx = 0; idx < x_length; idx++) {
        free(array[idx]);
    }
    free(array);

    return;
}

int _unravel_index(long* dims, int wdt, int kdz, int jdy, int idx) {
    return wdt * dims[1] * dims[2] * dims[3] + kdz * dims[2] * dims[3] + jdy * dims[3] + idx;
}

int _unravel_index(long* dims, int kdz, int jdy, int idx) {
    return kdz * dims[1] * dims[2] + jdy * dims[2] + idx;
}

int _unravel_index(long* dims, int jdy, int idx) {
    return jdy * dims[1] + idx;
}
