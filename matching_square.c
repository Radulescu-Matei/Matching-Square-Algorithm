// Author: APD team, except where source was noted

#include "helpers.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>

#define CONTOUR_CONFIG_COUNT    16
#define FILENAME_MAX_SIZE       50
#define STEP                    8
#define SIGMA                   200
#define RESCALE_X               2048
#define RESCALE_Y               2048

#define CLAMP(v, min, max) if(v < min) { v = min; } else if(v > max) { v = max; }

// Structure used to pass all of the needed arguments to the function that is
// parallelized.
typedef struct {
    ppm_image *image;
    ppm_image **contour_map;
    ppm_image *scaled_image;
    unsigned char **grid;
    int thread_id;
    int nr_threads;
    // Flag used to determine wether the image needs to be rescaled or not.
    int rescale_flag;
    pthread_barrier_t *barrier;
} arg_images;


// Updates a particular section of an image with the corresponding contour pixels.
// Used to create the complete contour image.
void update_image(ppm_image *image, ppm_image *contour, int x, int y) {
    for (int i = 0; i < contour->x; i++) {
        for (int j = 0; j < contour->y; j++) {
            int contour_pixel_index = contour->x * i + j;
            int image_pixel_index = (x + i) * image->y + y + j;

            image->data[image_pixel_index].red = contour->data[contour_pixel_index].red;
            image->data[image_pixel_index].green = contour->data[contour_pixel_index].green;
            image->data[image_pixel_index].blue = contour->data[contour_pixel_index].blue;
        }
    }
}

// Corresponds to step 1 of the marching squares algorithm, which focuses on sampling the image.
// Builds a p x q grid of points with values which can be either 0 or 1, depending on how the
// pixel values compare to the `sigma` reference value. The points are taken at equal distances
// in the original image, based on the `step_x` and `step_y` arguments.
void build_grid(arg_images *args) {
    int p = (*args).scaled_image->x / STEP;
    int q = (*args).scaled_image->y / STEP;

    int start = (*args).thread_id * (double)p / (*args).nr_threads;
    int end;
    
    // Devides the steps of the loop between the threads used.
    if(((*args).thread_id + 1) * (double)p / (*args).nr_threads < p){
        end = ((*args).thread_id + 1) * (double)p / (*args).nr_threads;
    }else end = p;

    for (int i = start; i < end; i++) {
        for (int j = 0; j < q; j++) {
            ppm_pixel curr_pixel = (*args).scaled_image->data[i * STEP * (*args).scaled_image->y + j * STEP];

            unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

            if (curr_color > SIGMA) {
                (*args).grid[i][j] = 0;
            } else {
                (*args).grid[i][j] = 1;
            }
        }
    }
    
    (*args).grid[p][q] = 0;

    // last sample points have no neighbors below / to the right, so we use pixels on the
    // last row / column of the input image for them
    for (int i = start; i < end; i++) {
        ppm_pixel curr_pixel = (*args).scaled_image->data[i * STEP * (*args).scaled_image->y + (*args).scaled_image->x - 1];

        unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

        if (curr_color > SIGMA) {
            (*args).grid[i][q] = 0;
        } else {
            (*args).grid[i][q] = 1;
        }
    }

    start = (*args).thread_id * (double)q / (*args).nr_threads;
    
    if(((*args).thread_id + 1) * (double)q / (*args).nr_threads < p){
        end = ((*args).thread_id + 1) * (double)q / (*args).nr_threads;
    }else end = q;
    for (int j = start; j < end; j++) {
        ppm_pixel curr_pixel = (*args).scaled_image->data[((*args).scaled_image->x - 1) * (*args).scaled_image->y + j * STEP];

        unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

        if (curr_color > SIGMA) {
            (*args).grid[p][j] = 0;
        } else {
            (*args).grid[p][j] = 1;
        }
    }
}

// Corresponds to step 2 of the marching squares algorithm, which focuses on identifying the
// type of contour which corresponds to each subgrid. It determines the binary value of each
// sample fragment of the original image and replaces the pixels in the original image with
// the pixels of the corresponding contour image accordingly.
void march(arg_images *args) {
    int p = (*args).scaled_image->x / STEP;
    int q = (*args).scaled_image->y / STEP;
    int start = (*args).thread_id * (double)p / (*args).nr_threads;
    int end;
    
    // Devides the steps of the loop between the threads used
    if(((*args).thread_id + 1) * (double)p / (*args).nr_threads < p){
        end = ((*args).thread_id + 1) * (double)p / (*args).nr_threads;
    }else end = p;

    for (int i = start; i < end; i++) {
        for (int j = 0; j < q; j++) {
            unsigned char k = 8 * (*args).grid[i][j] + 4 * (*args).grid[i][j + 1] + 2 * (*args).grid[i + 1][j + 1] + 1 * (*args).grid[i + 1][j];
            update_image((*args).scaled_image, (*args).contour_map[k], i * STEP, j * STEP);
        }
    }
}

// Calls `free` method on the utilized resources.
void free_resources(ppm_image *image, ppm_image **contour_map, unsigned char **grid, int step_x) {
    for (int i = 0; i < CONTOUR_CONFIG_COUNT; i++) {
        free(contour_map[i]->data);
        free(contour_map[i]);
    }
    free(contour_map);

    for (int i = 0; i <= image->x / step_x; i++) {
        free(grid[i]);
    }
    free(grid);

    free(image->data);
    free(image);
}

void init_contour(arg_images *args){
    int start_init = (*args).thread_id * (double)CONTOUR_CONFIG_COUNT / (*args).nr_threads;
    int end_init;
    
    // Devides the steps of the loop between the threads used.
    if(((*args).thread_id + 1) * (double)CONTOUR_CONFIG_COUNT / (*args).nr_threads < CONTOUR_CONFIG_COUNT){
        end_init = ((*args).thread_id + 1) * (double)CONTOUR_CONFIG_COUNT / (*args).nr_threads;
    }else end_init = CONTOUR_CONFIG_COUNT;

    for (int i = start_init; i < end_init; i++) {
        char filename[FILENAME_MAX_SIZE];
        sprintf(filename, "./contours/%d.ppm", i);
        (*args).contour_map[i] = read_ppm(filename);
    }
}

void rescale_image(arg_images *args){
    uint8_t sample[3];


    int start_re = (*args).thread_id * (double)(*args).scaled_image->x / (*args).nr_threads;
    int end_re;
    
    // Devides the steps of the loop between the threads used.
    if(((*args).thread_id + 1) * (double)(*args).scaled_image->x / (*args).nr_threads < (*args).scaled_image->x){
        end_re = ((*args).thread_id + 1) * (double)(*args).scaled_image->x / (*args).nr_threads;
    }else end_re = (*args).scaled_image->x;

    for (int i = start_re; i < end_re; i++) {
        for (int j = 0; j < (*args).scaled_image->y; j++) {
            float u = (float)i / (float)((*args).scaled_image->x - 1);
            float v = (float)j / (float)((*args).scaled_image->y - 1);
            sample_bicubic((*args).image, u, v, sample);

            (*args).scaled_image->data[i * (*args).scaled_image->y + j].red = sample[0];
            (*args).scaled_image->data[i * (*args).scaled_image->y + j].green = sample[1];
            (*args).scaled_image->data[i * (*args).scaled_image->y + j].blue = sample[2];
        }
    }
}

void *thread_function(void *arg)
{
	arg_images args= *(arg_images *)arg;

	init_contour(&args);
    pthread_barrier_wait(args.barrier);
    
    if(args.rescale_flag == 0){
        rescale_image(&args);
        // After the image is rescaled the barrier waits for all threads to finish before going to the next step.
        pthread_barrier_wait(args.barrier);
    }

    build_grid(&args);
    pthread_barrier_wait(args.barrier);

    march(&args);

	pthread_exit(NULL);

    return NULL;
}


int main(int argc, char *argv[]) {
    if (argc < 4) {
        fprintf(stderr, "Usage: ./tema1 <in_file> <out_file> <P>\n");
        return 1;
    }

    // Number of threads used given as argument
    int nr_threads = atoi(argv[3]);
    arg_images args[nr_threads];

    // Pointers used as shared memory between the diffrent threads
    ppm_image *image = read_ppm(argv[1]);
    ppm_image *new_image;
    ppm_image **contour_map;
    unsigned char **grid;

    contour_map = (ppm_image **)malloc(CONTOUR_CONFIG_COUNT * sizeof(ppm_image *));
    if (!contour_map) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    // If the image needs to be rescaled the flag is left on 0 if not it's set to 1
    int rescale_flag = 0;
    if (image->x <= RESCALE_X && image->y <= RESCALE_Y){
        rescale_flag = 1;
        // If it does not need to be rescaled the new_image will point to the initial image.
        new_image = image;
    }else {

        // alloc memory for image
        new_image = (ppm_image *)malloc(sizeof(ppm_image));
        if (!new_image) {
            fprintf(stderr, "Unable to allocate memory\n");
            exit(1);
        }
        new_image->x = RESCALE_X;
        new_image->y = RESCALE_Y;

        new_image->data = (ppm_pixel*)malloc(new_image->x * new_image->y * sizeof(ppm_pixel));
        if (!new_image) {
            fprintf(stderr, "Unable to allocate memory\n");
            exit(1);
        }
    }
  
    int p = new_image->x / STEP;
    int q = new_image->y / STEP;

    grid = (unsigned char **)malloc((p + 1) * sizeof(unsigned char*));
    for (int i = 0; i <= p; i++) {
        grid[i] = (unsigned char *)malloc((q + 1) * sizeof(unsigned char));
        if (!grid[i]) {
            fprintf(stderr, "Unable to allocate memory\n");
            exit(1);
        }
    }

    if (!grid) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }


    int r;
    pthread_t threads[nr_threads];
    long id;
    void *status;

    pthread_barrier_t *barrier = malloc(sizeof(pthread_barrier_t));
    pthread_barrier_init(barrier, NULL, nr_threads);

    for (id = 0; id < nr_threads; id++) {
        // The thread_id will indicate the number of the current thread,
        // all of the images points towards the same shared memory,
        // the rescale_flag, number of threads and barrier are the same for all
        // of them.
        args[id].thread_id = id;
        args[id].image = image;
        args[id].scaled_image = new_image;
        args[id].contour_map = contour_map;
        args[id].grid = grid;
        args[id].nr_threads = nr_threads;
        args[id].rescale_flag = rescale_flag;
        args[id].barrier = barrier;

        r = pthread_create(&threads[id], NULL, thread_function, &args[id]);

        if (r) {
        printf("Eroare la crearea thread-ului %ld\n", id);
        exit(-1);
        }
    }

    // Threads are joined after they all finish running the function.
    for (id = 0; id < nr_threads; id++) {
        r = pthread_join(threads[id], &status);

        if (r) {
        printf("Eroare la asteptarea thread-ului %ld\n", id);
        exit(-1);
        }
    }

    write_ppm(new_image, argv[2]);

    free_resources(new_image, contour_map, grid, STEP);

    // The barrier is destroyed at the end of the program.
    pthread_barrier_destroy(barrier);

    return 0;
}
