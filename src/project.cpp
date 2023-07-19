#include <cstdint>
#include <cstdio>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include <mpi.h>
#include <pthread.h>
#include <unistd.h>
#include <allegro.h>
#include "WorkPool.h"


#define v(r, c) ((r) * (N_COLS_PROC + 2) + (c))  //macro to access the matrices on which the workers operate
#define w(r, c) ((r) * (NCOLS) + (c))  //macro to access the matrix on which the display process operates


//Constants
#define DISPLAY_PROC_RANK 0
#define FIRST_WORKER_RANK ((!BENCHMARK_MODE) ? 1 : 0)
#define LAST_WORKER_RANK ((!BENCHMARK_MODE) ? N_WORKER_PROC : (N_WORKER_PROC - 1))

#define COMM_TAG 999

//Graphics constants
#define INIT_END_GRAPHICS_DELAY 3000000  //us
#define UPDATE_GRAPHICS_DELAY 0  //us
#define CELL_DIM 8

//Configuration filename
const char CONFIG_FILENAME[] = "Configuration.txt";



//Global variables
int RANK;
uint32_t RANK_LEFT, RANK_RIGHT;
uint8_t *readM, *writeM;

//Configuration variables
bool BENCHMARK_MODE;
uint32_t NROWS, NCOLS;
uint32_t N_COLS_PROC, N_CELLS_PROC;
uint32_t X_PART, Y_PART;
uint32_t N_THREADS, N_CELLS_PART;
uint64_t N_STEP;
std::string INPUT_FIRST_FILENAME;
std::string INPUT_LAST_FILENAME;

//MPI
int N_PROC;
uint32_t N_WORKER_PROC;
MPI_Datatype halo_col_type;
MPI_Datatype matrix_proc_partition_type;
MPI_Request send_request_left, send_request_right;
MPI_Request recv_request_left, recv_request_right;

//pthreads
pthread_t main_thread_id;
pthread_t *threads;
pthread_barrier_t barrier;

//Pool of works
WorkPool inside_pool;
WorkPool border_pool;

//Graphics
BITMAP *bitmap;
int32_t black_color, white_color, lime_color;

//Benchmark times
double start_time, end_time;



//Functions prototypes
bool read_configuration_file();
void print_config_parameters();
bool check_config();
uint32_t get_fixed_xpart();
bool check_num_workers();
void calculate_neighbors();

void init_display_proc();
void init_worker();
void finalize_workers();
void read_CA(const std::string);

void fill_work_pools();
bool is_halo_cell(uint32_t cell);
bool is_border_cell(uint32_t cell);

void init_threads();
void* compute_transitions(void *);
void transitionFunction(int64_t row, int64_t col);
void swap_matrix();

void send_borders();
void recv_borders();
void wait_for_send_borders();
void wait_for_recv_borders();
void send_matrix_proc_partition();
void recv_matrix_proc_partition();

void init_graphics();
void finalize_graphics();
void drawWithAllegro(uint64_t curr_step);





int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &N_PROC);
    MPI_Comm_rank(MPI_COMM_WORLD, &RANK);

    if (!read_configuration_file()) {
        printf("Error opening configuration file!\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    if (!check_config()) {
        printf("The configuration used is invalid.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if (!check_num_workers()) {
        printf("The configuration file used requires the use of %u MPI worker processes!\n", Y_PART);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if (!BENCHMARK_MODE && RANK == DISPLAY_PROC_RANK) {  //if i am the display process
        print_config_parameters();
        init_display_proc();
        

        recv_matrix_proc_partition();
        drawWithAllegro(0);
        usleep(INIT_END_GRAPHICS_DELAY);

        for (uint64_t i = 1; i <= N_STEP; ++i) {
            recv_matrix_proc_partition();
            drawWithAllegro(i);
            usleep(UPDATE_GRAPHICS_DELAY);
        }
        
        usleep(INIT_END_GRAPHICS_DELAY);
        finalize_graphics();
    }
    else {  //if i am the worker process
        init_worker();

        init_threads();
        finalize_workers();
    }


    if (BENCHMARK_MODE && RANK == DISPLAY_PROC_RANK) {
        double elapsed_time_ms = (end_time - start_time) * 1000;
        printf("Time spent computing %lu steps: %f ms\n", N_STEP, elapsed_time_ms);
        printf("Average time per single step: %f ms\n", elapsed_time_ms / N_STEP);
    }

    MPI_Finalize();
    return 0;
}
END_OF_MAIN()





bool read_configuration_file() {
    std::ifstream config_file(CONFIG_FILENAME);
    if (!config_file.is_open())
        return false;
    
    config_file >> BENCHMARK_MODE;
    config_file >> NROWS;
    config_file >> NCOLS;
    config_file >> X_PART;
    config_file >> Y_PART;
    config_file >> N_THREADS;
    config_file >> N_STEP;
    config_file >> INPUT_FIRST_FILENAME;
    config_file >> INPUT_LAST_FILENAME;


    N_WORKER_PROC = (!BENCHMARK_MODE) ? N_PROC - 1: N_PROC;
    N_COLS_PROC = (NCOLS / Y_PART);
    N_CELLS_PROC = NROWS * N_COLS_PROC;

    if (N_CELLS_PROC % X_PART)
        X_PART = get_fixed_xpart();

    if (X_PART < N_THREADS)
        N_THREADS = X_PART;
    
    N_CELLS_PART = N_CELLS_PROC / X_PART;

    config_file.close();
    return true;
}

void print_config_parameters() {
    printf("Configuration parameters:\n");
    printf("2d Matrix Sizes: %ux%u (rows * cols)\n", NROWS, NCOLS);
    printf("Number of Process partitons (along Y Axis): %u\n", Y_PART);
    printf("Number of Thread partitions (along X Axis): %u\n", X_PART);
    printf("Number of Threads per MPI process: %u\n", N_THREADS);
    printf("Number of steps to compute: %lu\n", N_STEP);
    printf("\n");
}

bool check_config() {
    //if the number of cols isn't divisible by the number of Y Partition
    if (NCOLS % Y_PART)
        return false;

    return true;
}

uint32_t get_fixed_xpart() {
    for(size_t var = 1; ; ++var) {
        if (N_CELLS_PROC % (X_PART + var) == 0)
            return X_PART + var;
        if (var < X_PART && N_CELLS_PROC % (X_PART - var) == 0)
            return X_PART - var;
    }
}

bool check_num_workers() {
    return N_WORKER_PROC == Y_PART;
}

void calculate_neighbors() {
    RANK_LEFT = (RANK != FIRST_WORKER_RANK) ? RANK - 1 : -1;
    RANK_RIGHT = (RANK != LAST_WORKER_RANK) ? RANK + 1 : -1;
}





void init_display_proc() {
    init_graphics();

    readM = new uint8_t[NROWS * NCOLS];

    MPI_Type_vector(NROWS, N_COLS_PROC, NCOLS, MPI_UINT8_T, &matrix_proc_partition_type);
    MPI_Type_commit(&matrix_proc_partition_type);
}

void init_worker() {
    calculate_neighbors();

    main_thread_id = pthread_self();
    pthread_barrier_init(&barrier, NULL, N_THREADS);

    //Init the read and write Matrix
    readM = new uint8_t[NROWS * (N_COLS_PROC + 2)];
    writeM = new uint8_t[NROWS * (N_COLS_PROC + 2)];

    //Init the 2d array
    for (uint32_t i = 0; i < NROWS; ++i) {
        for (uint32_t j = 0; j < N_COLS_PROC + 2; ++j) {
            readM[v(i, j)] = 0;
            writeM[v(i, j)] = 0;
        }
    }
    
    if (RANK == FIRST_WORKER_RANK)
        read_CA(INPUT_FIRST_FILENAME);

    //Init MPI data structures
    MPI_Type_vector(NROWS, 1, N_COLS_PROC + 2, MPI_UINT8_T, &halo_col_type);
    MPI_Type_commit(&halo_col_type);
    MPI_Type_vector(NROWS, N_COLS_PROC, N_COLS_PROC + 2, MPI_UINT8_T, &matrix_proc_partition_type);
    MPI_Type_commit(&matrix_proc_partition_type);

    fill_work_pools();
}

void finalize_workers() {
    MPI_Type_free(&halo_col_type);
    pthread_barrier_destroy(&barrier);

    delete []readM;
    delete []writeM;

    MPI_Type_free(&matrix_proc_partition_type);
}

void read_CA(const std::string filename) {
    std::ifstream file(filename);
    std::string line;

    uint32_t row = 0;
    while (std::getline(file, line) && row < NROWS) {
        std::istringstream iss(line);
        std::string value;
        uint32_t col = 1;
        while (std::getline(iss, value, '\t') && col < N_COLS_PROC + 1) {
            readM[v(row, col)] = std::stoi(value);
            ++col;
        }
        ++row;
    }
}





void fill_work_pools() {
    uint32_t current_cell = 1;
    
    for(uint32_t n_part = 0; n_part < X_PART; ++n_part) {
        border_pool.new_partition();
        inside_pool.new_partition();
        
        for(uint32_t i = 0; i < N_CELLS_PART; ++i) {
            while (is_halo_cell(current_cell))
                ++current_cell;
            
            if (is_border_cell(current_cell))
                border_pool.push_work(current_cell);
            else
                inside_pool.push_work(current_cell);
            
            ++current_cell;
        }
    }

    border_pool.pop_back_empty_partition();
    inside_pool.pop_back_empty_partition();
}

inline bool is_halo_cell(uint32_t cell) {
    uint32_t cell_col = cell % (N_COLS_PROC + 2);
    return cell_col == 0 || cell_col == N_COLS_PROC + 1;
}

inline bool is_border_cell(uint32_t cell) {
    uint32_t cell_col = cell % (N_COLS_PROC + 2);
    return (RANK_LEFT != -1 && cell_col == 1) || (RANK_RIGHT != -1 && cell_col == N_COLS_PROC);
}





void init_threads() {
    const uint32_t N_THREADS_EXT = N_THREADS - 1;
    if (N_THREADS_EXT > 0) {
        threads = new pthread_t[N_THREADS_EXT];

        for (uint32_t i = 0; i < N_THREADS_EXT; ++i)
            if (pthread_create(threads + i, NULL, &compute_transitions, NULL))
                MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    send_matrix_proc_partition();
    compute_transitions(NULL);

    for (uint32_t i = 0; i < N_THREADS_EXT; ++i)
        pthread_join(threads[i], NULL);
    
    delete []threads;
}

void* compute_transitions(void* args) {
    uint64_t steps = N_STEP;
    const std::vector<uint64_t> *my_work;

    if (BENCHMARK_MODE) {
        pthread_barrier_wait(&barrier);
        if (pthread_equal(pthread_self(), main_thread_id)) {
            MPI_Barrier(MPI_COMM_WORLD);
            start_time = MPI_Wtime();
        }
        pthread_barrier_wait(&barrier);
    }
        

    while (steps-- > 0) {
        if (pthread_equal(pthread_self(), main_thread_id)) {
            send_borders();
            recv_borders();
        }
        
        while ((my_work = inside_pool.get_work()) != nullptr) {
            for (size_t i = 0; i < my_work->size(); ++i) {
                uint64_t work_cell = my_work->at(i);
                uint32_t row = work_cell / (N_COLS_PROC + 2);
                uint32_t col = work_cell % (N_COLS_PROC + 2);

                transitionFunction(row, col);
            }
        }

        if (pthread_equal(pthread_self(), main_thread_id))
            wait_for_recv_borders();
        
        pthread_barrier_wait(&barrier);

        
        while ((my_work = border_pool.get_work()) != nullptr) {
            for (size_t i = 0; i < my_work->size(); ++i) {
                uint64_t work_cell = my_work->at(i);
                uint32_t row = work_cell / (N_COLS_PROC + 2);
                uint32_t col = work_cell % (N_COLS_PROC + 2);

                transitionFunction(row, col);
            }
        }

        if (pthread_equal(pthread_self(), main_thread_id))
            wait_for_send_borders();

        pthread_barrier_wait(&barrier);

        if (pthread_equal(pthread_self(), main_thread_id)) {
            swap_matrix();
            send_matrix_proc_partition();
            inside_pool.reset_work_completion();
            border_pool.reset_work_completion();
        }

        pthread_barrier_wait(&barrier);
    }

    if (BENCHMARK_MODE) {
        pthread_barrier_wait(&barrier);
        if (pthread_equal(pthread_self(), main_thread_id)) {
            MPI_Barrier(MPI_COMM_WORLD);
            end_time = MPI_Wtime();
        }
        pthread_barrier_wait(&barrier);
    }

    return NULL;
}

inline void transitionFunction(int64_t r, int64_t c) {
    uint32_t cont = 0;
    for (int8_t dr = -1; dr < 2; ++dr)
        for (int8_t dc = -1; dc < 2; ++dc)
            if ((dr != 0 || dc != 0) && r+dr >= 0 && r+dr < NROWS && c+dc >= 0 && c+dc < (N_COLS_PROC + 2) && readM[v(r+dr, c+dc)] == 1)
                ++cont;

    if (readM[v(r, c)] == 1)
        writeM[v(r, c)] = (cont == 2 || cont == 3) ? 1 : 0;
    else
        writeM[v(r, c)] = (cont == 3) ? 1 : 0;
}

inline void swap_matrix() {
    uint8_t *tmp = readM;
    readM = writeM;
    writeM = tmp;
}





inline void send_borders() {
    if (RANK_LEFT != -1)
        MPI_Isend(&readM[v(0, 1)], 1, halo_col_type, RANK_LEFT, COMM_TAG, MPI_COMM_WORLD, &send_request_left);
    
    if (RANK_RIGHT != -1)
        MPI_Isend(&readM[v(0, N_COLS_PROC)], 1, halo_col_type, RANK_RIGHT, COMM_TAG, MPI_COMM_WORLD, &send_request_right);
}

inline void recv_borders() {
    if (RANK_LEFT != -1)
        MPI_Irecv(&readM[v(0, 0)], 1, halo_col_type, RANK_LEFT, COMM_TAG, MPI_COMM_WORLD, &recv_request_left);
    
    if (RANK_RIGHT != -1)
        MPI_Irecv(&readM[v(0, N_COLS_PROC + 1)], 1, halo_col_type, RANK_RIGHT, COMM_TAG, MPI_COMM_WORLD, &recv_request_right);
}

inline void wait_for_send_borders() {
    MPI_Status status;

    if (RANK_LEFT != -1)
        MPI_Wait(&send_request_left, &status);
    
    if (RANK_RIGHT != -1)
        MPI_Wait(&send_request_right, &status);
}

inline void wait_for_recv_borders() {
    MPI_Status status;

    if (RANK_LEFT != -1)
        MPI_Wait(&recv_request_left, &status);
    
    if (RANK_RIGHT != -1)
        MPI_Wait(&recv_request_right, &status);
}

inline void send_matrix_proc_partition() {
    //blocking send
    if (BENCHMARK_MODE)
        return;
    MPI_Send(&readM[v(0, 1)], 1, matrix_proc_partition_type, DISPLAY_PROC_RANK, COMM_TAG, MPI_COMM_WORLD);
}

inline void recv_matrix_proc_partition() {
    MPI_Status status;
    for (uint32_t i = 1; i < N_PROC; ++i)
        MPI_Recv(&readM[w(0, (i-1) * N_COLS_PROC)], 1, matrix_proc_partition_type, i, COMM_TAG, MPI_COMM_WORLD, &status);
    fflush(stdout);
}





void init_graphics() {
    allegro_init();
    install_keyboard();
    set_color_depth(8);
    bitmap = create_bitmap(NCOLS * CELL_DIM, NROWS * CELL_DIM);
    set_gfx_mode(GFX_AUTODETECT_WINDOWED, NCOLS * CELL_DIM, NROWS * CELL_DIM, 0, 0);

    black_color = makecol(0, 0, 0);
    white_color = makecol(0xFF, 0xFF, 0xFF);
    lime_color = makecol(0x65, 0xFE, 0x08);
}

void finalize_graphics() {
    destroy_bitmap(bitmap);
    allegro_exit();
}

inline void drawWithAllegro(uint64_t step) {
    for (uint32_t i = 0; i < NROWS; ++i) {
        for (uint32_t j = 0; j < NCOLS; ++j) {
            uint64_t x = j * CELL_DIM;
            uint64_t y = i * CELL_DIM;
            rectfill(bitmap, x, y, x+CELL_DIM, y+CELL_DIM, (readM[w(i, j)] == 1) ? white_color : black_color);
        }
    }

    for (uint32_t i = 1; i < Y_PART; ++i)
        line(bitmap, i * N_COLS_PROC * CELL_DIM, 0, i * N_COLS_PROC * CELL_DIM, (NROWS - 1) * CELL_DIM + CELL_DIM - 1, lime_color);


    textprintf_ex(bitmap, font, 10, 10, lime_color, black_color, "Step: %lu", step);
    blit(bitmap, screen, 0, 0, 0, 0, NCOLS * CELL_DIM, NROWS * CELL_DIM);
}