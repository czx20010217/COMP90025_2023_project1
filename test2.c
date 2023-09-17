#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sched.h>
#include <omp.h>
#include <limits.h>
#include <float.h>      // for DBL_MAX
#include <time.h>       // for clock()
#include <sys/time.h>       // for clock()
#include <unistd.h>
#include <assert.h>

#define MAX_NODES 100000

#define _inline
//#define _inline inline

//#define DEBUG(x) x
#define DEBUG(x)

//#define DEBUG1(x) x
#define DEBUG1(x)

// From https://stackoverflow.com/questions/17432502/how-can-i-measure-cpu-time-and-wall-clock-time-on-both-linux-windows
double get_wall_time() {
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

int functionCallCount = 0;

void
assert_msg (int cond, char *msg)
{
    if (!cond) {
        fprintf (stderr, "%s\n", msg);
        exit(1);
    }
}

typedef struct{
    int **buckets_map;
    double **board;
    double **cost_board;
    int **deleted_map;
    int  in_bucket_current;
    int  in_buckets_total;
    int  current_bucket;
    int  processed_nodes;
    int  delta;
} status;

typedef struct {
    int    x, y;
    int    prev_x, prev_y;
    double cost;
    int bucket;
    int pos;
} node;

_inline int
is_equal (node* n, int x, int y)
{
    // Return true if n == NULL to ensure while loop terminates
    // even if pq_pop_min is modified to return NULL for an empty queue.
    return !n || (n->x == x && n->y == y);
}

_inline int
greater (node *n1, node *n2)
{
    return (n1->bucket > n2->bucket);
}

/******************************************************************************/
/*
 * Calculate the cost of each square in the grid, given its seed.
 * This is deliberately expensive so that overall program run-time is not
 * dominated by overheads.
 * More computationally expensive if res is smaller.
 * Wider range of costs if scale is larger.
 *
 * Based on Park and Miller's Oct 1988 CACM random number generator
 */

typedef struct {
    int par1, par2;
} params;

double
cell_cost (long int seed, params *par)
{
    #pragma omp critical
    functionCallCount++;
    const unsigned long a = 16807;
    const unsigned long m = 2147483647;
    /* For debugging only */
    // return (seed);

    /* Real code */

    seed = -seed;       // Make high bits non-zero
    int res   = par->par1;
    int scale = par->par2;

    int cost;
    // printf("current input: %ld, %d, %d\n", seed, res, scale);
    
    for (cost = 0; seed >> res != 0; cost++) {
        seed = (a * seed) % m;
    }
    // printf("finish\n");

    return (10 + (cost >> (8 * sizeof(unsigned long) - res - scale))) / 10.0;
    // return (seed);

}

double **
read_board (int x_size, int y_size)
{
    double **board = (double **)malloc (x_size * sizeof (*board));
    double *board_data = (double*)malloc (x_size*y_size * sizeof(*board_data));
    assert_msg(board != NULL && board_data != NULL, "Could not allocate board");

    for (int i = 0; i < x_size; i++) {
        board[i] = board_data + i*y_size;

        for (int j = 0; j < y_size; j++) {
            assert_msg (scanf ("%lf", &(board[i][j])) == 1, "Failed to read board");
        }
    }

    return board;
}

double **
init_cost_board (int x_size, int y_size)
{
    double **board = (double **)malloc (x_size * sizeof (*board));
    double *board_data = (double*)malloc (x_size*y_size * sizeof(*board_data));
    assert_msg(board != NULL && board_data != NULL, "Could not allocate board");

    for (int i = 0; i < x_size; i++) {
        board[i] = board_data + i*y_size;

        for (int j = 0; j < y_size; j++) {
            board[i][j] = -1;
        }
    }

    return board;
}

node **
init_cand (int x_size, int y_size)
{
    node **cand = (node **)malloc (x_size * sizeof (node*));
    node *cand_data = (node*)malloc (x_size * y_size * sizeof(node));
    assert_msg(cand != NULL && cand_data != NULL, "Could not allocate open");

    memset (cand_data, 0, y_size * x_size * sizeof(node));

    for (int i = 0; i < x_size; i++) {
        cand[i] = cand_data + i*y_size;
        for (int j = 0; j < y_size; j++){
            cand[i][j].bucket = -1;
            cand[i][j].pos = -1;
            cand[i][j].cost = DBL_MAX;
        }
            
    }

    return cand;
}

int relax(int x, int y, int new_x, int new_y, double new_cost, status* stats, node **cand){
    int bucket;
    int flag = 0;
    if(new_cost < cand[new_x][new_y].cost){
        //printf("in\n");
        
        if(stats->buckets_map[new_x][new_y] == -1){
            stats->in_buckets_total += 1;
            flag = 1;
        }
        

        cand[new_x][new_y].cost = new_cost;
        cand[new_x][new_y].x = new_x;
        cand[new_x][new_y].y = new_y;
        cand[new_x][new_y].prev_x = x;
        cand[new_x][new_y].prev_y = y;

        bucket = new_cost / stats->delta;
        // if (stats->buckets_map[new_x][new_y] == -1 || stats->buckets_map[new_x][new_y] > bucket){stats->buckets_map[new_x][new_y] = bucket;}
        stats->buckets_map[new_x][new_y] = bucket;
        
        
        if (bucket ==  stats->current_bucket) stats->in_bucket_current = 1;
        //printf("out\n");
    }
    return flag;
}



int get_min_bucket(int** node_list, int x_size, int y_size, status* stats){
    int min_bucket = 10000000;
    for (int x = 0; x < x_size; x++){
        for (int y = 0; y < y_size; y++){
            if (stats->buckets_map[x][y] != -1 && stats->buckets_map[x][y] < min_bucket){
                min_bucket = stats->buckets_map[x][y];
            }
            
        }
    }
    return min_bucket;
}

double get_cost(double** cost_board, double** seed_board, int x, int y, params *par){
    double cost;
    if (cost_board[x][y] == -1){
        cost = cell_cost(seed_board[x][y], par);
        cost_board[x][y] = cost;
    }else{
        cost = cost_board[x][y];
    }
}

int
main ()
{
    printf("count started: %d\n", functionCallCount);
    int x_size, y_size;
    double **board;
    double **cost_board;
    node **open;
    int i,j;
    params par;

    assert_msg (scanf ("%d %d %d %d", &x_size, &y_size, &(par.par1), &(par.par2)) == 4, "Failed to read size");
    board = read_board(x_size, y_size);
    cost_board = init_cost_board(x_size, y_size);

    int x_end = x_size - 1;
    int y_end = y_size - 1;
    status *stats = (status*) malloc(sizeof(status));

    clock_t t = clock();
    double w = get_wall_time ();

    // init bucket map and deleted_map
    stats->buckets_map = (int **)malloc(x_size * sizeof(int *));
    stats->deleted_map = (int **)malloc(x_size * sizeof(int *));
    for (int i = 0; i < x_size; i++) {
        stats->buckets_map[i] = (int *)malloc(y_size * sizeof(int));
        stats->deleted_map[i] = (int *)malloc(y_size * sizeof(int));
        if (stats->buckets_map[i] == NULL) {
            printf("Memory allocation failed.\n");
            return 1;
        }
    }

    stats->board = board;
    omp_lock_t lock[x_size][y_size];

    // #pragma omp parallel for
    for (int i = 0; i < x_size; i++){
        for (int j = 0; j < y_size; j++){
            stats->buckets_map[i][j] = -1;
            stats->deleted_map[i][j] = 0;
            omp_init_lock(&(lock[i][j]));
        }
            
    }

    stats->in_bucket_current = 0;
    stats->in_buckets_total = 1;
    stats->current_bucket = 0;
    stats->processed_nodes = 0;
    stats->delta = 100;
    int activeVertices = 1;
    int x, y, dx, dy;
    double cost, new_cost;
    int bucket;

    node **cand = init_cand (x_size, y_size);
    cost_board[0][0] = cell_cost(board[0][0], &par);
    cand[0][0].cost = cost_board[0][0];
    
    cost_board[x_end][0] = cell_cost(board[x_end][0], &par);
    cost_board[0][y_end] = cell_cost(board[0][y_end], &par);
    cost_board[x_end][y_end] = cell_cost(board[x_end][y_end], &par);
    
    stats->delta = (cost_board[x_end][y_end]+cost_board[0][0]+cost_board[x_end][0]+cost_board[0][y_end]) / (2*8);
    
    stats->buckets_map[0][0] = 0;
    printf("here\n");

    while (stats->in_buckets_total > 0){
        // for (int i = 0; i < y_size; i++) { for (int j = 0; j < x_size; j++) { printf ("(%lg)%5.1lg", board[i][j], cand[i][j].cost); } printf ("\n"); }
        // printf("current count: %d %d\n", functionCallCount, stats->in_buckets_total);
        // printf("current end cost: %5.1f\n", cand[x_end][y_end].cost);
        stats->processed_nodes += activeVertices;
        activeVertices = 0;
        stats->in_bucket_current = 1;
        while(stats->in_bucket_current){
            // for (int i = 0; i < y_size; i++) { for (int j = 0; j < x_size; j++) { printf ("(%lg)%5.1lg", board[i][j], cand[i][j].cost); } printf ("\n"); }
            stats->in_bucket_current = 0;
            int old_total = stats->in_buckets_total;
            int new_bucket = 0;
            // for (int i = 0; i < y_size; i++) { for (int j = 0; j < x_size; j++) { printf ("(%lg)%5.1lg", board[i][j], stats->buckets_map[i][j]); } printf ("\n"); }

            #pragma omp parallel for private(x, y, cost) shared(stats) collapse(2) schedule(dynamic) reduction(+ : new_bucket)
            for (x = 0; x < x_size; x++){
                for (y = 0; y < y_size; y++){
                    if(__sync_bool_compare_and_swap(&(stats->buckets_map[x][y]), stats->current_bucket, -1)){
                        stats->deleted_map[x][y] = 1;
                        // if (stats->in_buckets_total){printf("pos: %d %d\n", x, y);}



                        #pragma omp critical
                        {
                            stats->in_buckets_total--;
                            new_bucket += -1;
                        }
                        // printf("first time here\n");
                        // 8 direction visit
                        for (int dx = -1; dx <= 1; dx++) {
                            for (int dy = -1; dy <= 1; dy++) {
                                int new_x = x + dx;
                                int new_y = y + dy;
                                int result;
                                if (new_x < 0 || new_x > x_end || new_y < 0 || new_y > y_end
                                            || (dx == 0 && dy == 0))
                                    continue;


                                omp_set_lock(&(lock[new_x][new_y]));
                                if (cost_board[new_x][new_y] == -1){
                                    cost = cell_cost(board[new_x][new_y], &par);
                                    cost_board[new_x][new_y] = cost;
                                }else{
                                    cost = cost_board[new_x][new_y];
                                }
                                omp_unset_lock(&(lock[new_x][new_y]));

                                

                                if (cost <= stats->delta){
                                    // #pragma omp critical
                                    // {
                                    //     assert(cost <= stats->delta);
                                    //     result = relax(x, y, new_x, new_y, cand[x][y].cost + cost, stats, cand);
                                    // }
                                    omp_set_lock(&(lock[new_x][new_y]));
                                    #pragma omp critical
                                    {
                                        assert(cost <= stats->delta);
                                        result = relax(x, y, new_x, new_y, cand[x][y].cost + cost, stats, cand);
                                    }
                                    omp_unset_lock(&(lock[new_x][new_y]));


                                }

                                

                            }
                        }
                    }
                }
            }
            // for (x = 0; x < x_size; x++){
            //     for (y = 0; y < y_size; y++){
            //         assert(stats->buckets_map[x][y] == -1 || stats->buckets_map[x][y] >= stats->current_bucket);
            //     }
            // }
            // assert(stats->in_buckets_total - old_total == new_bucket);
        }

        int old_total = stats->in_buckets_total;
        int new_bucket = 0;
        #pragma omp parallel for private(x, y, cost) shared(stats) collapse(2)
        for (x = 0; x < x_size; x++){
            for (y = 0; y < y_size; y++){
                if (__sync_bool_compare_and_swap(&(stats->deleted_map[x][y]), 1, -1)){
                    for (int dx = -1; dx <= 1; dx++) {
                        for (int dy = -1; dy <= 1; dy++) {
                            int new_x = x + dx;
                            int new_y = y + dy;
                            int result;
                            if (new_x < 0 || new_x > x_end || new_y < 0 || new_y > y_end
                                        || (dx == 0 && dy == 0))
                                continue;

                            if (cost_board[new_x][new_y] == -1){
                                cost = cell_cost(board[new_x][new_y], &par);
                                cost_board[new_x][new_y] = cost;
                            }else{
                                cost = cost_board[new_x][new_y];
                            }
                            // check heavy
                            if (cost > stats->delta){
                                #pragma omp critical
                                {
                                    assert(cost > stats->delta);
                                    result = relax(x, y, new_x, new_y, cand[x][y].cost + cost, stats, cand);
                                }

                                if (result){
                                    new_bucket += 1;
                                }
                            }
                            
                        }
                    }
                }
                
            }
        }
        // //assert(stats->in_buckets_total - old_total == new_bucket);
        // assert(stats->in_bucket_current == 0);
        // for (x = 0; x < x_size; x++){
        //     for (y = 0; y < y_size; y++){
        //         assert(stats->buckets_map[x][y] == -1 || stats->buckets_map[x][y] > stats->current_bucket);
        //     }
        // }
        stats->current_bucket++;
        
    }
    node *p = &cand[x_end][y_end];
    printf("finalcost: %f\n", p->cost);


    while (!is_equal(p, 0, 0)) {
        printf ("%d %d %g %g\n", p->x, p->y, board[p->x][p->y], p->cost);
        // printf("current cost: %f\n", p->cost);
        p = &(cand[p->prev_x][p->prev_y]);
    }
    printf ("%d %d %g %g\n", 0, 0, board[0][0], p->cost);
    printf ("Time: %ld %lf\n", clock() - t, get_wall_time() - w);

    printf("count ended: %d\n", functionCallCount);
    //for (int i = 0; i < y_size; i++) { for (int j = 0; j < x_size; j++) { printf ("(%lg)%5.1lg", board[i][j], stats->buckets_map[i][j]); } printf ("\n"); }

    return 0;
}