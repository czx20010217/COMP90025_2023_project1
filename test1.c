#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sched.h>
#include <omp.h>
#include <float.h>      // for DBL_MAX
#include <time.h>       // for clock()
#include<unistd.h>

#define MAX_NODES 100000

#define _inline
//#define _inline inline

//#define DEBUG(x) x
#define DEBUG(x)

//#define DEBUG1(x) x
#define DEBUG1(x)

/* UNSEEN must be 0, as nodes initialized using memset */
#define UNSEEN 0
#define OPEN   1
#define CLOSED 2

void
assert_msg (int cond, char *msg)
{
    if (!cond) {
        fprintf (stderr, "%s\n", msg);
        exit(1);
    }
}

typedef struct {
    int    x, y;
    int    prev_x, prev_y;
    int    is_closed;
    double cost;
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
    return (n1->cost > n2->cost);
}

/******************************************************************************/
/* Priority queue */
/* Entries are of type *node. */
/* Ordering is specified by function  greater(node *n1, node *n2) */
/* that returns 1 if *n1 < *n2 and 0 otherwise. */

/* left and right child offsets in priority queue */
#define PQ_LEFT  1
#define PQ_RIGHT 2

node *pq_array[MAX_NODES];
int pq_last;

void
pq_init ()
{
    pq_last = -1;       // More efficient to record last entry, not size
}

_inline int
pq_parent (int i)
{
    return (i-1) / 2;
}

_inline int
pq_child (int i, int left_right)
{
    return 2*i + left_right;
}

_inline void
pq_swap (node **n1, node **n2)
{
    node *tmp = *n1;
    *n1 = *n2;
    *n2 = tmp;
}

void
pq_up_heap (int i)
{
    while (i > 0 && greater (pq_array[pq_parent(i)], pq_array[i])) {
        pq_swap (&(pq_array[pq_parent(i)]), &(pq_array[i]));
        DEBUG1(printf ("U %d:", i); for (int j = 0; j <= pq_last; j++) { printf ("%d(%d,%d)%lg ", j, pq_array[j]->x, pq_array[j]->y, pq_array[j]->cost); } printf ("\n");)
        i = pq_parent (i);
    }
    DEBUG1(printf ("parent cost %lg my cost %lg\n", pq_array[pq_parent(i)]->cost, pq_array[i]->cost);)
    DEBUG1(printf ("After U %d:", i); for (int j = 0; j <= pq_last; j++) { printf ("(%d,%d)%lg ", pq_array[j]->x, pq_array[j]->y, pq_array[j]->cost); } printf ("\n");)
}


/* Find correct location for newly inserted node */
void
pq_down_heap (int i)
{
    int max_idx = i;

    int left = pq_child (i, PQ_LEFT);
    if (left <= pq_last && greater (pq_array[max_idx], pq_array[left])) {
        max_idx = left;
    }

    int right = pq_child (i, PQ_RIGHT);
    if (right <= pq_last && greater (pq_array[max_idx], pq_array[right])) {
        max_idx = right;
    }

    DEBUG1(printf ("D %d %d:", i, max_idx); for (int j = 0; j <= pq_last; j++) { printf ("%d(%d,%d)%lg ", j, pq_array[j]->x, pq_array[j]->y, pq_array[j]->cost); } printf ("\n");)
    if (i != max_idx) {
        pq_swap (&(pq_array[i]), &(pq_array[max_idx]));
        pq_down_heap (max_idx);
    }
}

void
pq_insert (node *n)
{
    DEBUG(printf ("<---(%d, %d):(%d,%d) %lg\n", n->x, n->y, n->prev_x, n->prev_y, n->cost);)

    pq_array[++pq_last] = n;
    pq_up_heap(pq_last);

    DEBUG1(printf ("I "); for (int i = 0; i <= pq_last; i++) { printf ("(%d,%d)%lg ", pq_array[i]->x, pq_array[i]->y, pq_array[i]->cost); } printf ("\n");)
}

node *
pq_pop_min ()
{
    node *retval = pq_array[0];

    DEBUG(printf (" -> (%d, %d):(%d,%d) %lg\n", retval->x, retval->y, retval->prev_x, retval->prev_y, retval->cost);)
    assert_msg (pq_last>=0, "Cannot pop from an empty priority queue");

    pq_array[0] = pq_array[pq_last--];
    pq_down_heap(0);

    DEBUG1(printf (" -> (%d, %d) %lf\n", retval->x, retval->y, retval->cost);)

    return retval;
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

int functionCallCount = 0;

double
cell_cost (long int seed, params *par)
{
    functionCallCount++;
    //const unsigned long a = 16807;
    //const unsigned long m = 2147483647;
    /* For debugging only */
    // return (seed);

    /* Real code */
    long a = 16807;
    long m = 2147483647;

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
        for (int j = 0; j < y_size; j++)
            cand[i][j].cost = DBL_MAX;
    }

    pq_init();
    pq_insert(&(cand[0][0]));

    return cand;
}

int
main ()
{
    printf ("statrted: \n");
    int x_size, y_size;
    double **board;
    double **cost_board;
    // node **open;
    // int i,j;
    params par;
    int some = 10;

    assert_msg (scanf ("%d %d %d %d", &x_size, &y_size, &(par.par1), &(par.par2)) == 4, "Failed to read size");
    board = read_board(x_size, y_size);
    cost_board = init_cost_board(x_size, y_size);

    node **cand = init_cand (x_size, y_size);
    cand[0][0].cost = board[0][0];
    int x_end = x_size - 1;
    int y_end = y_size - 1;
    node *pivot;

    for (int x = 0; x < x_size; x++){
        for (int y = 0; y < y_size; y++){
            double result = cell_cost(board[x][y], &par);
            printf("seed: %3.1f, result: %3.5f\n", board[x][y], result);

        }
    }
    

    printf("count ended: %d\n", functionCallCount);
    return 0;
}
