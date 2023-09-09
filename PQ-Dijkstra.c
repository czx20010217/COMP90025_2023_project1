#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>      // for DBL_MAX
#include <time.h>       // for clock()
#include <sys/time.h>       // for clock()

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
    const unsigned long a = 16807;
    const unsigned long m = 2147483647;

    /* For debugging only */
    // return (seed);

    /* Real code */
    seed = -seed;       // Make high bits non-zero
    int res   = par->par1;
    int scale = par->par2;

    int cost;
    
    for (cost = 0; seed >> res != 0; cost++) {
        seed = (a * seed) % m;
    }

    return (10 + (cost >> (8 * sizeof(unsigned long) - res - scale))) / 10.0;
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

/******************************************************************************/

/* UNSEEN must be 0, as nodes initialized using memset */
#define UNSEEN 0
#define OPEN   1
#define CLOSED 2

void
a_star (double **board, int x_size, int y_size, params par)
{
    int x_end = x_size - 1;
    int y_end = y_size - 1;

    node *pivot;
    node **cand = init_cand (x_size, y_size);
    cand[0][0].cost = cell_cost (board[0][0], &par);

    while (!is_equal(pivot = pq_pop_min(), x_end, y_end)) {
        pivot->is_closed = CLOSED;

        /* Expand all neighbours */
        for (int dx = -1; dx <= 1; dx++) {
            for (int dy = -1; dy <= 1; dy++) {
                int new_x = pivot->x + dx;
                int new_y = pivot->y + dy;
                if (new_x < 0 || new_x > x_end || new_y < 0 || new_y > y_end
                              || (dx == 0 && dy == 0))
                    continue;
                if (!cand[new_x][new_y].is_closed) {
                    /* Note: this calculates costs multiple times */
                    /* You will probably want to avoid that, */
                    /* but this version is easy to parallelize. */
                    double node_cost = cell_cost(board[new_x][new_y], &par);
                    if (pivot->cost + node_cost < cand[new_x][new_y].cost) {
                        cand[new_x][new_y].cost = pivot->cost + node_cost;
                        cand[new_x][new_y].x = new_x;
                        cand[new_x][new_y].y = new_y;
                        cand[new_x][new_y].prev_x = pivot->x;
                        cand[new_x][new_y].prev_y = pivot->y;
                        /* Here we simply insert a better path into the PQ. */
                        /* It is more efficient to change the weight of */
                        /* the old entry, but this also works. */
                        pq_insert (&(cand[new_x][new_y]));
                    }
                }
            }
        }
        DEBUG(for (int i = 0; i < y_size; i++) { for (int j = 0; j < x_size; j++) { printf ("(%lg)%lg%c ", board[i][j], cand[i][j].cost, " _*"[cand[i][j].is_closed]); } printf ("\n"); })
    }

    node *p = &cand[x_end][y_end];
    while (!is_equal(p, 0, 0)) {
        printf ("%d %d %g %g\n", p->x, p->y, board[p->x][p->y], p->cost);
        p = &(cand[p->prev_x][p->prev_y]);
    }
    printf ("%d %d %g %g\n", 0, 0, board[0][0], p->cost);
}

/******************************************************************************/

int
main ()
{
    int x_size, y_size;
    double **board;
    node **open;
    int i,j;
    params par;

    assert_msg (scanf ("%d %d %d %d", &x_size, &y_size, &(par.par1), &(par.par2)) == 4, "Failed to read size");
    board = read_board(x_size, y_size);


    clock_t t = clock();
    double w = get_wall_time ();
    a_star (board, x_size, y_size, par);
    printf ("Time: %ld %lf\n", clock() - t, get_wall_time() - w);

    return 0;
}
