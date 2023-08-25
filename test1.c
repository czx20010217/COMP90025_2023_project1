#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sched.h>
#include <omp.h>
#include <float.h>      // for DBL_MAX
#include <time.h>       // for clock()

#define MAX_NODES 100000

#define _inline
//#define _inline inline

//#define DEBUG(x) x
#define DEBUG(x)

//#define DEBUG1(x) x
#define DEBUG1(x)

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

int main() {
    pq_init();
    #pragma omp parallel num_threads(16)
    {
        int thread_id = omp_get_thread_num();
        int cpu_num = sched_getcpu();
        int socket_id = thread_id % 2; // Simulated socket assignment

        #pragma omp critical
        {
            printf("Thread %d is on socket %d\n", thread_id, cpu_num);
        }
    }
    
    #pragma omp parallel num_threads(16)
    {
        int thread_id = omp_get_thread_num();
        node *cand_data = (node*)malloc (sizeof(node));

        #pragma omp critical
        {
            if (cand_data == NULL) {
                printf("Memory allocation failed on thread %d %d.\n", thread_id);
            }
            
            cand_data->cost = (double)thread_id;
            pq_insert(cand_data);
            printf("Thread %d is inserting number %f\n", thread_id, cand_data->cost);
        }
    }
    
    #pragma omp parallel num_threads(16)
    {
        int thread_id = omp_get_thread_num();

        #pragma omp critical
        {
            
            node* pivot = pq_pop_min();
            printf("Thread %d is poping number %f\n", thread_id, pivot->cost);
        }
    }

    return 0;
}
