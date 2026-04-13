/*
 *  bpt.c
 */
#define Version "1.16.1"
/*
 *
 *  bpt:  B+ Tree Implementation
 *
 *  Copyright (c) 2018  Amittai Aviram  http://www.amittai.com
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  1. Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.

 *  3. The name of the copyright holder may not be used to endorse
 *  or promote products derived from this software without specific
 *  prior written permission.

 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.

 *  Author:  Amittai Aviram
 *    http://www.amittai.com
 *    amittai.aviram@gmail.com or afa13@columbia.edu
 *  Original Date:  26 June 2010
 *  Last modified: 02 September 2018
 *
 *  This implementation demonstrates the B+ tree data structure
 *  for educational purposes, including insertion, deletion, search, and display
 *  of the search path, the leaves, or the whole tree.
 *
 *  Must be compiled with a C99-compliant C compiler such as the latest GCC.
 *
 */

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <sys/time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define DEFAULT_ORDER       4
#define MIN_ORDER           3
#define MAX_ORDER           20
#define BUFFER_SIZE         256
#define ALIGNMET            (1UL << 21)

#define DEFAULT_NELEMENTS   (3400UL << 20)
#define DEFAULT_TREE_ORDER  16

#ifdef _OPENMP
#define DEFAULT_NLOOKUP     2000000000UL
#else
#define DEFAULT_NLOOKUP     50000000UL
#endif

size_t allocator_stat = 0;

static inline void *allocate(size_t size)
{
	void *memptr;
	if (posix_memalign(&memptr, ALIGNMET, size)) {
		printf("ENOMEM\n");
		exit(1);
	}
	allocator_stat += size;
	memset(memptr, 0, size);
	return memptr;
}

static inline void *allocate_align64(size_t size)
{
	void *memptr;
	if (posix_memalign(&memptr, 64, size)) {
		printf("ENOMEM\n");
		exit(1);
	}
	allocator_stat += size;
	memset(memptr, 0, size);
	return memptr;
}

typedef struct record {
	union {
		uint64_t value;
		struct record *next;
	};
	uint64_t flags;
} record;

typedef struct node {
	void        **pointers;
	uint64_t     *keys;
	struct node  *parent;
	bool          is_leaf;
	uint64_t      num_keys;
	struct node  *next;
} node;

struct element {
	uint64_t payload;
	uint64_t payload2;
};

uint64_t order        = DEFAULT_ORDER;
node    *queue        = NULL;
bool     verbose_output = false;

void     license_notice(void);
void     usage_1(void);
void     usage_2(void);
void     usage_3(void);
void     enqueue(node *new_node);
node    *dequeue(void);
uint64_t height(node * const root);
uint64_t path_to_root(node * const root, node *child);
void     print_leaves(node * const root);
void     print_tree(node * const root);
void     find_and_print(node * const root, uint64_t key, bool verbose);
void     find_and_print_range(node * const root, uint64_t range1, uint64_t range2, bool verbose);
uint64_t find_range(node * const root, uint64_t key_start, uint64_t key_end, bool verbose,
                    uint64_t returned_keys[], void *returned_pointers[]);
node    *find_leaf(node * const root, uint64_t key, bool verbose);
record  *find(node *root, uint64_t key, bool verbose, node **leaf_out);
uint64_t cut(uint64_t length);

record  *make_record(uint64_t value);
node    *make_node(void);
node    *make_leaf(void);
uint64_t get_left_index(node *parent, node *left);
node    *insert_into_leaf(node *leaf, uint64_t key, record *pointer);
node    *insert_into_leaf_after_splitting(node *root, node *leaf, uint64_t key, record *pointer);
node    *insert_into_node(node *root, node *parent, uint64_t left_index, uint64_t key, node *right);
node    *insert_into_node_after_splitting(node *root, node *parent, uint64_t left_index, uint64_t key, node *right);
node    *insert_into_parent(node *root, node *left, uint64_t key, node *right);
node    *insert_into_new_root(node *left, uint64_t key, node *right);
node    *start_new_tree(uint64_t key, record *pointer);
node    *insert(node *root, uint64_t key, uint64_t value);

uint64_t get_neighbor_index(node *n);
node    *adjust_root(node *root);
node    *coalesce_nodes(node *root, node *n, node *neighbor, uint64_t neighbor_index, uint64_t k_prime);
node    *redistribute_nodes(node *root, node *n, node *neighbor, uint64_t neighbor_index,
                            uint64_t k_prime_index, uint64_t k_prime);
node    *delete_entry(node *root, node *n, uint64_t key, void *pointer);
node    *delete(node *root, uint64_t key);

void license_notice(void) {
	printf("bpt version %s -- Copyright (c) 2018  Amittai Aviram "
	       "http://www.amittai.com\n", Version);
	printf("This program comes with ABSOLUTELY NO WARRANTY.\n"
	       "This is free software, and you are welcome to redistribute it\n"
	       "under certain conditions.\n"
	       "Please see the headnote in the source code for details.\n");
}

void usage_1(void) {
	printf("B+ Tree of Order %ld.\n", order);
	printf("Following Silberschatz, Korth, Sidarshan, Database Concepts, 5th ed.\n\n");
	printf("(%d <= order <= %d).\n", MIN_ORDER, MAX_ORDER);
}

void usage_2(void) {
	printf("Enter any of the following commands after the prompt > :\n"
	       "\ti <k>  -- Insert <k> as both key and value.\n"
	       "\ti <k> <v> -- Insert value <v> as the value of key <k>.\n"
	       "\tf <k>  -- Find the value under key <k>.\n"
	       "\tp <k>  -- Print the path from the root to key <k>.\n"
	       "\tr <k1> <k2> -- Print keys and values in range [<k1>, <k2>].\n"
	       "\td <k>  -- Delete key <k> and its associated value.\n"
	       "\tx -- Destroy the whole tree.\n"
	       "\tt -- Print the B+ tree.\n"
	       "\tl -- Print the keys of the leaves.\n"
	       "\tv -- Toggle verbose output.\n"
	       "\tq -- Quit.\n"
	       "\t? -- Print this help message.\n");
}

void usage_3(void) {
	printf("Usage: ./btree [--nelements <n>] [--nlookup <n>] [--order <n>]\n");
	printf("\twhere %d <= order <= %d\n", MIN_ORDER, MAX_ORDER);
	printf("\tnelements: number of elements to insert (default: %zu)\n", DEFAULT_NELEMENTS);
	printf("\tnlookup:   number of lookups to perform (default: %zu)\n", DEFAULT_NLOOKUP);
	printf("\torder:     B+ tree order                (default: %d)\n",  DEFAULT_TREE_ORDER);
}

void enqueue(node *new_node) {
	node *c;
	if (queue == NULL) {
		queue = new_node;
		queue->next = NULL;
	} else {
		c = queue;
		while (c->next != NULL)
			c = c->next;
		c->next = new_node;
		new_node->next = NULL;
	}
}

node *dequeue(void) {
	node *n = queue;
	queue = queue->next;
	n->next = NULL;
	return n;
}

void print_leaves(node * const root) {
	if (root == NULL) { printf("Empty tree.\n"); return; }
	uint64_t i;
	node *c = root;
	while (!c->is_leaf) c = c->pointers[0];
	while (true) {
		for (i = 0; i < c->num_keys; i++) {
			if (verbose_output) printf("%p ", c->pointers[i]);
			printf("%ld ", c->keys[i]);
		}
		if (verbose_output) printf("%p ", c->pointers[order - 1]);
		if (c->pointers[order - 1] != NULL) { printf(" | "); c = c->pointers[order - 1]; }
		else break;
	}
	printf("\n");
}

uint64_t height(node * const root) {
	uint64_t h = 0;
	node *c = root;
	while (!c->is_leaf) { c = c->pointers[0]; h++; }
	return h;
}

uint64_t path_to_root(node * const root, node *child) {
	uint64_t length = 0;
	node *c = child;
	while (c != root) { c = c->parent; length++; }
	return length;
}

void print_tree(node * const root) {
	node *n = NULL;
	uint64_t i = 0, rank = 0, new_rank = 0;
	if (root == NULL) { printf("Empty tree.\n"); return; }
	queue = NULL;
	enqueue(root);
	while (queue != NULL) {
		n = dequeue();
		if (n->parent != NULL && n == n->parent->pointers[0]) {
			new_rank = path_to_root(root, n);
			if (new_rank != rank) { rank = new_rank; printf("\n"); }
		}
		if (verbose_output) printf("(%p)", n);
		for (i = 0; i < n->num_keys; i++) {
			if (verbose_output) printf("%p ", n->pointers[i]);
			printf("%ld ", n->keys[i]);
		}
		if (!n->is_leaf)
			for (i = 0; i <= n->num_keys; i++)
				enqueue(n->pointers[i]);
		if (verbose_output) {
			if (n->is_leaf) printf("%p ", n->pointers[order - 1]);
			else            printf("%p ", n->pointers[n->num_keys]);
		}
		printf("| ");
	}
	printf("\n");
}

void find_and_print(node * const root, uint64_t key, bool verbose) {
	record *r = find(root, key, verbose, NULL);
	if (r == NULL) printf("Record not found under key %ld.\n", key);
	else           printf("Record at %p -- key %ld, value %ld.\n", r, key, r->value);
}

void find_and_print_range(node * const root, uint64_t key_start, uint64_t key_end, bool verbose) {
	uint64_t i;
	uint64_t array_size = key_end - key_start + 1;
	uint64_t returned_keys[array_size];
	void *returned_pointers[array_size];
	uint64_t num_found = find_range(root, key_start, key_end, verbose,
	                                returned_keys, returned_pointers);
	if (!num_found) printf("None found.\n");
	else for (i = 0; i < num_found; i++)
		printf("Key: %ld   Location: %p  Value: %ld\n",
		       returned_keys[i], returned_pointers[i],
		       ((record *)returned_pointers[i])->value);
}

uint64_t find_range(node * const root, uint64_t key_start, uint64_t key_end, bool verbose,
                    uint64_t returned_keys[], void *returned_pointers[]) {
	uint64_t i, num_found = 0;
	node *n = find_leaf(root, key_start, verbose);
	if (n == NULL) return 0;
	for (i = 0; i < n->num_keys && n->keys[i] < key_start; i++) ;
	if (i == n->num_keys) return 0;
	while (n != NULL) {
		for (; i < n->num_keys && n->keys[i] <= key_end; i++) {
			returned_keys[num_found]    = n->keys[i];
			returned_pointers[num_found] = n->pointers[i];
			num_found++;
		}
		n = n->pointers[order - 1];
		i = 0;
	}
	return num_found;
}

node *find_leaf(node * const root, uint64_t key, bool verbose) {
	if (root == NULL) {
		if (verbose) printf("Empty tree.\n");
		return root;
	}
	uint64_t i = 0;
	node *c = root;
	while (!c->is_leaf) {
		if (verbose) {
			printf("[");
			for (i = 0; i < c->num_keys - 1; i++) printf("%ld ", c->keys[i]);
			printf("%ld] ", c->keys[i]);
		}
		i = 0;
		while (i < c->num_keys) {
			if (key >= c->keys[i]) i++;
			else break;
		}
		if (verbose) printf("%ld ->\n", i);
		c = (node *)c->pointers[i];
	}
	if (verbose) {
		printf("Leaf [");
		for (i = 0; i < c->num_keys - 1; i++) printf("%ld ", c->keys[i]);
		printf("%ld] ->\n", c->keys[i]);
	}
	return c;
}

record *find(node *root, uint64_t key, bool verbose, node **leaf_out) {
	if (root == NULL) {
		if (leaf_out != NULL) *leaf_out = NULL;
		return NULL;
	}
	uint64_t i = 0;
	node *leaf = find_leaf(root, key, verbose);
	for (i = 0; i < leaf->num_keys; i++)
		if (leaf->keys[i] == key) break;
	if (leaf_out != NULL) *leaf_out = leaf;
	if (i == leaf->num_keys) return NULL;
	return (record *)leaf->pointers[i];
}

uint64_t cut(uint64_t length) {
	return (length % 2 == 0) ? length / 2 : length / 2 + 1;
}

#define NODE_SLAB_GROW   (1 << 18)
#define RECORD_SLAB_GROW (1 << 18)

struct node   *free_nodes = NULL;
struct record *free_recs  = NULL;

node *alloc_node(void) {
	if (!free_nodes) {
		node *n = allocate(NODE_SLAB_GROW);
		for (size_t i = 0; i < NODE_SLAB_GROW / sizeof(struct node); i++) {
			n[i].next = free_nodes;
			free_nodes = &n[i];
		}
	}
	node *nd = free_nodes;
	free_nodes = nd->next;
	return nd;
}

void free_node(node *n) {
	n->next = free_nodes;
	free_nodes = n;
}

record *alloc_record(void) {
	if (!free_recs) {
		record *r = allocate(RECORD_SLAB_GROW);
		for (size_t i = 0; i < RECORD_SLAB_GROW / sizeof(struct record); i++) {
			r[i].next = free_recs;
			free_recs = &r[i];
		}
	}
	record *rec = free_recs;
	free_recs = rec->next;
	return rec;
}

void free_record(record *r) {
	r->next = free_recs;
	free_recs = r;
}

record *make_record(uint64_t value) {
	record *new_record = alloc_record();
	if (new_record == NULL) { perror("Record creation."); exit(EXIT_FAILURE); }
	new_record->value = value;
	return new_record;
}

node *make_node(void) {
	node *new_node = alloc_node();
	if (new_node == NULL) { perror("Node creation."); exit(EXIT_FAILURE); }
	new_node->keys = allocate_align64((order - 1) * sizeof(uint64_t));
	if (new_node->keys == NULL) { perror("New node keys array."); exit(EXIT_FAILURE); }
	new_node->pointers = allocate_align64(order * sizeof(void *));
	if (new_node->pointers == NULL) { perror("New node pointers array."); exit(EXIT_FAILURE); }
	new_node->is_leaf  = false;
	new_node->num_keys = 0;
	new_node->parent   = NULL;
	new_node->next     = NULL;
	return new_node;
}

node *make_leaf(void) {
	node *leaf = make_node();
	leaf->is_leaf = true;
	return leaf;
}

uint64_t get_left_index(node *parent, node *left) {
	uint64_t left_index = 0;
	while (left_index <= parent->num_keys && parent->pointers[left_index] != left)
		left_index++;
	return left_index;
}

node *insert_into_leaf(node *leaf, uint64_t key, record *pointer) {
	uint64_t i, insertion_point = 0;
	while (insertion_point < leaf->num_keys && leaf->keys[insertion_point] < key)
		insertion_point++;
	for (i = leaf->num_keys; i > insertion_point; i--) {
		leaf->keys[i]     = leaf->keys[i - 1];
		leaf->pointers[i] = leaf->pointers[i - 1];
	}
	leaf->keys[insertion_point]     = key;
	leaf->pointers[insertion_point] = pointer;
	leaf->num_keys++;
	return leaf;
}

node *insert_into_leaf_after_splitting(node *root, node *leaf, uint64_t key, record *pointer) {
	node *new_leaf;
	uint64_t *temp_keys;
	void **temp_pointers;
	uint64_t insertion_index, split, new_key, i, j;

	new_leaf = make_leaf();
	temp_keys = malloc((order + 1) * sizeof(uint64_t));
	if (temp_keys == NULL) { perror("Temporary keys array."); exit(EXIT_FAILURE); }
	temp_pointers = malloc((order + 1) * sizeof(void *));
	if (temp_pointers == NULL) { perror("Temporary pointers array."); exit(EXIT_FAILURE); }

	insertion_index = 0;
	while (insertion_index < order - 1 && leaf->keys[insertion_index] < key)
		insertion_index++;
	for (i = 0, j = 0; i < leaf->num_keys; i++, j++) {
		if (j == insertion_index) j++;
		temp_keys[j]     = leaf->keys[i];
		temp_pointers[j] = leaf->pointers[i];
	}
	temp_keys[insertion_index]     = key;
	temp_pointers[insertion_index] = pointer;

	leaf->num_keys = 0;
	split = cut(order - 1);
	for (i = 0; i < split; i++) {
		leaf->pointers[i] = temp_pointers[i];
		leaf->keys[i]     = temp_keys[i];
		leaf->num_keys++;
	}
	for (i = split, j = 0; i < order; i++, j++) {
		new_leaf->pointers[j] = temp_pointers[i];
		new_leaf->keys[j]     = temp_keys[i];
		new_leaf->num_keys++;
	}
	free(temp_pointers);
	free(temp_keys);

	new_leaf->pointers[order - 1] = leaf->pointers[order - 1];
	leaf->pointers[order - 1]     = new_leaf;
	for (i = leaf->num_keys;     i < order - 1; i++) leaf->pointers[i]     = NULL;
	for (i = new_leaf->num_keys; i < order - 1; i++) new_leaf->pointers[i] = NULL;

	new_leaf->parent = leaf->parent;
	new_key = new_leaf->keys[0];
	return insert_into_parent(root, leaf, new_key, new_leaf);
}

node *insert_into_node(node *root, node *n, uint64_t left_index, uint64_t key, node *right) {
	uint64_t i;
	for (i = n->num_keys; i > left_index; i--) {
		n->pointers[i + 1] = n->pointers[i];
		n->keys[i]         = n->keys[i - 1];
	}
	n->pointers[left_index + 1] = right;
	n->keys[left_index]         = key;
	n->num_keys++;
	return root;
}

node *insert_into_node_after_splitting(node *root, node *old_node, uint64_t left_index,
                                       uint64_t key, node *right) {
	uint64_t i, j, split, k_prime;
	node *new_node, *child;
	uint64_t *temp_keys;
	node **temp_pointers;

	temp_pointers = malloc((order + 1) * sizeof(node *));
	if (temp_pointers == NULL) { perror("Temporary pointers array for splitting nodes."); exit(EXIT_FAILURE); }
	temp_keys = malloc(order * sizeof(uint64_t));
	if (temp_keys == NULL) { perror("Temporary keys array for splitting nodes."); exit(EXIT_FAILURE); }

	for (i = 0, j = 0; i < old_node->num_keys + 1; i++, j++) {
		if (j == left_index + 1) j++;
		temp_pointers[j] = old_node->pointers[i];
	}
	for (i = 0, j = 0; i < old_node->num_keys; i++, j++) {
		if (j == left_index) j++;
		temp_keys[j] = old_node->keys[i];
	}
	temp_pointers[left_index + 1] = right;
	temp_keys[left_index]         = key;

	split = cut(order);
	new_node = make_node();
	old_node->num_keys = 0;
	for (i = 0; i < split - 1; i++) {
		old_node->pointers[i] = temp_pointers[i];
		old_node->keys[i]     = temp_keys[i];
		old_node->num_keys++;
	}
	old_node->pointers[i] = temp_pointers[i];
	k_prime = temp_keys[split - 1];
	for (++i, j = 0; i < order; i++, j++) {
		new_node->pointers[j] = temp_pointers[i];
		new_node->keys[j]     = temp_keys[i];
		new_node->num_keys++;
	}
	new_node->pointers[j] = temp_pointers[i];
	free(temp_pointers);
	free(temp_keys);

	new_node->parent = old_node->parent;
	for (i = 0; i <= new_node->num_keys; i++) {
		child = new_node->pointers[i];
		child->parent = new_node;
	}
	return insert_into_parent(root, old_node, k_prime, new_node);
}

node *insert_into_parent(node *root, node *left, uint64_t key, node *right) {
	uint64_t left_index;
	node *parent = left->parent;
	if (parent == NULL) return insert_into_new_root(left, key, right);
	left_index = get_left_index(parent, left);
	if (parent->num_keys < order - 1)
		return insert_into_node(root, parent, left_index, key, right);
	return insert_into_node_after_splitting(root, parent, left_index, key, right);
}

node *insert_into_new_root(node *left, uint64_t key, node *right) {
	node *root = make_node();
	root->keys[0]     = key;
	root->pointers[0] = left;
	root->pointers[1] = right;
	root->num_keys++;
	root->parent  = NULL;
	left->parent  = root;
	right->parent = root;
	return root;
}

node *start_new_tree(uint64_t key, record *pointer) {
	node *root = make_leaf();
	root->keys[0]             = key;
	root->pointers[0]         = pointer;
	root->pointers[order - 1] = NULL;
	root->parent              = NULL;
	root->num_keys++;
	return root;
}

node *insert(node *root, uint64_t key, uint64_t value) {
	record *record_pointer = find(root, key, false, NULL);
	if (record_pointer != NULL) { record_pointer->value = value; return root; }
	record_pointer = make_record(value);
	if (root == NULL) return start_new_tree(key, record_pointer);
	node *leaf = find_leaf(root, key, false);
	if (leaf->num_keys < order - 1) {
		leaf = insert_into_leaf(leaf, key, record_pointer);
		return root;
	}
	return insert_into_leaf_after_splitting(root, leaf, key, record_pointer);
}

uint64_t get_neighbor_index(node *n) {
	uint64_t i;
	for (i = 0; i <= n->parent->num_keys; i++)
		if (n->parent->pointers[i] == n) return i - 1;
	printf("Search for nonexistent pointer to node in parent.\n");
	printf("Node:  %#lx\n", (unsigned long)n);
	exit(EXIT_FAILURE);
}

node *remove_entry_from_node(node *n, uint64_t key, node *pointer) {
	uint64_t i, num_pointers;
	i = 0;
	while (n->keys[i] != key) i++;
	for (++i; i < n->num_keys; i++) n->keys[i - 1] = n->keys[i];
	num_pointers = n->is_leaf ? n->num_keys : n->num_keys + 1;
	i = 0;
	while (n->pointers[i] != pointer) i++;
	for (++i; i < num_pointers; i++) n->pointers[i - 1] = n->pointers[i];
	n->num_keys--;
	if (n->is_leaf)
		for (i = n->num_keys; i < order - 1; i++) n->pointers[i] = NULL;
	else
		for (i = n->num_keys + 1; i < order; i++) n->pointers[i] = NULL;
	return n;
}

node *adjust_root(node *root) {
	node *new_root;
	if (root->num_keys > 0) return root;
	if (!root->is_leaf) { new_root = root->pointers[0]; new_root->parent = NULL; }
	else                  new_root = NULL;
	free(root->keys);
	free(root->pointers);
	free_node(root);
	return new_root;
}

node *coalesce_nodes(node *root, node *n, node *neighbor, uint64_t neighbor_index, uint64_t k_prime) {
	uint64_t i, j, neighbor_insertion_index, n_end;
	node *tmp;

	if (neighbor_index == -1) { tmp = n; n = neighbor; neighbor = tmp; }

	neighbor_insertion_index = neighbor->num_keys;

	if (!n->is_leaf) {
		neighbor->keys[neighbor_insertion_index] = k_prime;
		neighbor->num_keys++;
		n_end = n->num_keys;
		for (i = neighbor_insertion_index + 1, j = 0; j < n_end; i++, j++) {
			neighbor->keys[i]     = n->keys[j];
			neighbor->pointers[i] = n->pointers[j];
			neighbor->num_keys++;
			n->num_keys--;
		}
		neighbor->pointers[i] = n->pointers[j];
		for (i = 0; i < neighbor->num_keys + 1; i++) {
			tmp = (node *)neighbor->pointers[i];
			tmp->parent = neighbor;
		}
	} else {
		for (i = neighbor_insertion_index, j = 0; j < n->num_keys; i++, j++) {
			neighbor->keys[i]     = n->keys[j];
			neighbor->pointers[i] = n->pointers[j];
			neighbor->num_keys++;
		}
		neighbor->pointers[order - 1] = n->pointers[order - 1];
	}
	root = delete_entry(root, n->parent, k_prime, n);
	free(n->keys);
	free(n->pointers);
	free_node(n);
	return root;
}

node *redistribute_nodes(node *root, node *n, node *neighbor, uint64_t neighbor_index,
                         uint64_t k_prime_index, uint64_t k_prime) {
	uint64_t i;
	node *tmp;

	if (neighbor_index != -1) {
		if (!n->is_leaf) n->pointers[n->num_keys + 1] = n->pointers[n->num_keys];
		for (i = n->num_keys; i > 0; i--) {
			n->keys[i]     = n->keys[i - 1];
			n->pointers[i] = n->pointers[i - 1];
		}
		if (!n->is_leaf) {
			n->pointers[0] = neighbor->pointers[neighbor->num_keys];
			tmp = (node *)n->pointers[0];
			tmp->parent = n;
			neighbor->pointers[neighbor->num_keys] = NULL;
			n->keys[0] = k_prime;
			n->parent->keys[k_prime_index] = neighbor->keys[neighbor->num_keys - 1];
		} else {
			n->pointers[0] = neighbor->pointers[neighbor->num_keys - 1];
			neighbor->pointers[neighbor->num_keys - 1] = NULL;
			n->keys[0] = neighbor->keys[neighbor->num_keys - 1];
			n->parent->keys[k_prime_index] = n->keys[0];
		}
	} else {
		if (n->is_leaf) {
			n->keys[n->num_keys]     = neighbor->keys[0];
			n->pointers[n->num_keys] = neighbor->pointers[0];
			n->parent->keys[k_prime_index] = neighbor->keys[1];
		} else {
			n->keys[n->num_keys]         = k_prime;
			n->pointers[n->num_keys + 1] = neighbor->pointers[0];
			tmp = (node *)n->pointers[n->num_keys + 1];
			tmp->parent = n;
			n->parent->keys[k_prime_index] = neighbor->keys[0];
		}
		for (i = 0; i < neighbor->num_keys - 1; i++) {
			neighbor->keys[i]     = neighbor->keys[i + 1];
			neighbor->pointers[i] = neighbor->pointers[i + 1];
		}
		if (!n->is_leaf) neighbor->pointers[i] = neighbor->pointers[i + 1];
	}
	n->num_keys++;
	neighbor->num_keys--;
	return root;
}

node *delete_entry(node *root, node *n, uint64_t key, void *pointer) {
	uint64_t min_keys, neighbor_index, k_prime_index, k_prime, capacity;
	node *neighbor;

	n = remove_entry_from_node(n, key, pointer);
	if (n == root) return adjust_root(root);

	min_keys = n->is_leaf ? cut(order - 1) : cut(order) - 1;
	if (n->num_keys >= min_keys) return root;

	neighbor_index = get_neighbor_index(n);
	k_prime_index  = neighbor_index == -1 ? 0 : neighbor_index;
	k_prime        = n->parent->keys[k_prime_index];
	neighbor       = neighbor_index == -1 ? n->parent->pointers[1] :
	                                        n->parent->pointers[neighbor_index];
	capacity = n->is_leaf ? order : order - 1;

	if (neighbor->num_keys + n->num_keys < capacity)
		return coalesce_nodes(root, n, neighbor, neighbor_index, k_prime);
	return redistribute_nodes(root, n, neighbor, neighbor_index, k_prime_index, k_prime);
}

node *delete(node *root, uint64_t key) {
	node   *key_leaf   = NULL;
	record *key_record = find(root, key, false, &key_leaf);
	if (key_record != NULL && key_leaf != NULL) {
		root = delete_entry(root, key_leaf, key, key_record);
		free_record(key_record);
	}
	return root;
}

void destroy_tree_nodes(node *root) {
	uint64_t i;
	if (root->is_leaf)
		for (i = 0; i < root->num_keys; i++) free(root->pointers[i]);
	else
		for (i = 0; i < root->num_keys + 1; i++) destroy_tree_nodes(root->pointers[i]);
	free(root->pointers);
	free(root->keys);
	free_node(root);
}

node *destroy_tree(node *root) {
	destroy_tree_nodes(root);
	return NULL;
}

#define N       16
#define MASK    ((1 << (N - 1)) + (1 << (N - 1)) - 1)
#define LOW(x)  ((unsigned)(x) & MASK)
#define HIGH(x) LOW((x) >> N)
#define MUL(x, y, z)  { int32_t l = (long)(x) * (long)(y); \
	(z)[0] = LOW(l); (z)[1] = HIGH(l); }
#define CARRY(x, y)   ((int32_t)(x) + (long)(y) > MASK)
#define ADDEQU(x, y, z) (z = CARRY(x, (y)), x = LOW(x + (y)))
#define X0  0x330E
#define X1  0xABCD
#define X2  0x1234
#define A0  0xE66D
#define A1  0xDEEC
#define A2  0x5
#define C   0xB
#define SET3(x, x0, x1, x2) ((x)[0] = (x0), (x)[1] = (x1), (x)[2] = (x2))
#define SEED(x0, x1, x2) (SET3(x, x0, x1, x2), SET3(a, A0, A1, A2), c = C)

static uint64_t x[3] = { X0, X1, X2 }, a[3] = { A0, A1, A2 }, c = C;
static void next(void);

uint64_t redisLrand48(void) {
	next();
	return (((uint64_t)x[2] << (N - 1)) + (x[1] >> 1));
}

void redisSrand48(int32_t seedval) {
	SEED(X0, LOW(seedval), HIGH(seedval));
}

static void next(void) {
	uint64_t p[2], q[2], r[2], carry0, carry1;
	MUL(a[0], x[0], p);
	ADDEQU(p[0], c, carry0);
	ADDEQU(p[1], carry0, carry1);
	MUL(a[0], x[1], q);
	ADDEQU(p[1], q[0], carry0);
	MUL(a[1], x[0], r);
	x[2] = LOW(carry0 + carry1 + CARRY(p[1], r[0]) + q[1] + r[1] +
	           a[0] * x[2] + a[1] * x[1] + a[2] * x[0]);
	x[1] = LOW(p[1] + r[0]);
	x[0] = LOW(p[0]);
}

static size_t nelements_from_memory(size_t target_bytes, int tree_order)
{
	double S_e = (double)sizeof(struct element);
	double S_r = (double)sizeof(record);
	double S_n = (double)(sizeof(node)
	             + (size_t)(tree_order - 1) * sizeof(uint64_t)
	             + (size_t) tree_order      * sizeof(void *));
	double leaf_fill   = (double)(((tree_order - 1) + 1) / 2);
	double node_factor = S_n * tree_order / (leaf_fill * (tree_order - 1));
	double per_element = S_e + S_r + node_factor;
	return (size_t)((double)target_bytes / per_element);
}

static void print_memory_prediction(size_t nelements, int tree_order)
{
	size_t elem_bytes = nelements * sizeof(struct element);
	size_t rec_bytes  = nelements * sizeof(record);
	size_t leaf_fill  = (size_t)(((tree_order - 1) + 1) / 2);
	size_t num_leaves = (nelements + leaf_fill - 1) / leaf_fill;
	size_t num_nodes  = num_leaves * (size_t)tree_order / (size_t)(tree_order - 1);
	size_t node_bytes = num_nodes * (sizeof(node)
	                    + (size_t)(tree_order - 1) * sizeof(uint64_t)
	                    + (size_t) tree_order      * sizeof(void *));
	size_t total = elem_bytes + rec_bytes + node_bytes;
	printf("Predicted memory: %zu MB  (elements: %zu MB, records: %zu MB, nodes: %zu MB)\n",
	       total >> 20, elem_bytes >> 20, rec_bytes >> 20, node_bytes >> 20);
}

static void print_usage(const char *prog) {
	fprintf(stderr, "Usage: %s [--nelements <n>] [--nlookup <n>] [--order <n>]\n", prog);
	fprintf(stderr, "  --nelements <n>  number of elements to insert        (default: %zu)\n", DEFAULT_NELEMENTS);
	fprintf(stderr, "  --nlookup   <n>  number of lookup iterations         (default: %zu)\n", DEFAULT_NLOOKUP);
	fprintf(stderr, "  --order     <n>  B+ tree order, %d-%d                (default: %d)\n",
	        MIN_ORDER, MAX_ORDER, DEFAULT_TREE_ORDER);
}

int main(int argc, char **argv)
{
	size_t nelements  = DEFAULT_NELEMENTS;
	size_t nlookup    = DEFAULT_NLOOKUP;
	int    tree_order = DEFAULT_TREE_ORDER;

	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "--nelements") == 0 && i + 1 < argc) {
			nelements = (size_t)strtoull(argv[++i], NULL, 0);
		} else if (strcmp(argv[i], "--nlookup") == 0 && i + 1 < argc) {
			nlookup = (size_t)strtoull(argv[++i], NULL, 0);
		} else if (strcmp(argv[i], "--order") == 0 && i + 1 < argc) {
			tree_order = (int)strtol(argv[++i], NULL, 0);
			if (tree_order < MIN_ORDER || tree_order > MAX_ORDER) {
				fprintf(stderr, "order must be between %d and %d\n", MIN_ORDER, MAX_ORDER);
				return EXIT_FAILURE;
			}
		} else if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
			print_usage(argv[0]);
			return EXIT_SUCCESS;
		} else {
			fprintf(stderr, "Unknown argument: %s\n", argv[i]);
			print_usage(argv[0]);
			return EXIT_FAILURE;
		}
	}

	order = tree_order;

	printf("BTree Elements: %zuM\n", nelements / 1000000);
	printf("BTree #Lookups: %zuM\n", nlookup   / 1000000);
	print_memory_prediction(nelements, tree_order);

	node *root = NULL;
	verbose_output = false;

	redisSrand48(0xcafebabe);

	struct element *elms  = allocate((nelements / 2) * sizeof(struct element));
	struct element *elms2 = allocate((nelements / 2) * sizeof(struct element));

	for (size_t i = 0; i < nelements / 2; i++) {
		elms[i].payload  = i + 1;
		elms[i].payload2 = i + 1;
		elms2[i].payload  = i + 1;
		elms2[i].payload2 = i + 1;
	}

	for (size_t i = 0; i < nelements; i += 2) {
		root = insert(root, i,                 (uint64_t)&elms[i / 2]);
		root = insert(root, nelements - i - 1, (uint64_t)&elms2[i / 2]);
	}

	printf("Actual memory:   %zu MB\n", allocator_stat >> 20);

	uint64_t sum = 0;
	uint64_t hits = 0;
	struct timeval start, end;
	gettimeofday(&start, NULL);

#ifdef _OPENMP
#pragma omp parallel for reduction(+:sum) reduction(+:hits)
#endif
	for (size_t i = 0; i < nlookup; i++) {
		size_t rdn = (size_t)(i * 6364136223846793005ULL + 1442695040888963407ULL);
		record *r = find(root, rdn % nelements, false, NULL);
		if (r) {
			struct element *e = (struct element *)r->value;
			if (e) { sum += e->payload; hits++; }
		}
		r = find(root, ((rdn + 1) << 2) % nelements, false, NULL);
		if (r) {
			struct element *e = (struct element *)r->value;
			if (e) { sum += e->payload; hits++; }
		}
	}

	gettimeofday(&end, NULL);
	printf("hits: %zu  sum: %zu  time: %zu seconds\n", hits, sum, end.tv_sec - start.tv_sec);

	return EXIT_SUCCESS;
}
