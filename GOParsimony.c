/******************************************************************/
/*                                                                */ 
/* GOParsimony.c, a program to analyze Gene Ontology annotations  */  
/* parsimoniously, assigning precedence to terms in the digraphs  */ 
/*                                                                */ 
/* Copyright 2026. Paul Martin Harrison                           */ 
/* Licensed with a 3-clause BSD license.                          */ 
/*                                                                */ 
/* To compile: gcc -O2 -o GOParsimony GOParsimony.c -lm           */ 
/*                                                                */ 
/* To run and get help: ./GOParsimony -h                          */ 
/*                                                                */ 
/* This software is part of the GOParsimony package available on  */ 
/* Github.                                                        */ 
/*                                                                */ 
/******************************************************************/ 

#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <ctype.h> 
#include <unistd.h> 
#include <math.h> 

// --- Constants and Global Map Simulation ---   

#define MAX_TERM_ID_LEN 15
#define MAX_NAME_LEN 256
#define MAX_NAMESPACE_LEN 64
#define MAX_REL_TYPE_LEN 64
#define MAX_DEF_LEN 2048 
#define MAX_TERMS 60000 
#define INITIAL_CAPACITY 4 

typedef struct GoTerm GoTerm; 

typedef struct {
    char id[MAX_TERM_ID_LEN];
    GoTerm *term;
} MapEntry;

MapEntry term_map[MAX_TERMS];
int term_count = 0;
int global_session_id = 0; 
int global_cluster_id = 1;

// --- Data Structures ---

typedef struct Relationship {
    char type[MAX_REL_TYPE_LEN];    
    char target_id[MAX_TERM_ID_LEN]; 
    GoTerm *target_ptr;             
} Relationship;

struct GoTerm {
    char id[MAX_TERM_ID_LEN];
    char name[MAX_NAME_LEN];
    char namespace[MAX_NAMESPACE_LEN];
    char definition[MAX_DEF_LEN];
    
    Relationship *parents;
    int num_parents;
    int parents_capacity;

    Relationship *children;
    int num_children;
    int children_capacity;

    Relationship *reg_parents;
    int num_reg_parents;
    int reg_parents_capacity;

    Relationship *reg_children;
    int num_reg_children;
    int reg_children_capacity;

    int visited_marker; 
    bool is_slim; 
    
    int bg_count;
    int sample_count;

    void *enrich_res;
    void *deplet_res;
};

typedef struct {
    char prot_id[64];
    char go_id[MAX_TERM_ID_LEN];
} Annotation;

typedef struct EnrichmentResult {
    GoTerm *term;
    int k_count;          
    int K_count;          
    double expected_val;
    double p_value;
    double adj_p_value; 
    char res_type[16];  
    bool is_significant; 
    
    int precedence[2]; 
    bool prec_locked[2];  
    bool has_equal_pval_ancestor[2]; 
    int cluster_id[2]; 

    GoTerm **ancestors[2];
    int num_ancestors[2];
    GoTerm **descendants[2];
    int num_descendants[2];
} EnrichmentResult;

// --- Function Prototypes ---
void print_help();
void build_graph(char *obo_content);
void parse_goslim(char *slim_content);
void link_children();
void cleanup_graph();
GoTerm* get_term(const char *term_id);
void insert_term(GoTerm *term);
void trim_whitespace(char *str);
void init_relationship_arrays(GoTerm *term);
void grow_and_add_relationship(Relationship **array_ptr, int *count_ptr, int *capacity_ptr, const char *type, const char *target_id, GoTerm *target_ptr);
char* read_file_to_buffer(const char *filename);
bool is_regulatory(const char *type);

Annotation* load_background_annotations(const char *filename, int *out_count);
Annotation* load_sample_annotations(const char *filename, Annotation *bg_annots, int bg_count, int *out_count);
int expand_and_write(Annotation *annots, int count, const char *filename, const char *file_label, bool is_bg, const char *prefix, char *out_expanded_name);

int cmp_annot(const void *a, const void *b);
void expand_annotation(GoTerm *term, FILE *out, const char *prot_id, bool is_bg, int session_id);

double log_choose(int n, int k);
double hypergeom_pmf(int k, int N, int K, int n);
double hypergeom_pval_enrich(int k, int N, int K, int n);
double hypergeom_pval_deplet(int k, int N, int K, int n);
void calculate_enrichment(int bg_N, int sample_n, double alpha_threshold, bool use_bonferroni, const char *prefix, bool concise_mode, int concise_limit, bool include_depletions, const char *sample_expanded_fname);
int cmp_results_pval(const void *a, const void *b);

const char* get_short_ns(const char* ns);
void calculate_precedence(EnrichmentResult *list, int total_items, const char *target_type, const char *target_ns, int mode);
void assign_clusters(EnrichmentResult *list, int total_items, const char *target_type, const char *target_ns, int mode);
void print_cluster_graph(EnrichmentResult *node, int cluster_id, int depth, int *offsets, bool is_first_child, FILE *out, const char *target_type, int mode);

void collect_ancestors(GoTerm *t, GoTerm ***arr, int *cnt, int *cap, int session, int mode);
void collect_descendants(GoTerm *t, GoTerm ***arr, int *cnt, int *cap, int session, int mode);
int cmp_goterm_ptrs(const void *a, const void *b);

// --- Help Statement ---
void print_help() {
    printf("\n====================================================================================================\n");
    printf("         GOParsimony: parsimonious saltation for efficient Gene Ontology Enrichment Analysis \n");
    printf("====================================================================================================\n");
    printf("\nUsage: ./GOParsimony [OPTIONS]\n\n");

    printf("\nRequired Options:\n");
    printf("  -b <file>   Background population annotations (Simple 2-col or GAF)\n");
    printf("  -s <file>   Sample population IDs (single column list of Protein IDs or accessions)\n");

    printf("\nOntology Options:\n");
    printf("  -o <file>   Specify the OBO ontology file (Default: go-basic.obo)\n");
    printf("  -q <id>     Query a specific GO term ID (e.g., GO:0006915)\n");

    printf("\nStatistical Correction Options:\n");
    printf("  -d <float>  Set the threshold for significance (default: 0.05)\n");
    printf("              (Acts as FDR for BH procedure, or Alpha for Bonferroni)\n");
    printf("  -B          Use Bonferroni correction instead of Benjamini-Hochberg\n");
    printf("  -D          Enable analysis and output of Depletions (default is Enrichments only)\n");

    printf("\nMiscellaneous:\n");
    printf("  -e          Output the expanded annotation lists to files (default: do not output)\n");
    printf("  -p <string> Add a prefix to all output files (e.g., -p my_run)\n");
    printf("  -c [int]    Concise mode: Output ONLY the structural table, filtered to Precedence <= [int] (default: 1)\n");
    printf("  -m <file>   Specify a GO slims OBO file to tag specific terms\n");
    printf("  -h          Print this help statement and exit\n");
    printf("=======================================================================\n\n");
}

// --- Hypergeometric Math Functions ---

double log_choose(int n, int k) {
    if (k < 0 || k > n) return -INFINITY; 
    return lgamma(n + 1.0) - lgamma(k + 1.0) - lgamma(n - k + 1.0);
}

double hypergeom_pmf(int k, int N, int K, int n) {
    if (k < 0 || k > K || k > n || (n - k) > (N - K)) return 0.0;
    double lp = log_choose(K, k) + log_choose(N - K, n - k) - log_choose(N, n);
    return exp(lp);
}

double hypergeom_pval_enrich(int k, int N, int K, int n) {
    double p = 0.0;
    int max_k = (n < K) ? n : K;
    for (int idx = k; idx <= max_k; idx++) {
        p += hypergeom_pmf(idx, N, K, n);
    }
    return (p > 1.0) ? 1.0 : p;
}

double hypergeom_pval_deplet(int k, int N, int K, int n) {
    double p = 0.0;
    int min_k = (0 > n - (N - K)) ? 0 : n - (N - K);
    for (int idx = min_k; idx <= k; idx++) {
        p += hypergeom_pmf(idx, N, K, n);
    }
    return (p > 1.0) ? 1.0 : p;
}

// --- Helper Functions ---

const char* get_short_ns(const char* ns) {
    if (strcmp(ns, "biological_process") == 0) return "BP";
    if (strcmp(ns, "cellular_component") == 0) return "CC";
    if (strcmp(ns, "molecular_function") == 0) return "MF";
    return "UN";
}

bool is_regulatory(const char *type) {
    return (strcmp(type, "regulates") == 0 || 
            strcmp(type, "negatively_regulates") == 0 || 
            strcmp(type, "positively_regulates") == 0);
}

void grow_and_add_relationship(Relationship **array_ptr, int *count_ptr, int *capacity_ptr, 
                               const char *type, const char *target_id, GoTerm *target_ptr) {
    if (*count_ptr >= *capacity_ptr) {
        int new_capacity = (*capacity_ptr == 0) ? INITIAL_CAPACITY : (*capacity_ptr * 2);
        Relationship *new_array = (Relationship *)realloc(*array_ptr, new_capacity * sizeof(Relationship));
        if (!new_array) exit(EXIT_FAILURE);
        *array_ptr = new_array;
        *capacity_ptr = new_capacity;
    }
    Relationship *new_rel = &((*array_ptr)[*count_ptr]);
    strncpy(new_rel->type, type, MAX_REL_TYPE_LEN - 1);
    new_rel->type[MAX_REL_TYPE_LEN - 1] = '\0';
    strncpy(new_rel->target_id, target_id, MAX_TERM_ID_LEN - 1);
    new_rel->target_id[MAX_TERM_ID_LEN - 1] = '\0';
    new_rel->target_ptr = target_ptr; 
    (*count_ptr)++;
}

void init_relationship_arrays(GoTerm *term) {
    term->parents = NULL; term->num_parents = 0; term->parents_capacity = 0;
    term->children = NULL; term->num_children = 0; term->children_capacity = 0;
    term->reg_parents = NULL; term->num_reg_parents = 0; term->reg_parents_capacity = 0;
    term->reg_children = NULL; term->num_reg_children = 0; term->reg_children_capacity = 0;
    
    term->visited_marker = 0;
    term->is_slim = false;
    term->bg_count = 0;
    term->sample_count = 0;
    term->enrich_res = NULL;
    term->deplet_res = NULL;
    term->id[0] = term->name[0] = term->namespace[0] = term->definition[0] = '\0';
}

void insert_term(GoTerm *term) {
    if (term_count >= MAX_TERMS) return;
    strncpy(term_map[term_count].id, term->id, MAX_TERM_ID_LEN - 1);
    term_map[term_count].id[MAX_TERM_ID_LEN - 1] = '\0';
    term_map[term_count].term = term;
    term_count++;
}

GoTerm* get_term(const char *term_id) {
    for (int i = 0; i < term_count; i++) {
        if (strcmp(term_map[i].id, term_id) == 0) return term_map[i].term;
    }
    return NULL;
}

void trim_whitespace(char *str) {
    if (!str) return;
    char *start = str;
    while(isspace((unsigned char)*start)) start++;
    if(*start == 0) { *str = '\0'; return; }
    char *end = start + strlen(start) - 1;
    while(end > start && isspace((unsigned char)*end)) end--;
    *(end+1) = '\0';
    if (start > str) memmove(str, start, end - start + 2);
}

char* read_file_to_buffer(const char *filename) {
    FILE *f = fopen(filename, "rb");
    if (!f) { perror("Could not open file"); return NULL; }
    fseek(f, 0, SEEK_END);
    long length = ftell(f);
    fseek(f, 0, SEEK_SET);
    char *buffer = malloc(length + 1);
    if (buffer) {
        size_t read_len = fread(buffer, 1, length, f);
        buffer[read_len] = '\0';
    }
    fclose(f);
    return buffer;
}

// --- Query Traversal Helpers ---

int cmp_goterm_ptrs(const void *a, const void *b) {
    GoTerm *ta = *(GoTerm **)a;
    GoTerm *tb = *(GoTerm **)b;
    return strcmp(ta->id, tb->id);
}

void collect_ancestors(GoTerm *t, GoTerm ***arr, int *cnt, int *cap, int session, int mode) {
    if (!t) return;
    int n_parents = (mode == 0) ? t->num_parents : t->num_reg_parents;
    Relationship *rels = (mode == 0) ? t->parents : t->reg_parents;
    for(int i = 0; i < n_parents; i++) {
        GoTerm *p = rels[i].target_ptr;
        if(p && p->visited_marker != session) {
            p->visited_marker = session;
            if(*cnt >= *cap) {
                *cap = (*cap == 0) ? 16 : (*cap * 2);
                *arr = realloc(*arr, *cap * sizeof(GoTerm*));
            }
            (*arr)[(*cnt)++] = p;
            collect_ancestors(p, arr, cnt, cap, session, mode);
        }
    }
}

void collect_descendants(GoTerm *t, GoTerm ***arr, int *cnt, int *cap, int session, int mode) {
    if (!t) return;
    int n_children = (mode == 0) ? t->num_children : t->num_reg_children;
    Relationship *rels = (mode == 0) ? t->children : t->reg_children;
    for(int i = 0; i < n_children; i++) {
        GoTerm *c = rels[i].target_ptr;
        if(c && c->visited_marker != session) {
            c->visited_marker = session;
            if(*cnt >= *cap) {
                *cap = (*cap == 0) ? 16 : (*cap * 2);
                *arr = realloc(*arr, *cap * sizeof(GoTerm*));
            }
            (*arr)[(*cnt)++] = c;
            collect_descendants(c, arr, cnt, cap, session, mode);
        }
    }
}

// --- Core Logic ---

void build_graph(char *obo_content) {
    printf("Starting OBO file parsing...\n");
    char *content_copy = strdup(obo_content);
    if (!content_copy) return;

    char *line = strtok(content_copy, "\n");
    GoTerm *current_term = NULL;
    bool in_term_block = false;

    while (line != NULL) {
        trim_whitespace(line);

        if (line[0] == '[') {
            if (current_term) { insert_term(current_term); current_term = NULL; }
            if (strcmp(line, "[Term]") == 0) {
                in_term_block = true;
                current_term = (GoTerm *)malloc(sizeof(GoTerm));
                init_relationship_arrays(current_term);
            } else { in_term_block = false; }
        }
        else if (in_term_block && current_term) {
            char *key_end = strstr(line, ": ");
            if (key_end) {
                char key[MAX_REL_TYPE_LEN], value[MAX_DEF_LEN];
                int key_len = key_end - line;
                strncpy(key, line, (key_len < MAX_REL_TYPE_LEN) ? key_len : MAX_REL_TYPE_LEN - 1);
                key[(key_len < MAX_REL_TYPE_LEN) ? key_len : MAX_REL_TYPE_LEN - 1] = '\0';
                
                char *value_raw = key_end + 2;
                if (strcmp(key, "def") == 0) {
                    char *def_end = strchr(value_raw, '[');
                    int val_len = def_end ? (def_end - value_raw) : strlen(value_raw);
                    strncpy(value, value_raw, (val_len < MAX_DEF_LEN) ? val_len : MAX_DEF_LEN - 1);
                    value[(val_len < MAX_DEF_LEN) ? val_len : MAX_DEF_LEN - 1] = '\0';
                } else {
                    strncpy(value, value_raw, MAX_DEF_LEN - 1);
                    value[MAX_DEF_LEN - 1] = '\0';
                }
                trim_whitespace(value);

                if (strcmp(key, "id") == 0) strncpy(current_term->id, value, MAX_TERM_ID_LEN - 1);
                else if (strcmp(key, "name") == 0) strncpy(current_term->name, value, MAX_NAME_LEN - 1);
                else if (strcmp(key, "namespace") == 0) strncpy(current_term->namespace, value, MAX_NAMESPACE_LEN - 1);
                else if (strcmp(key, "def") == 0) strncpy(current_term->definition, value, MAX_DEF_LEN - 1);
                else if (strcmp(key, "is_a") == 0) {
                    char *excl = strchr(value, '!');
                    if (excl) *excl = '\0';
                    trim_whitespace(value);
                    if (strlen(value) > 0) grow_and_add_relationship(&(current_term->parents), &(current_term->num_parents), &(current_term->parents_capacity), "is_a", value, NULL);
                } else if (strcmp(key, "relationship") == 0) {
                    char *excl = strchr(value, '!');
                    if (excl) *excl = '\0';
                    trim_whitespace(value);
                    char *space = strchr(value, ' ');
                    if (space) {
                        *space = '\0';
                        char *rel_type = value;
                        char *parent_id = space + 1;
                        trim_whitespace(parent_id);
                        if (strlen(rel_type) > 0 && strlen(parent_id) > 0) {
                            if (is_regulatory(rel_type)) grow_and_add_relationship(&(current_term->reg_parents), &(current_term->num_reg_parents), &(current_term->reg_parents_capacity), rel_type, parent_id, NULL);
                            else grow_and_add_relationship(&(current_term->parents), &(current_term->num_parents), &(current_term->parents_capacity), rel_type, parent_id, NULL);
                        }
                    }
                }
            }
        }
        line = strtok(NULL, "\n"); 
    }
    if (current_term) insert_term(current_term);
    free(content_copy);
    printf("Parsing complete. Found %d terms.\n", term_count);
}

void parse_goslim(char *slim_content) {
    printf("Parsing GO slims file to tag terms...\n");
    char *content_copy = strdup(slim_content);
    if (!content_copy) return;

    char *line = strtok(content_copy, "\n");
    bool in_term_block = false;

    while (line != NULL) {
        trim_whitespace(line);

        if (line[0] == '[') {
            if (strcmp(line, "[Term]") == 0) in_term_block = true;
            else in_term_block = false;
        } 
        else if (in_term_block) {
            if (strncmp(line, "id: ", 4) == 0) {
                char *id_val = line + 4;
                trim_whitespace(id_val);
                
                GoTerm *t = get_term(id_val);
                if (t) t->is_slim = true;
            }
        }
        line = strtok(NULL, "\n");
    }
    free(content_copy);
    printf("Slim parsing complete.\n");
}

void link_children() {
    for (int i = 0; i < term_count; i++) {
        GoTerm *term = term_map[i].term;
        for(int j = 0; j < term->num_parents; j++) {
            Relationship *parent_rel = &(term->parents[j]);
            GoTerm *parent_term = get_term(parent_rel->target_id);
            if (parent_term) {
                parent_rel->target_ptr = parent_term;
                grow_and_add_relationship(&(parent_term->children), &(parent_term->num_children), &(parent_term->children_capacity), parent_rel->type, term->id, term);
            }
        }
        for(int j = 0; j < term->num_reg_parents; j++) {
            Relationship *parent_rel = &(term->reg_parents[j]);
            GoTerm *parent_term = get_term(parent_rel->target_id);
            if (parent_term) {
                parent_rel->target_ptr = parent_term;
                grow_and_add_relationship(&(parent_term->reg_children), &(parent_term->num_reg_children), &(parent_term->reg_children_capacity), parent_rel->type, term->id, term);
            }
        }
    }
}

void cleanup_graph() {
    for (int i = 0; i < term_count; i++) {
        free(term_map[i].term->parents);
        free(term_map[i].term->children);
        free(term_map[i].term->reg_parents);
        free(term_map[i].term->reg_children);
        free(term_map[i].term);
    }
}

// --- Annotation Loading and Expansion Logic ---

int cmp_annot(const void *a, const void *b) {
    return strcmp(((Annotation*)a)->prot_id, ((Annotation*)b)->prot_id);
}

int expand_session_id = 0;

void expand_annotation(GoTerm *term, FILE *out, const char *prot_id, bool is_bg, int session_id) {
    if (!term) return;
    
    if (term->enrich_res == (void *)(long)session_id) return; 
    term->enrich_res = (void *)(long)session_id;

    fprintf(out, "%s\t%s\n", prot_id, term->id);
    
    if (is_bg) term->bg_count++;
    else term->sample_count++;

    for (int i = 0; i < term->num_parents; i++) {
        expand_annotation(term->parents[i].target_ptr, out, prot_id, is_bg, session_id);
    }
}

Annotation* load_background_annotations(const char *filename, int *out_count) {
    printf("\nProcessing background file: %s\n", filename);
    FILE *f = fopen(filename, "r");
    if (!f) { perror("Failed to open background file"); return NULL; }
    
    // --- Auto-Detect Format ---
    bool is_gaf = false;
    char line[4096];
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '!') continue; 
        int tab_count = 0;
        for (int i = 0; line[i] != '\0'; i++) {
            if (line[i] == '\t') tab_count++;
        }
        if (tab_count >= 10) {
            is_gaf = true;
        }
        break;
    }
    rewind(f); 
    
    if (is_gaf) {
        printf(" -> Auto-detected GAF format.\n");
    } else {
        printf(" -> Auto-detected simple 2-column format.\n");
    }

    int capacity = 10000;
    int count = 0;
    Annotation *annots = malloc(capacity * sizeof(Annotation));
    
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '!') continue;
        
        char prot[128] = {0};
        char go[MAX_TERM_ID_LEN] = {0};
        
        if (is_gaf) {
            char *col[20];
            int col_count = 0;
            char *ptr = line;
            col[col_count++] = ptr;
            
            while (*ptr != '\0' && *ptr != '\n' && col_count < 20) {
                if (*ptr == '\t') {
                    *ptr = '\0';
                    col[col_count++] = ptr + 1;
                }
                ptr++;
            }
            
            if (col_count >= 5) {
                if (strstr(col[3], "NOT") != NULL) continue;
                
                strncpy(prot, col[1], 127);
                prot[127] = '\0';
                
                strncpy(go, col[4], MAX_TERM_ID_LEN - 1);
                go[MAX_TERM_ID_LEN - 1] = '\0';
                
                trim_whitespace(prot);
                trim_whitespace(go);
                
                if (strlen(prot) == 0 || strlen(go) == 0) continue;
            } else {
                continue;
            }
        } else {
            trim_whitespace(line);
            if (strlen(line) == 0) continue;
            if (sscanf(line, "%127s %14s", prot, go) != 2) continue;
        }
        
        if (count >= capacity) {
            capacity *= 2;
            Annotation *temp = realloc(annots, capacity * sizeof(Annotation));
            if (!temp) { free(annots); fclose(f); return NULL; }
            annots = temp;
        }
        strncpy(annots[count].prot_id, prot, 63);
        annots[count].prot_id[63] = '\0';
        strncpy(annots[count].go_id, go, MAX_TERM_ID_LEN - 1);
        annots[count].go_id[MAX_TERM_ID_LEN - 1] = '\0';
        count++;
    }
    fclose(f);
    
    qsort(annots, count, sizeof(Annotation), cmp_annot);
    *out_count = count;
    return annots;
}

Annotation* load_sample_annotations(const char *filename, Annotation *bg_annots, int bg_count, int *out_count) {
    printf("\nProcessing sample IDs file: %s\n", filename);
    FILE *f = fopen(filename, "r");
    if (!f) { perror("Failed to open sample file"); return NULL; }

    int capacity = 10000;
    int count = 0;
    Annotation *annots = malloc(capacity * sizeof(Annotation));

    char line[1024];
    int missing_count = 0;
    while (fgets(line, sizeof(line), f)) {
        trim_whitespace(line);
        if (strlen(line) == 0) continue;
        
        char prot[128] = {0};
        if (sscanf(line, "%127s", prot) == 1) {
            // Find prot in bg_annots using binary search (find first occurrence)
            int idx = -1;
            int low = 0, high = bg_count - 1;
            while (low <= high) {
                int mid = low + (high - low) / 2;
                int cmp = strcmp(bg_annots[mid].prot_id, prot);
                if (cmp == 0) {
                    idx = mid;
                    high = mid - 1; 
                } else if (cmp < 0) {
                    low = mid + 1;
                } else {
                    high = mid - 1;
                }
            }

            if (idx != -1) {
                while (idx < bg_count && strcmp(bg_annots[idx].prot_id, prot) == 0) {
                    if (count >= capacity) {
                        capacity *= 2;
                        Annotation *temp = realloc(annots, capacity * sizeof(Annotation));
                        if (!temp) { free(annots); fclose(f); return NULL; }
                        annots = temp;
                    }
                    annots[count++] = bg_annots[idx];
                    idx++;
                }
            } else {
                missing_count++;
            }
        }
    }
    fclose(f);
    
    if (missing_count > 0) {
        printf(" -> Note: %d IDs in sample list were not found in the background annotations.\n", missing_count);
    }

    qsort(annots, count, sizeof(Annotation), cmp_annot);
    *out_count = count;
    return annots;
}

int expand_and_write(Annotation *annots, int count, const char *filename, const char *file_label, bool is_bg, const char *prefix, char *out_expanded_name) {
    const char *basename = strrchr(filename, '/');
    if (basename) basename++; else basename = filename;
    
    char base_no_ext[1024];
    strncpy(base_no_ext, basename, 1023);
    base_no_ext[1023] = '\0';
    char *ext = strrchr(base_no_ext, '.');
    if (ext && (strcmp(ext, ".gaf") == 0 || strcmp(ext, ".txt") == 0 || strcmp(ext, ".tsv") == 0)) {
        *ext = '\0';
    }
    
    if (prefix) {
        snprintf(out_expanded_name, 2048, "%s_Expanded_%s.txt", prefix, base_no_ext);
    } else {
        snprintf(out_expanded_name, 2048, "Expanded_%s.txt", base_no_ext);
    }
    
    FILE *out = fopen(out_expanded_name, "w");
    if (!out) { perror("Failed to create output file"); return 0; }
    
    int unique_proteins = 0;
    char current_prot[64] = "";
    
    for (int i = 0; i < count; i++) {
        if (strcmp(annots[i].prot_id, current_prot) != 0) {
            expand_session_id++;
            strcpy(current_prot, annots[i].prot_id);
            unique_proteins++;
        }
        
        GoTerm *term = get_term(annots[i].go_id);
        if (term) expand_annotation(term, out, current_prot, is_bg, expand_session_id);
    }
    
    fclose(out);
    printf("Successfully expanded %s list (%d unique items).\n", file_label, unique_proteins);
    
    return unique_proteins; 
}


// --- Central Precedence Logic mapping ---
void calculate_precedence(EnrichmentResult *list, int total_items, const char *target_type, const char *target_ns, int mode) {
    for (int i = 0; i < total_items; i++) {
        if (list[i].is_significant && strcmp(list[i].res_type, target_type) == 0 && strcmp(list[i].term->namespace, target_ns) == 0) {
            list[i].precedence[mode] = 1;
            list[i].prec_locked[mode] = false;
        }
    }

    for (int i = 0; i < total_items; i++) {   
        if (list[i].is_significant && strcmp(list[i].res_type, target_type) == 0 && strcmp(list[i].term->namespace, target_ns) == 0) {
            if(list[i].prec_locked[mode] == false) { 
                list[i].prec_locked[mode] = true;
             
                for(int j = 0; j < total_items; j++) {
                    if (i == j) continue;
                    if (list[j].is_significant && strcmp(list[j].res_type, target_type) == 0 && strcmp(list[j].term->namespace, target_ns) == 0) {

                        if (list[j].prec_locked[mode] == false) {
                            
                            bool is_relative = false;
                            
                            for (int p = 0; p < list[i].num_ancestors[mode]; p++) {
                                if (list[i].ancestors[mode][p] == list[j].term) {
                                    is_relative = true;
                                    break;
                                }
                            }
                            
                            if (!is_relative) {
                                for (int c = 0; c < list[i].num_descendants[mode]; c++) {
                                    if (list[i].descendants[mode][c] == list[j].term) {
                                        is_relative = true;
                                        break;
                                    }
                                }
                            }

                            if (is_relative) {
                                if (list[j].p_value > list[i].p_value) {
                                    if (list[j].precedence[mode] < list[i].precedence[mode] + 1) {
                                        list[j].precedence[mode] = list[i].precedence[mode] + 1;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    } 
    
    // Secondary sweep to catch and flag Equal P-Values (*)
    for (int i = 0; i < total_items; i++) {
        if (list[i].is_significant && strcmp(list[i].res_type, target_type) == 0 && strcmp(list[i].term->namespace, target_ns) == 0) {
            for (int p = 0; p < list[i].num_ancestors[mode]; p++) {
                GoTerm *anc_term = list[i].ancestors[mode][p];
                EnrichmentResult *anc_res = (strcmp(target_type, "Enr") == 0) ? (EnrichmentResult*)anc_term->enrich_res : (EnrichmentResult*)anc_term->deplet_res;
                if (anc_res && anc_res->is_significant && anc_res->precedence[mode] == list[i].precedence[mode]) {
                    if (anc_res->p_value == list[i].p_value) {
                        list[i].has_equal_pval_ancestor[mode] = true;
                        break;
                    }
                }
            }
        }
    }
} 

void assign_clusters(EnrichmentResult *list, int total_items, const char *target_type, const char *target_ns, int mode) {
    for (int i = 0; i < total_items; i++) {
        if (list[i].is_significant && strcmp(list[i].res_type, target_type) == 0 && strcmp(list[i].term->namespace, target_ns) == 0) {
            list[i].cluster_id[mode] = 0;
        }
    }
    for (int i = 0; i < total_items; i++) {
        if (list[i].is_significant && strcmp(list[i].res_type, target_type) == 0 && strcmp(list[i].term->namespace, target_ns) == 0 && list[i].cluster_id[mode] == 0) {
            
            int cid = global_cluster_id++;
            list[i].cluster_id[mode] = cid;
            
            EnrichmentResult **queue = malloc(total_items * sizeof(EnrichmentResult *));
            int head = 0, tail = 0;
            queue[tail++] = &list[i];

            while (head < tail) {
                EnrichmentResult *curr = queue[head++];
                
                int n_p = (mode == 0) ? curr->term->num_parents : curr->term->num_reg_parents;
                Relationship *p_rels = (mode == 0) ? curr->term->parents : curr->term->reg_parents;
                
                for (int p = 0; p < n_p; p++) {
                    EnrichmentResult *p_res = (strcmp(target_type, "Enr") == 0) ? (EnrichmentResult *)p_rels[p].target_ptr->enrich_res : (EnrichmentResult *)p_rels[p].target_ptr->deplet_res;
                    if (p_res && p_res->is_significant && p_res->cluster_id[mode] == 0 && strcmp(p_res->term->namespace, target_ns) == 0) {
                        p_res->cluster_id[mode] = cid;
                        queue[tail++] = p_res;
                    }
                }
                
                int n_c = (mode == 0) ? curr->term->num_children : curr->term->num_reg_children;
                Relationship *c_rels = (mode == 0) ? curr->term->children : curr->term->reg_children;
                
                for (int c = 0; c < n_c; c++) {
                    EnrichmentResult *c_res = (strcmp(target_type, "Enr") == 0) ? (EnrichmentResult *)c_rels[c].target_ptr->enrich_res : (EnrichmentResult *)c_rels[c].target_ptr->deplet_res;
                    if (c_res && c_res->is_significant && c_res->cluster_id[mode] == 0 && strcmp(c_res->term->namespace, target_ns) == 0) {
                        c_res->cluster_id[mode] = cid;
                        queue[tail++] = c_res;
                    }
                }
            }
            free(queue);
        }
    }
}

void print_cluster_graph(EnrichmentResult *node, int cluster_id, int depth, int *offsets, bool is_first_child, FILE *out, const char *target_type, int mode) {
    if (depth >= 99) return; 

    char buffer[512];
    int len;
    
    char prec_mark[16];
    if (node->has_equal_pval_ancestor[mode]) sprintf(prec_mark, "%d*", node->precedence[mode]);
    else sprintf(prec_mark, "%d", node->precedence[mode]);
    
    if (depth == 0) {
        len = sprintf(buffer, "%s (Prec %s)", node->term->id, prec_mark);
        if (out) fprintf(out, "%s", buffer);
        offsets[depth + 1] = offsets[depth] + len;
    } else {
        if (!is_first_child) {
            if (out) fprintf(out, "\n");
            for (int i = 0; i < offsets[depth]; i++) {
                if (out) fprintf(out, " ");
            }
        }
        len = sprintf(buffer, " --> %s (Prec %s)", node->term->id, prec_mark);
        if (out) fprintf(out, "%s", buffer);
        offsets[depth + 1] = offsets[depth] + len;
    }

    bool first = true;
    int n_c = (mode == 0) ? node->term->num_children : node->term->num_reg_children;
    Relationship *c_rels = (mode == 0) ? node->term->children : node->term->reg_children;
    
    for (int i = 0; i < n_c; i++) {
        GoTerm *child = c_rels[i].target_ptr;
        EnrichmentResult *c_res = (strcmp(target_type, "Enr") == 0) ? (EnrichmentResult*)child->enrich_res : (EnrichmentResult*)child->deplet_res;
        
        if (c_res && c_res->is_significant && c_res->cluster_id[mode] == cluster_id) {
            print_cluster_graph(c_res, cluster_id, depth + 1, offsets, first, out, target_type, mode);
            first = false;
        }
    }
}

// --- Statistics Output ---

int cmp_results_pval(const void *a, const void *b) {
    EnrichmentResult *r1 = (EnrichmentResult *)a;
    EnrichmentResult *r2 = (EnrichmentResult *)b;
    if (r1->p_value < r2->p_value) return -1;
    if (r1->p_value > r2->p_value) return 1;
    return 0;
}

void calculate_enrichment(int bg_N, int sample_n, double alpha_threshold, bool use_bonferroni, const char *prefix, bool concise_mode, int concise_limit, bool include_depletions, const char *sample_expanded_fname) {
    if (bg_N == 0 || sample_n == 0) {
        printf("\nError: Cannot calculate enrichment. Background or sample population is empty.\n");
        return;
    }

    // --- STEP 1: COUNT TOTAL TESTS FOR STATISTICAL CORRECTION ---
    int num_tests = 0;
    for (int i = 0; i < term_count; i++) {
        GoTerm *term = term_map[i].term;
        if (strcmp(term->id, "GO:0008150") == 0 || strcmp(term->id, "GO:0003674") == 0 || strcmp(term->id, "GO:0005575") == 0) continue;
        if (term->bg_count > 0) num_tests++;
    }

    if (num_tests == 0) {
        printf("\nNo valid tests to run after applying filters.\n");
        return;
    }

    EnrichmentResult *results = malloc(num_tests * sizeof(EnrichmentResult));
    if (!results) { perror("Memory fail for results"); return; }
    
    // --- STEP 2: CALCULATE STATISTICS FOR ALL HYPOTHESES ---
    int r_idx = 0;
    for (int i = 0; i < term_count; i++) {
        GoTerm *term = term_map[i].term;
        if (strcmp(term->id, "GO:0008150") == 0 || strcmp(term->id, "GO:0003674") == 0 || strcmp(term->id, "GO:0005575") == 0) continue;
        int k_val = term->sample_count;
        int K_val = term->bg_count;

        if (K_val == 0) continue; 

        double expected = (double)sample_n * K_val / bg_N;
        double p_val = 1.0;
        
        if (k_val >= expected) {
            p_val = hypergeom_pval_enrich(k_val, bg_N, K_val, sample_n);
            strcpy(results[r_idx].res_type, "Enr"); 
        } else {
            p_val = hypergeom_pval_deplet(k_val, bg_N, K_val, sample_n);
            strcpy(results[r_idx].res_type, "Dep"); 
        }

        results[r_idx].term = term;
        results[r_idx].k_count = k_val;
        results[r_idx].K_count = K_val;
        results[r_idx].expected_val = expected;
        results[r_idx].p_value = p_val;
        results[r_idx].is_significant = false;
        
        for(int m=0; m<2; m++) {
            results[r_idx].precedence[m] = 0; 
            results[r_idx].cluster_id[m] = 0;
            results[r_idx].has_equal_pval_ancestor[m] = false;
            results[r_idx].ancestors[m] = NULL;
            results[r_idx].num_ancestors[m] = 0;
            results[r_idx].descendants[m] = NULL;
            results[r_idx].num_descendants[m] = 0;
        }
        r_idx++;
    }

    // --- STEP 3: APPLY MULTIPLE HYPOTHESIS CORRECTIONS ---
    qsort(results, num_tests, sizeof(EnrichmentResult), cmp_results_pval); 

    if (use_bonferroni) {
        for (int i = 0; i < num_tests; i++) {
            double p_adj = results[i].p_value * num_tests;
            results[i].adj_p_value = (p_adj > 1.0) ? 1.0 : p_adj;
            if (results[i].adj_p_value <= alpha_threshold) results[i].is_significant = true;
        }
    } else {
        for (int i = 0; i < num_tests; i++) {
            double q = results[i].p_value * num_tests / (i + 1);
            results[i].adj_p_value = (q > 1.0) ? 1.0 : q;
        }
        for (int i = num_tests - 2; i >= 0; i--) {
            if (results[i+1].adj_p_value < results[i].adj_p_value) results[i].adj_p_value = results[i+1].adj_p_value;
        }
        for (int i = 0; i < num_tests; i++) {
            if (results[i].adj_p_value <= alpha_threshold) results[i].is_significant = true;
        }
    }
    
    // --- STEP 4: APPLY DEPLETION FILTERING (-D flag check) ---
    if (!include_depletions) {
        for (int i = 0; i < num_tests; i++) {
            if (strcmp(results[i].res_type, "Dep") == 0) results[i].is_significant = false;
        }
    }

    // --- POPULATE TOTAL ANCESTORS AND DESCENDANTS FOR SIGNIFICANT TERMS ---
    for (int i = 0; i < num_tests; i++) {
        if (results[i].is_significant) {
            for (int m=0; m<2; m++) {
                if (concise_mode && m == 1) continue; 
                int anc_cap = 0, desc_cap = 0;
                global_session_id++;
                collect_ancestors(results[i].term, &results[i].ancestors[m], &results[i].num_ancestors[m], &anc_cap, global_session_id, m);
                global_session_id++;
                collect_descendants(results[i].term, &results[i].descendants[m], &results[i].num_descendants[m], &desc_cap, global_session_id, m);
            }
        }
    }

    for (int i = 0; i < term_count; i++) {
        term_map[i].term->enrich_res = NULL;
        term_map[i].term->deplet_res = NULL;
    }
    for (int i = 0; i < num_tests; i++) {
        if (strcmp(results[i].res_type, "Enr") == 0) results[i].term->enrich_res = &results[i];
        else results[i].term->deplet_res = &results[i];
    }

    const char *ns_types[] = {"biological_process", "cellular_component", "molecular_function"};
    const char *res_types[] = {"Enr", "Dep"};
    const char *adj_p_label = use_bonferroni ? "Bonf_P" : "Q_Value";
    
    int num_res_types = include_depletions ? 2 : 1;

    // --- PERFORM CALCULATIONS FOR BOTH NETWORKS ---
    for (int mode = 0; mode < 2; mode++) {
        if (concise_mode && mode == 1) continue;

        for (int t = 0; t < num_res_types; t++) {
            for (int n = 0; n < 3; n++) {
                calculate_precedence(results, num_tests, res_types[t], ns_types[n], mode);
                assign_clusters(results, num_tests, res_types[t], ns_types[n], mode);
                
                if (mode == 1) {
                    for (int cid = 1; cid < global_cluster_id; cid++) {
                        int c_size = 0;
                        for (int k = 0; k < num_tests; k++) {
                            if (results[k].is_significant && strcmp(results[k].res_type, res_types[t]) == 0 && strcmp(results[k].term->namespace, ns_types[n]) == 0 && results[k].cluster_id[mode] == cid) {
                                c_size++;
                            }
                        }
                        if (c_size == 1) {
                            for (int k = 0; k < num_tests; k++) {
                                if (results[k].is_significant && strcmp(results[k].res_type, res_types[t]) == 0 && strcmp(results[k].term->namespace, ns_types[n]) == 0 && results[k].cluster_id[mode] == cid) {
                                    results[k].cluster_id[mode] = 0;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // --- OUTPUT PHASE ---
    char results_fname[2048], graphs_fname[2048], edges_fname[2048];
    if (prefix) {
        snprintf(results_fname, sizeof(results_fname), "%s_Enrichment_Results.tsv", prefix);
        snprintf(graphs_fname, sizeof(graphs_fname), "%s_Precedence_Graphs.txt", prefix);
        snprintf(edges_fname, sizeof(edges_fname), "%s_Precedence_Edges.tsv", prefix);
    } else {
        strcpy(results_fname, "Enrichment_Results.tsv");
        strcpy(graphs_fname, "Precedence_Graphs.txt");
        strcpy(edges_fname, "Precedence_Edges.tsv");
    }

    printf("\n==========================================================================================================================\n");
    printf("GO TERM STATISTICAL ANALYSIS\n");
    if (use_bonferroni) printf("(Bonferroni Correction)\n");
    else printf("(Benjamini-Hochberg FDR)\n");
    if (!include_depletions) printf("** ENRICHMENTS ONLY MODE ACTIVE (Default: Depletions calculated for FDR/Bonf penalty, but hidden) **\n");
    if (concise_mode) printf("** CONCISE MODE: Displaying ONLY Significant terms with Precedence <= %d **\n", concise_limit);
    printf("==========================================================================================================================\n");
    
    // Terminal progress log
    printf("-> Processing statistics...\n");
    printf("-> Assembling and routing tables...\n");

    FILE *out = fopen(results_fname, "w");
    if (out) {
        if (concise_mode) {
            fprintf(out, "Namespace\tk_val\tK_val\tExpect\tP_Value\t%s\tEnr/Dep\tPrec\tClusterID\tSig_Level\tIs_Slim\tGOTerm\tDescription\n", adj_p_label);
        } else {
            fprintf(out, "Namespace\tk_val\tK_val\tExpect\tP_Value\t%s\tEnr/Dep\tPrec\tClusterID\tReg_Prec\tReg_ClstID\tSig_Level\tIs_Slim\tGOTerm\tDescription\n", adj_p_label);
        }

        for (int i = 0; i < num_tests; i++) {
            if (!include_depletions && strcmp(results[i].res_type, "Dep") == 0) continue;
            
            if (results[i].p_value <= 0.05) {
                if (concise_mode) {
                    if (!results[i].is_significant || results[i].precedence[0] > concise_limit) continue;
                }

                char sig_str[64] = "Not Sig";
                char prec_str[16] = "-", reg_prec_str[16] = "-";
                
                if (results[i].is_significant) {
                    if (use_bonferroni) strcpy(sig_str, "Sig (Bonf)");
                    else strcpy(sig_str, "Sig (FDR)");
                    
                    if (results[i].has_equal_pval_ancestor[0]) sprintf(prec_str, "%d*", results[i].precedence[0]);
                    else sprintf(prec_str, "%d", results[i].precedence[0]);
                    
                    if (results[i].has_equal_pval_ancestor[1]) sprintf(reg_prec_str, "%d*", results[i].precedence[1]);
                    else sprintf(reg_prec_str, "%d", results[i].precedence[1]);
                } else {
                    strcpy(sig_str, "Sig (Raw)");
                }

                const char *ns_short = get_short_ns(results[i].term->namespace);
                const char *is_slim_str = results[i].term->is_slim ? "Yes" : "No";

                if (concise_mode) {
                    fprintf(out, "%s\t%d\t%d\t%.2f\t%.4e\t%.4e\t%s\t%s\t%d\t%s\t%s\t%s\t%s\n",
                            ns_short, results[i].k_count, results[i].K_count, 
                            results[i].expected_val, results[i].p_value, results[i].adj_p_value, results[i].res_type, prec_str, results[i].cluster_id[0], sig_str, is_slim_str, results[i].term->id, results[i].term->name);
                } else {
                    fprintf(out, "%s\t%d\t%d\t%.2f\t%.4e\t%.4e\t%s\t%s\t%d\t%s\t%d\t%s\t%s\t%s\t%s\n",
                            ns_short, results[i].k_count, results[i].K_count, 
                            results[i].expected_val, results[i].p_value, results[i].adj_p_value, results[i].res_type, prec_str, results[i].cluster_id[0], reg_prec_str, results[i].cluster_id[1], sig_str, is_slim_str, results[i].term->id, results[i].term->name);
                }
            }
        }
        fclose(out);
        printf("=> Full statistical results saved to '%s'\n", results_fname);
    }
        
    // Write Precedence Graphs and Edges (Structural Only) to Files
    if (!concise_mode) {
        FILE *graph_out = fopen(graphs_fname, "w");
        if (graph_out) {
            fprintf(graph_out, "==============================================================================================================\n");
            fprintf(graph_out, "GO TERM PRECEDENCE GRAPHS\n");
            fprintf(graph_out, "==============================================================================================================\n\n");

            const char *ns_labels[] = {"BIOLOGICAL_PROCESS", "CELLULAR_COMPONENT", "MOLECULAR_FUNCTION"};

            for (int t = 0; t < num_res_types; t++) {
                for (int n = 0; n < 3; n++) {
                    bool has_graphs = false;
                    for (int cid = 1; cid < global_cluster_id; cid++) {
                        for (int i = 0; i < num_tests; i++) {
                            if (results[i].is_significant && strcmp(results[i].res_type, res_types[t]) == 0 && strcmp(results[i].term->namespace, ns_types[n]) == 0 && results[i].cluster_id[0] == cid) {
                                
                                bool is_cluster_root = true;
                                for (int p = 0; p < results[i].term->num_parents; p++) {
                                    GoTerm *parent_term = results[i].term->parents[p].target_ptr;
                                    EnrichmentResult *p_res = (strcmp(res_types[t], "Enr") == 0) ? (EnrichmentResult*)parent_term->enrich_res : (EnrichmentResult*)parent_term->deplet_res;
                                    if (p_res && p_res->cluster_id[0] == cid) {
                                        is_cluster_root = false;
                                        break;
                                    }
                                }
                                
                                if (is_cluster_root) {
                                    if (!has_graphs) {
                                        fprintf(graph_out, "--- %s: %s Graphs ---\n\n", ns_labels[n], res_types[t]);
                                        has_graphs = true;
                                    }
                                    int offsets[100] = {0};
                                    print_cluster_graph(&results[i], cid, 0, offsets, true, graph_out, res_types[t], 0); 
                                    fprintf(graph_out, "\n\n");
                                }
                            }
                        }
                    }
                }
            }
            fclose(graph_out);
            printf("=> Text representation of Precedence Graphs saved to '%s'\n", graphs_fname);
        }

        FILE *edge_out = fopen(edges_fname, "w");
        if (edge_out) {
            fprintf(edge_out, "Parent\tChild\tClusterID\tNamespace\tEnr/Dep\tEdgeStyle\n"); 
            
            GoTerm **bfs_queue = malloc(MAX_TERMS * sizeof(GoTerm*));
            
            for (int i = 0; i < num_tests; i++) {
                if (!include_depletions && strcmp(results[i].res_type, "Dep") == 0) continue;
                
                if (results[i].is_significant && results[i].cluster_id[0] > 0) {
                    
                    bool is_head_node = true;
                    
                    int n_c = results[i].term->num_children;
                    Relationship *c_rels = results[i].term->children;
                    
                    for (int c = 0; c < n_c; c++) {
                        GoTerm *child = c_rels[c].target_ptr;
                        EnrichmentResult *c_res = (strcmp(results[i].res_type, "Enr") == 0) ? (EnrichmentResult*)child->enrich_res : (EnrichmentResult*)child->deplet_res;
                        
                        if (c_res && c_res->is_significant && c_res->cluster_id[0] == results[i].cluster_id[0]) {
                            fprintf(edge_out, "%s\t%s\t%d\t%s\t%s\tsolid\n", 
                                results[i].term->id, child->id, results[i].cluster_id[0], results[i].term->namespace, results[i].res_type);
                        }
                    }
                    
                    int n_p = results[i].term->num_parents;
                    Relationship *p_rels = results[i].term->parents;
                    for (int p = 0; p < n_p; p++) {
                        GoTerm *parent = p_rels[p].target_ptr;
                        EnrichmentResult *p_res = (strcmp(results[i].res_type, "Enr") == 0) ? (EnrichmentResult*)parent->enrich_res : (EnrichmentResult*)parent->deplet_res;
                        if (p_res && p_res->is_significant && strcmp(parent->namespace, results[i].term->namespace) == 0) {
                            is_head_node = false;
                            break;
                        }
                    }
                    
                    if (is_head_node) {
                        int q_head = 0, q_tail = 0;
                        global_session_id++; 
                        
                        for (int p = 0; p < n_p; p++) {
                            GoTerm *parent = p_rels[p].target_ptr;
                            if (parent) {
                                parent->visited_marker = global_session_id;
                                bfs_queue[q_tail++] = parent;
                            }
                        }
                        
                        while (q_head < q_tail) {
                            GoTerm *curr = bfs_queue[q_head++];
                            EnrichmentResult *curr_res = (strcmp(results[i].res_type, "Enr") == 0) ? (EnrichmentResult*)curr->enrich_res : (EnrichmentResult*)curr->deplet_res;
                            
                            if (curr_res && curr_res->is_significant && strcmp(curr->namespace, results[i].term->namespace) == 0) {
                                fprintf(edge_out, "%s\t%s\t%d\t%s\t%s\tdotted\n", 
                                    curr->id, results[i].term->id, results[i].cluster_id[0], results[i].term->namespace, results[i].res_type);
                            } else {
                                int c_n_p = curr->num_parents;
                                Relationship *c_p_rels = curr->parents;
                                for (int p2 = 0; p2 < c_n_p; p2++) {
                                    GoTerm *gp = c_p_rels[p2].target_ptr;
                                    if (gp && gp->visited_marker != global_session_id) {
                                        gp->visited_marker = global_session_id;
                                        bfs_queue[q_tail++] = gp;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            free(bfs_queue);
            fclose(edge_out);
            printf("=> Raw DAG Edges saved to '%s'\n", edges_fname);
        }
    } else {
         printf("=> Concise mode enabled: Skipped generating Precedence Graphs and Edges files.\n");
    }

    // --- STEP 5: OUTPUT SIGNIFICANTLY ANNOTATED LIST CROSS-REFERENCE ---
    FILE *sf = fopen(sample_expanded_fname, "r");
    if (sf) {
        char sig_prot_fname[2048];
        if (prefix) snprintf(sig_prot_fname, sizeof(sig_prot_fname), "%s_Significantly_Annotated_List.tsv", prefix);
        else strcpy(sig_prot_fname, "Significantly_Annotated_List.tsv");
        
        FILE *spf = fopen(sig_prot_fname, "w");
        if (spf) {
            fprintf(spf, "Protein_ID\tGOTerm\tDescription\tEnr/Dep\tPrec\tP_Value\n");
            
            char line[1024];
            while (fgets(line, sizeof(line), sf)) {
                char prot[256], go[32];
                if (sscanf(line, "%255s %31s", prot, go) == 2) {
                    GoTerm *t = get_term(go);
                    if (t) {
                        EnrichmentResult *res = NULL;
                        if (t->enrich_res && ((EnrichmentResult*)t->enrich_res)->is_significant) res = (EnrichmentResult*)t->enrich_res;
                        else if (t->deplet_res && ((EnrichmentResult*)t->deplet_res)->is_significant && include_depletions) res = (EnrichmentResult*)t->deplet_res;
                        
                        if (res) {
                            if (concise_mode && res->precedence[0] > concise_limit) continue;
                            
                            char p_str[16];
                            if (res->has_equal_pval_ancestor[0]) sprintf(p_str, "%d*", res->precedence[0]);
                            else sprintf(p_str, "%d", res->precedence[0]);
                            
                            fprintf(spf, "%s\t%s\t%s\t%s\t%s\t%.4e\n", prot, t->id, t->name, res->res_type, p_str, res->p_value);
                        }
                    }
                }
            }
            fclose(spf);
            printf("=> Significantly annotated list saved to '%s'\n", sig_prot_fname);
        }
        fclose(sf);
    }
    
    printf("==========================================================================================================================\n\n");

    for (int i = 0; i < num_tests; i++) {
        for(int m=0; m<2; m++) {
            if (results[i].ancestors[m]) free(results[i].ancestors[m]);
            if (results[i].descendants[m]) free(results[i].descendants[m]);
        }
    }
    free(results);
}

// --- Main ---

int main(int argc, char **argv) {
    char *obo_filename = "go-basic.obo";
    bool custom_obo_provided = false;
    char *obo_data = NULL;
    char *query_id = NULL;
    char *bg_filename = NULL;
    char *sample_filename = NULL;
    char *prefix = NULL;
    char *slim_filename = NULL;
    double threshold = 0.05; 
    bool use_bonferroni = false;
    bool concise_mode = false;
    bool include_depletions = false;
    bool output_expanded = false;
    int concise_limit = 1;
    int opt;

    if (argc == 1) {
        print_help();
        return 0;
    }

    while ((opt = getopt(argc, argv, "o:q:b:s:d:Bhp:cm:De")) != -1) {
        if (opt == 'o') {
            obo_filename = optarg;
            custom_obo_provided = true;
        }
        else if (opt == 'q') query_id = optarg;
        else if (opt == 'b') bg_filename = optarg;
        else if (opt == 's') sample_filename = optarg;
        else if (opt == 'd') threshold = atof(optarg);
        else if (opt == 'B') use_bonferroni = true;
        else if (opt == 'p') prefix = optarg;
        else if (opt == 'm') slim_filename = optarg;
        else if (opt == 'D') include_depletions = true;
        else if (opt == 'e') output_expanded = true;
        else if (opt == 'c') {
            concise_mode = true;
            if (optind < argc && argv[optind][0] != '-') {
                concise_limit = atoi(argv[optind]);
                if (concise_limit < 1) concise_limit = 1;
                optind++;
            }
        }
        else if (opt == 'h') {
            print_help();
            return 0;
        }
    }

    if (access(obo_filename, F_OK) == -1) {
        if (!custom_obo_provided) {
            printf("\nError: Default OBO file 'go-basic.obo' not found in the current directory.\n");
            printf("Please supply an OBO ontology file using the -o option.\n");
            printf("Example: ./GOParsimony -o my_ontology.obo\n\n");
        } else {
            printf("\nError: The specified OBO file '%s' could not be found.\n\n", obo_filename);
        }
        return 1;
    }

    obo_data = read_file_to_buffer(obo_filename);
    if (!obo_data) {
        printf("Error: Failed to read OBO file '%s'.\n", obo_filename);
        return 1;
    }

    build_graph(obo_data);
    
    if (slim_filename) {
        char *slim_data = read_file_to_buffer(slim_filename);
        if (slim_data) {
            parse_goslim(slim_data);
            free(slim_data);
        } else {
            printf("Error: Could not read GO slims file '%s'.\n", slim_filename);
        }
    }

    link_children();
    
    if (query_id) {
        GoTerm *term = get_term(query_id);
        if (term) {
            printf("\n=======================================================================\n");
            printf(">>> Query Result for %s\n", term->id);
            printf("=======================================================================\n");
            printf("Name:       %s\n", term->name);
            printf("Namespace:  %s\n", term->namespace);
            
            printf("\n--- Direct Parents (%d) ---\n", term->num_parents);
            for (int i = 0; i < term->num_parents; i++) {
                printf("  [%s] -> %s (%s)\n", term->parents[i].type, term->parents[i].target_id, 
                       term->parents[i].target_ptr ? term->parents[i].target_ptr->name : "Unknown");
            }
            if (term->num_parents == 0) printf("  (None)\n");

            printf("\n--- Direct Children (%d) ---\n", term->num_children);
            for (int i = 0; i < term->num_children; i++) {
                printf("  [%s] -> %s (%s)\n", term->children[i].type, term->children[i].target_id, 
                       term->children[i].target_ptr ? term->children[i].target_ptr->name : "Unknown");
            }
            if (term->num_children == 0) printf("  (None)\n");

            if (term->num_reg_parents > 0) {
                printf("\n--- Regulatory Parents (%d) ---\n", term->num_reg_parents);
                for (int i = 0; i < term->num_reg_parents; i++) {
                    printf("  [%s] -> %s (%s)\n", term->reg_parents[i].type, term->reg_parents[i].target_id, 
                           term->reg_parents[i].target_ptr ? term->reg_parents[i].target_ptr->name : "Unknown");
                }
            }

            if (term->num_reg_children > 0) {
                printf("\n--- Regulatory Children (%d) ---\n", term->num_reg_children);
                for (int i = 0; i < term->num_reg_children; i++) {
                    printf("  [%s] -> %s (%s)\n", term->reg_children[i].type, term->reg_children[i].target_id, 
                           term->reg_children[i].target_ptr ? term->reg_children[i].target_ptr->name : "Unknown");
                }
            }
            
            global_session_id++;
            GoTerm **ancestors = NULL;
            int anc_count = 0, anc_cap = 0;
            collect_ancestors(term, &ancestors, &anc_count, &anc_cap, global_session_id, 0);
            qsort(ancestors, anc_count, sizeof(GoTerm*), cmp_goterm_ptrs);

            printf("\n--- Total Structural Ancestors (%d) ---\n", anc_count);
            for(int i = 0; i < anc_count; i++) {
                printf("  %s (%s)\n", ancestors[i]->id, ancestors[i]->name);
            }
            if(anc_count == 0) printf("  (None)\n");
            free(ancestors);

            global_session_id++;
            GoTerm **descendants = NULL;
            int desc_count = 0, desc_cap = 0;
            collect_descendants(term, &descendants, &desc_count, &desc_cap, global_session_id, 0);
            qsort(descendants, desc_count, sizeof(GoTerm*), cmp_goterm_ptrs);

            printf("\n--- Total Structural Descendants (%d) ---\n", desc_count);
            for(int i = 0; i < desc_count; i++) {
                printf("  %s (%s)\n", descendants[i]->id, descendants[i]->name);
            }
            if(desc_count == 0) printf("  (None)\n");
            free(descendants);

            global_session_id++;
            GoTerm **reg_ancestors = NULL;
            int reg_anc_count = 0, reg_anc_cap = 0;
            collect_ancestors(term, &reg_ancestors, &reg_anc_count, &reg_anc_cap, global_session_id, 1);
            qsort(reg_ancestors, reg_anc_count, sizeof(GoTerm*), cmp_goterm_ptrs);

            printf("\n--- Total Regulatory Ancestors (%d) ---\n", reg_anc_count);
            for(int i = 0; i < reg_anc_count; i++) {
                printf("  %s (%s)\n", reg_ancestors[i]->id, reg_ancestors[i]->name);
            }
            if(reg_anc_count == 0) printf("  (None)\n");
            free(reg_ancestors);

            global_session_id++;
            GoTerm **reg_descendants = NULL;
            int reg_desc_count = 0, reg_desc_cap = 0;
            collect_descendants(term, &reg_descendants, &reg_desc_count, &reg_desc_cap, global_session_id, 1);
            qsort(reg_descendants, reg_desc_count, sizeof(GoTerm*), cmp_goterm_ptrs);

            printf("\n--- Total Regulatory Descendants (%d) ---\n", reg_desc_count);
            for(int i = 0; i < reg_desc_count; i++) {
                printf("  %s (%s)\n", reg_descendants[i]->id, reg_descendants[i]->name);
            }
            if(reg_desc_count == 0) printf("  (None)\n");
            free(reg_descendants);

            printf("=======================================================================\n\n");
        }
        else {
            printf("\nError: Query term '%s' was not found in the OBO file.\n\n", query_id);
        }
    }

    int bg_N = 0;
    int sample_n = 0;

    char bg_expanded_fname[2048] = {0};
    char sample_expanded_fname[2048] = {0};

    Annotation *bg_annots = NULL;
    int bg_annot_count = 0;

    if (bg_filename) {
        bg_annots = load_background_annotations(bg_filename, &bg_annot_count);
        if (bg_annots) {
            bg_N = expand_and_write(bg_annots, bg_annot_count, bg_filename, "background", true, prefix, bg_expanded_fname);
        }
    }

    Annotation *sample_annots = NULL;
    int sample_annot_count = 0;

    if (sample_filename) {
        if (bg_annots) {
            sample_annots = load_sample_annotations(sample_filename, bg_annots, bg_annot_count, &sample_annot_count);
            if (sample_annots) {
                sample_n = expand_and_write(sample_annots, sample_annot_count, sample_filename, "sample", false, prefix, sample_expanded_fname);
            }
        } else {
            printf("\nError: Background annotations must be provided to resolve sample IDs.\n");
        }
    }

    // Safely free the intermediate annotation structures as their components (bg_count and sample_count 
    // mapped to the global term graph) have already been securely logged via expand_and_write
    if (bg_annots) free(bg_annots);
    if (sample_annots) free(sample_annots);

    if (bg_N > 0 && sample_n > 0) {
        calculate_enrichment(bg_N, sample_n, threshold, use_bonferroni, prefix, concise_mode, concise_limit, include_depletions, sample_expanded_fname);
    }

    if (!output_expanded) {
        if (bg_filename) remove(bg_expanded_fname);
        if (sample_filename) remove(sample_expanded_fname);
    }

    free(obo_data);
    cleanup_graph();
    return 0;
}

