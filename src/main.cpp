/*
 *	Author: Harry
 *	Email: khchanak@cse.ust.hk
 */

#include "yoo_alg.h"
#include "irtree.h"
#include "frac.h"
#include "data_utility.h"
#include <iostream>
#include <fstream>
#include <random>
#include "SmallestEnclosingCircle.h"
#include <cstddef>
#include <cmath>
#include "stringgenc.h"

using namespace std;
using std::size_t;
using std::vector;

IRTree_t IRTree_v;
bst_t* IF_v;

colocation_stat_t stat_v;

bool debug_mode = false;

//----------------------------------
//user parameter
float dist_thr;
float min_sup;
int cost_tag; //1= participation based; 2=fraction based
double fea_highest_freq;

//------------------------------------
void colocation_patterns_support(int cost);

void colocation();

void build_IF(data_t* data_v);

void gen_syn();

double get_highest_freq(bst_node_t* bst_node_v);

fsi_set_t* read_patterns();

int main(int argc, char* argv[])
{

  //     batch_gen_syn_data2();
    colocation();

//    colocation_patterns_support(3);
    return 0;
}


void colocation()
{
    int i;
    colocation_config_t* cfg;
    FILE *r_fp, *r_fp2;
    data_t* data_v;
    fsi_set_t** result;

    memset(&stat_v, 0, sizeof(colocation_stat_t));

    //Read the cofig.
    printf("Reading configuration ...\n");
    cfg = read_config_colocation();
    cfg->dim = 2;

    cost_tag = cfg->cost;
    dist_thr = cfg->dist_thr;
    min_sup = cfg->min_sup;
  //  min_conf = cfg->min_conf;

    //Read the data.
    printf("Reading data ...\n");
    data_v = read_data_colocation(cfg);

#ifndef WIN32
    float sys_t, usr_t, usr_t_sum = 0;
    struct rusage IR_tree_sta, IR_tree_end;

    GetCurTime(&IR_tree_sta);
#endif

    //Option 1: Build the tree from scratch.
    //Build the IR-tree.
    if (cfg->tree_tag == 0) {
        printf("Building IR-tree ...\n");
        build_IRTree(data_v);

        print_and_check_tree(1, cfg->tree_file);
        //check_IF( );
    } else {
        //Option 2: Read an existing tree.
        printf("Reading IR-Tree ...\n");
        read_tree(cfg->tree_file);
    }

#ifndef WIN32
    GetCurTime(&IR_tree_end);
    GetTime(&IR_tree_sta, &IR_tree_end, &stat_v.irtree_build_time, &sys_t);
#endif

    //Build Inverted file
    build_IF(data_v);

    //---
    if (cost_tag == 1) {
        printf("participation based: min_sup:%lf\n", min_sup);
    } else {
        fea_highest_freq = get_highest_freq(IF_v->root);
        printf("Fraction based: Highest freq:%lf\t min_sup:%lf\n", fea_highest_freq, min_sup);
    }
    //---

    //    //Get the whole range.
    //    MBR = get_MBR_node( IRTree_v.root, IRTree_v.dim);

    printf("Running apriori (Method %d):\n", cfg->alg_opt);

#ifndef WIN32
	
	float pre_t = 0;
	struct rusage pre_sta, pre_end;
	GetCurTime(&pre_sta);
#endif
	
	//Pre-computation of N_o_f
	precomputation(data_v, cfg->dist_thr);
	
#ifndef WIN32
	GetCurTime(&pre_end);
	GetTime(&pre_sta, &pre_end, &pre_t, &sys_t);
	printf("pre\ttime:%0.5lf\n", pre_t);
#endif
	
    for (i = 0; i < cfg->query_n; i++) //i==1
    {

#ifndef WIN32
        struct rusage query_sta, query_end;
        GetCurTime(&query_sta);

#endif

        printf("Query #%i ...\n", i + 1);

        if (cfg->alg_opt == 5)
            result = joinless_mining(data_v, cfg->obj_n, cfg->key_n);
        else
            result = apriori(cfg->alg_opt, cfg->obj_n, cfg->key_n);

#ifndef WIN32
        GetCurTime(&query_end);

        GetTime(&query_sta, &query_end, &usr_t, &sys_t);
        usr_t_sum += usr_t;

		printf("Time: %f\n\n", usr_t + pre_t);
		//pre-computation time should also be included in the query time
		stat_v.q_time = usr_t_sum / (i + 1) + pre_t;

#endif

        //Print the query results.
        if (i == 0) {
            if ((r_fp = fopen(COLOCATION_RESULT_FILE, "w")) == NULL) {
                fprintf(stderr, "Cannot open the colocation_result file.\n");
                exit(0);
            }

            //Print the query result.
            print_fsi_set(result, cfg->key_n, r_fp);

//            if ((r_fp2 = fopen(COLOCATION_RESULT_FILE2, "w")) == NULL) {
//                fprintf(stderr, "Cannot open the colocation_result2 file.\n");
//                exit(0);
//            }
			
           // print_maximal_fsi_set(result, cfg->key_n, r_fp2);

			fclose(r_fp);
			//fclose(r_fp2);

        }
        //release the result memory
        for (int k = 0; k < cfg->key_n; k++)
            if (result[k] != NULL)
                release_fsi_set(result[k]);

        print_colocation_stat(cfg, i + 1);

        //reset some statistcs for the next query
        stat_v.S1_sum = 0.0;
        stat_v.S2_sum = 0.0;
        stat_v.S3_sum = 0.0;
        stat_v.S5_sum = 0.0;

        stat_v.S1_time = 0.0;
        stat_v.S2_time = 0.0;
        stat_v.S3_time = 0.0;
        stat_v.S4_time = 0.0;
        stat_v.S5_time = 0.0;
    }

    free(cfg);
}

/*
 * Construct the inverted file @IF_v based on the data
 */
void build_IF(data_t* data_v)
{

    bst_node_t* bst_node_v;

    IF_v = bst_ini();

    //Insert all the objects to construct the IF
    for (int i = 0; i < data_v->obj_n; i++) {
        //        if(i%100==0)
        {

            bst_node_v = bst_search(IF_v, data_v->obj_v[i].fea);

            if (bst_node_v != NULL) {
                add_obj_set_entry(&data_v->obj_v[i], bst_node_v->p_list_obj);
            } else //bst_node_v = NULL
            {
                //Insert a new keyword in the binary tree.
                bst_node_v = (bst_node_t*)malloc(sizeof(bst_node_t));
                memset(bst_node_v, 0, sizeof(bst_node_t));

                /*s*/
                stat_v.memory_v += sizeof(bst_node_t);
                if (stat_v.memory_v > stat_v.memory_max)
                    stat_v.memory_max = stat_v.memory_v;
                /*s*/

                //Update the posting list.
                bst_node_v->key = data_v->obj_v[i].fea;
                bst_node_v->p_list_obj = alloc_obj_set();

                add_obj_set_entry(&data_v->obj_v[i], bst_node_v->p_list_obj);
                bst_insert(IF_v, bst_node_v);
            }
        }
    }
}

double get_highest_freq(bst_node_t* bst_node_v)
{

    if (bst_node_v == NULL)
        return 0;

    //    double freq = bst_node_v->p_list_obj->obj_n;
    //    double l = get_highest_freq(bst_node_v->left);
    //    double r = get_highest_freq(bst_node_v->right);
    //
    //    return fmax(freq,fmax(l,r));

    return fmax(bst_node_v->p_list_obj->obj_n,
        fmax(get_highest_freq(bst_node_v->left), get_highest_freq(bst_node_v->right)));
}

fsi_set_t* read_patterns()
{
    fsi_set_t* result;
    fsi_t *fsi_v, *fsi_cur;
    ifstream readConfig;

    string STRING;
    vector<string> STRING2, arg;

    result = alloc_fsi_set();

    readConfig.open("pattern.txt");

    //fsi_cur always at the end of the list
    fsi_cur = result->head;
    int i = 0;
    while (!readConfig.eof()) {
        getline(readConfig, STRING);
        //		cout << "STRING: "<< STRING << endl;
        if (STRING.size() == 0)
            break;

        //putting config into vars.
        STRING2 = split(STRING, "  ");

        int size = (int)STRING2.size() - 1;
        fsi_v = alloc_fsi(size);

        for (int j = 0; j < STRING2.size() - 1; j++) {
            //			cout << "STRING2  " << j << ": "<< STRING2[j] << endl;
            fsi_v->feaset[j] = atoi(STRING2[j].c_str());
        }

        //		print_fsi(fsi_v, stdout);

        fsi_cur->next = fsi_v;
        fsi_cur = fsi_cur->next;
        i++;
    }

    return result;
}


/*
 * Additional method for computing support of given fsi_sets
 * It reads "pattern.txt" line by line
 * Each line corrsponds to one fsi_set, the 0 to n-1 numbers are the features
 * Ouput result in console
 * Note that d and min-sup need to be set in config.txt
 //	1: participation
 //	2: fraction-score
 //	3: partitioning
 //	4: construction
 //	5: enumeration
 // 6: gen_scalability_data2 (put in here since it need to read data first)
 */
void colocation_patterns_support(int cost)
{
	colocation_config_t* cfg;
	data_t* data_v;
	
	memset(&stat_v, 0, sizeof(colocation_stat_t));
	
	//Read the cofig.
	printf("Reading configuration ...\n");
	cfg = read_config_colocation();
	cfg->dim = 2;
	
	//cost_tag = cfg->cost;
	cost_tag = cost;
	dist_thr = cfg->dist_thr;
	min_sup = cfg->min_sup;
	//min_conf = cfg->min_conf;
	
	//Read the data.
	printf("Reading data ...\n");
	data_v = read_data_colocation(cfg);
	
	//Option 1: Build the tree from scratch.
	//Build the IR-tree.
	if (cfg->tree_tag == 0) {
		printf("Building IR-tree ...\n");
		build_IRTree(data_v);
		
		print_and_check_tree(1, cfg->tree_file);
		//check_IF( );
	} else {
		//Option 2: Read an existing tree.
		printf("Reading IR-Tree ...\n");
		read_tree(cfg->tree_file);
	}
	
	
	//Build Inverted file
	build_IF(data_v);
	fea_highest_freq = get_highest_freq(IF_v->root);
	printf("Highest freq:%lf\n", fea_highest_freq);
	
	//---
	if (cost_tag == 1)
		printf("participation based:\n");
	else if(cost_tag==2){
		printf("Fraction based:\n");
		//Pre-computation of N_o_f
		precomputation(data_v, cfg->dist_thr);
		
		
	}else if(cost_tag==3)
		printf("partitioning based: \n");
	else if(cost_tag==4)
		printf("construction based: \n");
	else if(cost_tag==5)
		printf("enumeration based: \n");
	
	//---
	if(cost==6)
		gen_scalability_data2(data_v);
	//---
	
	fsi_set_t* patterns = read_patterns();
	
	fsi_t* fsi_cur = patterns->head->next;
	while (fsi_cur != NULL) {
		comp_support(cfg->alg_opt, fsi_cur, cfg->obj_n);
		print_fsi(fsi_cur, stdout);
		fsi_cur = fsi_cur->next;
	}
	
	release_fsi_set(patterns);
	free(cfg);
}
