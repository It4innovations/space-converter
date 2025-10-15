#include <stdio.h>
#include <stdlib.h>

#ifndef __NUMA_H__
#define __NUMA_H__


/*
  General infromation on the machine where the simulation is running.
  So far we only need processes/node and the number of processes;
*/
struct NumaArchitecture{
  int processes_per_node;
  int processes;
  int nodes;
};

/*
  information of this process:
  our process index (mpi rank) and node index (mpi rank % processes per node)
*/
struct NumaInfo{
  int process_idx;
  int node_idx;
};

/*
  static datastructes: each compile units will have its own data,
  but this is ok as long as they are read-only.
 */
static NumaArchitecture *numa_architecture = NULL;
static NumaInfo *numa_info = NULL;

/*
  grab numa info from global variables and env. 
  Indeed ATM you need to specify GADGET_PROCESSS_PER_NODE, as it is too complicated to find it out pragmatically.
*/
static void numa_init(){
  if (numa_architecture == NULL){
    //
    numa_architecture = (NumaArchitecture*) calloc(1, sizeof(NumaArchitecture));
    char *GADGET_PROCESSS_PER_NODE = getenv("GADGET_PROCESSS_PER_NODE");
    PANIC_IF(GADGET_PROCESSS_PER_NODE==NULL,"Environmnet variable GADGET_PROCESSS_PER_NODE should be specified to use NUMA");
    numa_architecture->processes_per_node = atoi(GADGET_PROCESSS_PER_NODE);
    numa_architecture->processes = NTask;
    PANIC_IF( (numa_architecture->processes%numa_architecture->processes_per_node )!=0,
	      "Numa Library so far will only work with if all nodes have the same number of processes. "
	      "procs = %d; procs_per_node = %d",
	      numa_architecture->processes, numa_architecture->processes_per_node);
    numa_architecture->nodes = numa_architecture->processes/numa_architecture->processes_per_node;
  }

  if (numa_info == NULL){
    numa_info = (NumaInfo*) calloc(1, sizeof(NumaInfo));
    numa_info->process_idx = ThisTask;
    numa_info->node_idx = numa_info->process_idx % numa_architecture->processes_per_node;
  }
}


static struct NumaArchitecture* get_numa_architecture(){
  numa_init();
  return numa_architecture; 
}

static struct NumaInfo get_numa_of_process_info(int process_idx){
  numa_init();
  struct NumaInfo numa_info;
  numa_info.process_idx = ThisTask;
  numa_info.node_idx = numa_info.process_idx % numa_architecture->processes_per_node;
  return  numa_info;
}

static NumaInfo* get_numa_info(){
  numa_init();
  return numa_info;
}

//don't forget to free resources if you used numa data.
static void numa_finalize(){
  free(numa_info);
  free(numa_architecture);
}

#endif
