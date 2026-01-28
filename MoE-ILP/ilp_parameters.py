from enum import Enum

class OBJ(Enum):
    LAMBDA = 1
    FEASABILITY = 2
    HRELATION = 3
    LATENCY = 4
    WORKLOAD_BALANCE = 5


class instanceParameter:
    
    def __init__(self):
        
        self.num_partitions = 6
        self.objective = OBJ.LAMBDA
        
        self.memory_constraint = True
        self.memory_bound = 27

        self.time_limit = 24 #hours

        # self.min_freq_pairs = 2000
        self.min_freq_trans = 2000
        self.per_layer_top_q = 0.98
        self.per_trans_top_q = 0.90
        self.weight_trans = 0.2

        self.layer_threshold = 10
        self.cost_mla = 1

        self.eps =1e-4
        self.load_threshold  = 0.3

        #MLAs
        self.lb_mla = 1
        self.ub_mla = 3
        self.max_mlas_in_partition = 3
        self.next_mla = True

        #initial solution path
        self.initial_solution_file = ""

    