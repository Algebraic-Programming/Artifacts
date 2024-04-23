/*
Copyright 2024 Huawei Technologies Co., Ltd.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

@author Aikaterini Karanasiou
*/

/*---------------------------------------------------------------
 |    Heuristic NO2		                                 |		
 *---------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <bits/stdc++.h>
#include <iostream>
#include <fstream>
#include<string>
#include <iostream>
#include <fstream>
#define MAXLINE       1000   /* max length of input line */

using namespace std;

struct topological_var {
	int number;
	list<int> t_ordering;
} ;

class hyperDag {
		
public:
	int* cluster_sources(topological_var var2, int *m_cl);
	int *sorting (topological_var var1, int *m_id, int *m_value) ;
	topological_var topological (topological_var my_variables, int *incoming);
	void scheduling (fstream& out_file_two);
	void transform_hyperdag_to_dag ();
	void readHyperDag(ifstream& infile, fstream& out_file);
	
private:
	int num_of_hyperedges, num_of_vertices, num_of_pins; /* number of hyperedges, vertices, lines*/
	int DAG_edges;// number of edges of the constructed DAG
	
	int *input_file; // stores input arcs
	int *workload_weight; // stores the weight of every node
	int *communication_delay_cost; // stores the communication dealy
	
	int n_pr1;
	int n_pr2;
	int n_cost;
	int g, l;
	int pr; //number of processors
	int cluster;
	int superstep; // number of total supersteps
	
	int *matrix_processor; // in which processor is executed every vertex
	int *matrix_superstep; // in which superstep is executed every vertex
	int *cost_of_processor; // counts the workload of every processor
	int *incoming_edges ; //number of incoming edges of every vertex
	int *outgoing_edges ; //number of outgoing edges of every vertex
	list <int> *destination_vertices; //the outgoing vertices of every vertex
	list <int> *incoming_vertices_destinations;// the incoming edges of every vertex 

};

/*cluster the input sources*/

int*  hyperDag::cluster_sources(topological_var var2, int *m_cl){
	list<int>::iterator it, it2;
	bool flag ;
	cluster = -1;
	for (it = var2.t_ordering.begin(); it !=  var2.t_ordering.end(); ++it){
		if(m_cl[*it] <0){
			flag = false;
			//cout<<"we check source:	"<<*it<<endl;
			for (auto element : destination_vertices[*it]) {
				//cout<<"destination vertices from source	:"<<*it<<endl;
				for (auto element2 : incoming_vertices_destinations[element]) {
					//cout<<"check vertex\t"<<element2<<"\tthat is incoming from \t"<<element<<"\t that is destination from \t"<<*it<<endl;
					for (it2 = var2.t_ordering.begin(); it2 != var2.t_ordering.end(); ++it2){
						if (element2 == *it2  && *it2 != *it) {
							flag = true;
							//cout<<endl<<"Vertex:	 "<< element <<"	has incoming edges from source:	"<<*it<<" and source:	"<< *it2<<"--> assign them to the same cluster:	"<<endl;
							
							//check if they are already assigned to a cluster
							if (m_cl[*it]== -1 && m_cl[*it2] == -1) {
								cluster++;
								//cout<<"they are unasigned"<<endl;
								//cout<<"put	"<<*it<<"\t and \t"<<*it2<<"\t in cluster\t"<<cluster<<endl;
								m_cl[*it] = m_cl[*it2]= cluster;
								
							}
							if (m_cl[*it] != -1) {
								m_cl[*it2] = m_cl[*it] ;
								//cout<<"we assign\t"<<*it2<<"\t to the cluster\t"<<m_cl[*it]<<"\t that belongs \t"<< *it<<endl;
							} 
							
							if (m_cl[*it2] != -1) {
								m_cl[*it] = m_cl[*it2] ;
								//cout<<"we assign\t"<<*it<<"\t to the cluster\t"<<m_cl[*it2]<<"\t that belongs \t"<< *it2<<endl;
								
							}
							if (m_cl[*it] == -2 || m_cl[*it2] == -2) {
								cout<< "error in clustering"<<endl;
								exit(1);
							}
							

						}
						
					}	
		
				}
			}
			if (!flag) {
				cluster++;
				m_cl[*it] = cluster;
				//cout<<"assign :  "<<m_cl[*it]<<"\t to cluster:\t"<<cluster<<endl; 
				
			}
		
		}
	
	}
	//cout<<endl<<"we have:	"<<cluster<< "	clusters and	"<< var2.number<<"	vertices to order"<<endl;	
	return m_cl;

}


int * hyperDag::sorting (topological_var var1, int *v_id, int *v_weight) {
		list<int>::iterator it;
		int k =0;
		for (it = var1.t_ordering.begin(); it != var1.t_ordering.end(); ++it) {
			v_weight[k] = workload_weight[*it];
			//cout<<"assign this value workload["<<*it<<"] ="	<< workload_weight[*it]<<endl;
			v_id[k] = *it;
			k++;
		}
		
		if ((k)!=var1.number) {
			cout<<"ERROR_SORTING_FUNC"<<endl;
			cout<<"k=	"<<k-1<<"var1.number=	"<<var1.number<<endl;
			exit(1);
		}
		
		/*sort according to v_weight*/
		int key;
		int key2;
		for (int i = 1; i < var1.number; i++) {
			key = v_weight[i];
			key2 = v_id[i];
			int j = i - 1;

			/* Move elements of arr[0..i-1], that are
			greater than key, to one position ahead
			of their current position */
			while (j >= 0 && v_weight[j] > key) {
				v_weight[j + 1] = v_weight[j];
				v_id[j+1] = v_id[j];
				j = j - 1;
			}
			v_weight[j + 1] = key;
			v_id[j+1] = key2;
		}
	
		for (int i = 0; i < var1.number; i++) {
			//cout<<"vertex_weight["<<i<<"]=	"<<v_weight[i]<<"	superstep:"<<superstep<<endl;
		}
		return v_id;
}



/*topological sort of an input matrix*/ 
topological_var hyperDag:: topological (topological_var my_variables, int *incoming){
	my_variables.number = 0;
	for (int i =0;i<num_of_vertices;i++) {
		if (incoming[i] == 0 && matrix_processor[i]<0) {
			//if there is no incoming edge put into the list an this vertex has not been assesed to any processor;
			my_variables.t_ordering.push_back(i);
			//cout<<"take vertex:	"<<i<<endl;
			my_variables.number++;
		}
	}
	return my_variables;
}

/*HyperDag Scheduling*/

void hyperDag::scheduling (fstream& out_file_two) {

	int  j;
	list<int>::iterator it;
	list <int> commun_cost_of_superstep; // counts the total communication cost of every superstep
	list <int> maximum_workload_cost; // counts the total workload cost of every superstep 
	list<int> *processors = new list <int> [pr+1]; // list that include the vertices from each processor
	int *workload_processor = new int [pr]; //number of workload for each processor in each superstep
	int *vertex_weight ;
	int *vertex_id ;
	int sum = 0;
	superstep = -1;
	
	int *incoming_2 = new int [num_of_vertices];
	
	for (int i=0;i<num_of_vertices;i++) incoming_2[i] = incoming_edges[i];
	for (int i=0;i<pr;i++) workload_processor[i] = 0;

	
	/*put vertices with no incoming edges to the processors --schedule*/
	/*while: not all of the vertices of the DAG have not been scheduled*/
	
	while (sum < num_of_vertices) {	
		superstep++;
		topological_var my_var;
		my_var = topological (my_var, incoming_2);
		if (my_var.number == 0) cout<< "error:no sources any more"<< endl;
		
		if (superstep == 0) {
			//cout<<"superstep == 0 --> clustering"<<endl;
			/*cluster the input sources*/
			int *matrix_cluster = new int[num_of_vertices];
			for (int i =0;i<num_of_vertices;i++) matrix_cluster[i] = -2;
		
			int k =0;
			for (it = my_var.t_ordering.begin(); it != my_var.t_ordering.end(); ++it){
				matrix_cluster[*it] = -1;
				k++;
			}
			
			matrix_cluster = cluster_sources(my_var, matrix_cluster);
			
			//test
			if (cluster == 0) {
				cout<<"error: cluster == 0"<<endl;
				exit(1);
			}
			//create the clusters
			list <int> *final_clusters = new list<int> [cluster + 1];
			for (int i=0;i<num_of_vertices;i++){
				// cout<<"check vertex:\t"<<i<< "\t "<<endl;
				if(matrix_cluster[i]>=0 ) {
					final_clusters[matrix_cluster[i]].push_back(i);
					//cout<<"put vertex\t"<<i<<"\t in cluster \t"<<matrix_cluster[i]<<endl;
				}
			}
			
			//schedule according to the clusters
			int pr_id = pr-1;
			for (int i=0;i<=cluster;i++){
				if(pr_id<0 ) pr_id = pr-1;
				for (auto element:final_clusters[i]){
					
					matrix_processor[element] = pr_id ; 
					matrix_superstep[element] = superstep ; 
					
					processors[pr_id].push_back(element);
					workload_processor[pr_id]++;
					
					sum++;
					
					//cout<<"\t schedule the vertex\t"<<element<<"\t in the processor\t"<<pr_id<<endl;
					//traverse their outgoing and update the incoming
					for (auto element2 : destination_vertices[element]) {	
						incoming_2[element2]--;
					}
				}
				pr_id--;
			}
			
			
			
			
			for (int i=0;i<pr;i++){
				//cout<<"processor\t"<<i<<endl;
				for (auto element:processors[i]){
					//cout<<"\t we have in processor:\t"<<i<<"\t vertex:\t"<<element<<endl;
				}
			}
			
		}
		
		else {
			/*sort the sources and schedule according to weights*/
			vertex_weight = new int [my_var.number];
			vertex_id = new int [my_var.number];
			for (int i = 0; i < my_var.number; i++) vertex_id[i] = vertex_weight[i] = 0;
			
			vertex_id = sorting (my_var, vertex_id, vertex_weight) ;
			 
			for (int i=0;i<pr;i++) workload_processor[i] = 0;
			int processor_id = pr-1;
		
			for(int i=0; i<my_var.number;i++){
				

				if (processor_id<0) processor_id = pr-1;
				matrix_processor[vertex_id[i]] = processor_id; // vertex *it is assesed to the processor: processor_id
				matrix_superstep[vertex_id[i]] = superstep ; // vertex *it is assesed to superstep: superstep
				
				processors[processor_id].push_back(vertex_id[i]);
				workload_processor[processor_id]++; // how many vertices has "processor_id" in this superstep
				
				sum++;
				processor_id--;	
			}
			
			//update the incoming and outgoing edges of each vertex 
			//count the cost of communication delay of each superstep by adding the outgoing edges of each vertex

			for (it = my_var.t_ordering.begin(); it != my_var.t_ordering.end(); ++it){
				//traverse their outgoing
				for (auto element : destination_vertices[*it]) {	
					incoming_2[element]--;
				}		
			}

		}
		
		
			
		
		
		//printf("\n we have until now %d scheduled of %d\n", sum , num_of_vertices);
		/*take the next ones if they are not dependent*/
		/*first traverse to every processor*/
		/*Attention to dependencies*/
		/* we want to schedule "element"*/
		list<int>::iterator it3;
		for (it = my_var.t_ordering.begin(); it != my_var.t_ordering.end(); ++it){
			for (auto element : destination_vertices[*it]) {	
				if(matrix_processor[element]== -1) {		
				 	bool flag = true;
				//	printf("\n we have (%d, %d), we check if we can take %d in the processor %d \n",*it, element, element, matrix_processor[*it]);
					/*check the incoming edges of the vertex that we want to schedule*/
					for (it3 = incoming_vertices_destinations[element].begin(); it3!= incoming_vertices_destinations[element].end();++it3) {	
						//printf("\n We want to take %d in processor %d. %d has incoming edge from %d that is in processor %d \n", element, matrix_processor[*it], element, *it3, matrix_processor[*it3] );				
						if ( matrix_processor[*it3] != matrix_processor[*it] || matrix_processor[*it3] == -1 ){	
							//printf("\n we cannot take %d in the processor %d\n",element, matrix_processor[*it] );
							flag = false;
							if ( matrix_processor[*it] == -1) {
								cout<<"		vertex= "<< *it<<"	has matrix_processor["<<*it<<"] = "<< matrix_processor[*it]<<endl;
								printf("\n ERROR: it is unassigned\n");
								exit(1); 
							}	
						}	
					}
					
					/* schedule element*/
					if (flag == true){
						

						matrix_processor[element] = matrix_processor[*it];
						matrix_superstep[element] = superstep ;
						processors[matrix_processor[*it]].push_back(element);
						//printf("\n we take %d  in proccessor %d \n", element, matrix_processor[element] );
						sum++;
						//printf("\n we add %d to the processor %d\n ", element, matrix_processor[*it]);
						//update the matrix of incoming edges
						for (auto element2 : destination_vertices[element]) {	
							incoming_2[element2]--;
						}
					}	
				}
			}					

		}
		
		

	}
	

	/*Check the output: 1*/
	for (int i=0;i<pr; i++){
		for (auto element1 : processors[i]) {
			for (auto element2:destination_vertices[element1]) {
					//printf("\n we check vertex:%d (of processor %d and superstep:%d) that has outgoing to vertex:%d (pr=%d, super=%d)", element1,i,matrix_superstep[element1],element2,matrix_processor[element2],matrix_superstep[element2]);
				if (matrix_superstep[element1] == matrix_superstep[element2]) {
					if (matrix_processor[element1] != matrix_processor[element2]) {
						cout<<"ERROR: check the output 1"<<endl;
						printf("\n Vertex %d (scheduled in processor:%d, at superstep:%d) has outgoing edge to vertex:%d scheduled at processor:%d at superstep:%d \n",element1, matrix_processor[element1], matrix_superstep[element1], element2,  matrix_processor[element2], matrix_superstep[element2]);
					}
					
				}
			}
		}
	}
	
	/*Check the output: 2*/
	for (int i=0;i<num_of_vertices;i++) {
		for (auto element:destination_vertices[i]) {
			if (matrix_superstep[i]<matrix_superstep[element]) {
			}
			if (matrix_superstep[i] == matrix_superstep[element]) {
				if(matrix_processor[i]!= matrix_processor[element]){
					cout<<"ERROR:check the output 2"<<endl;
					cout<< "they belong to the same superstep, they should be in the same processor\t("<<i<<"\t"<< element<<"\t"<<endl;
				}
			}
		}
	}
	
	/*Check the output: 3*/
	for (int i=0;i<num_of_vertices;i++) {
		if (matrix_processor[i] == -1 ) {
			cout<<"\t ERROR vertex \t" <<i <<"\t is assigned in -1 . (processor) \t"<<endl;
		}
		if (matrix_superstep[i] == -1 ) {
			cout<<"\t ERROR vertex  \t" <<i <<"\t is assigned in -1 . (superstep)\t"<<endl;
		}
	}
	
	/*Print the output*/
	for (int i=0;i<num_of_vertices;i++) {
		/*NODE PROCESSOR SUPERSTEP*/
		out_file_two<< i<<"\t"<<matrix_processor[i]<<"\t"<<matrix_superstep[i]<<"\n";
	}
	
	
	
}

/*Transform the input HyperDAG to a common DAG*/
void hyperDag:: transform_hyperdag_to_dag () {
	matrix_processor =  new int [num_of_vertices+1]; // in which processor is executed every vertex
	matrix_superstep =  new int [num_of_vertices +1]; // in which superstep is executed every vertex
	incoming_vertices_destinations = new list<int> [num_of_vertices +1];
	
	
	for (int i=0;i<=num_of_vertices;i++) matrix_processor[i] = matrix_superstep [i] = -1;
	destination_vertices = new list <int>[num_of_vertices +1];
	
	incoming_edges = new int [num_of_vertices]; //count the incoming edges of every vertex
	outgoing_edges = new int [num_of_vertices]; //count the outgoing edges of every vertex
	for (int i=0;i<num_of_vertices;i++) incoming_edges[i] = outgoing_edges [i] = 0;

	
	int matrix_hyperdags[num_of_pins][2];
	for(int i=0; i<num_of_pins; i++) {
		matrix_hyperdags[i][0] = -1;
		matrix_hyperdags[i][1] = -1;
	}
	
	//store the hyperdag in the matrix_hyperdags[hyperedge][Vertex]
	int i, j;
	int p = 0;
	
	for(i=0; i<num_of_pins; i++) {
		for(j=0;j<2;j++) {
			matrix_hyperdags[i][j] = input_file[p];
			p++;
		}
	}
	
	//find the largest id of hyperedge in order to build the sourcve_vertex matrix
	int hyperedge_id = 0;
	for(i=0; i<num_of_pins; i++) {
			if(matrix_hyperdags[i][0] > hyperedge_id)
			hyperedge_id = matrix_hyperdags[i][0] ;
	}
	
	int *hyperedges_source = new int [hyperedge_id+1];
	for (int i = 0;i<=hyperedge_id;i++) hyperedges_source[i] = -1;
	
	int *source_vertex = new int [num_of_vertices];
	for (int i = 0;i<num_of_vertices;i++) source_vertex[i] = 0;

	
	//identifying the my_var.t_ordering and the destination vertices
	for(i=0; i<num_of_pins; i++) {
		if (hyperedges_source[matrix_hyperdags[i][0]] == -1) {
			//add a source vertex-- no source yet
			hyperedges_source[matrix_hyperdags[i][0]] = matrix_hyperdags[i][1];
			source_vertex[matrix_hyperdags[i][1]]++;
			
		}
		else {
			//count the destination vertices#
			destination_vertices[hyperedges_source[matrix_hyperdags[i][0]]].push_back(matrix_hyperdags[i][1]);
		}	
	}
	
	
	//construct the DAG -- build matrix[i][j] i-->j
	p = 0;
	for (int i = 0; i<num_of_vertices; i++) {
		if (source_vertex[i]!= 0){
			//second print the destinations
			for (auto element : destination_vertices[i]) {
				
				p++;
			}
		}		 
	}
	
	DAG_edges = p;
	int matrix_edges[DAG_edges][2];
	
	for(i=0; i<DAG_edges; i++) {
		matrix_edges[i][0] = -1;
		matrix_edges[i][1] = -1;
	}
	
	int k =0;
	for (int i = 0; i<num_of_vertices; i++) {
		if (source_vertex[i]!= 0){
			//second print the destinationspr
			for (auto element : destination_vertices[i]) {
				//fprintf(stdout,"%d %d \n", i, element );
				incoming_vertices_destinations[element].push_back(i);
				matrix_edges[k][0] = i;
				matrix_edges[k][1] = element;
				k++;
			}
		}		 
	}
	
	if (k!=p){
		printf("\n ERROR:1 !\n" ); 
		exit(1);
	}  
	
	for(int i=0; i<DAG_edges; i++) {
		incoming_edges[matrix_edges[i][1]] ++ ; 
		outgoing_edges[matrix_edges[i][0]] ++ ;
	}


}

void hyperDag::readHyperDag(ifstream& infile, fstream& outfile){
	int x, y, p = 0;
	int param_g, param_l;
	int node_id, workload_cost, comm_delay_cost;
	
	/*open a file in read mode*/
	string line;
	getline(infile, line);
        while(!infile.eof() && line.at(0)=='%')
           getline(infile, line);
	
	
        sscanf(line.c_str(), "%d %d %d", &num_of_hyperedges, &num_of_vertices, &num_of_pins);
        outfile << "%HyperDAG model of naive implementation of 7 iterations of the conjugate gradient method\n";
	outfile << "%The probability for having a nonzero in any cell of matrix A is 0.25\n";
        outfile<<num_of_hyperedges <<"\t"<<num_of_vertices << "\t" << num_of_pins << "\n";
      
        input_file = new int [2*num_of_pins + 1];

        if(num_of_hyperedges<=0 || num_of_vertices<=0 || num_of_pins<=0){
       		cout<<"ERROR Incorrect sizes.\n";
       		exit(1);
        }
        
	/*Read the hyperdag*/
	for(int i=0; i<num_of_pins; ++i){
		getline(infile, line);
		while(!infile.eof() && line.at(0)=='%')
			getline(infile, line);

		sscanf(line.c_str(), "%d %d", &x, &y);
		outfile << x<<"\t"<<y<<"\n";
		input_file[p++] = x;
		input_file[p++] = y;
	}
	
	if((p/2) > num_of_pins){
		cout<<"The number of pins in the header is smaller than the actual number of pins (the header says "<< num_of_pins<< "but already read" << p/2<<"\n";
		exit(1);
	}
	
	/*read the weights and the communication delay cost*/
	workload_weight = new int[num_of_vertices + 1];
	communication_delay_cost = new int [num_of_vertices + 1];
	outfile << "%Node Weight Data\n";
	for (int i =0;i<num_of_vertices;i++){
		if(infile.eof()) {
			cout<<"ERROR: file terminated too early).\n";
			exit(1);
           	}		

		getline(infile, line);
		while(!infile.eof() && line.at(0)=='%')
			getline(infile, line);
		sscanf(line.c_str(), "%d %d %d",&node_id, &workload_cost, &comm_delay_cost);
		workload_weight[node_id] = workload_cost;
		communication_delay_cost[node_id] = comm_delay_cost;
		outfile << node_id<<"\t"<< workload_cost<<"\t"<< comm_delay_cost<<"\n";
	}
	
	/* BSP parameters*/
	outfile<<"% BSP Data \n";
	getline(infile, line);
	while(!infile.eof() && line.at(0)=='%')
	   getline(infile, line);
	   
	sscanf(line.c_str(), "%d %d %d", &pr, &param_g, &param_l);
	outfile<<pr<< "\t"<< param_g<< "\t"<<param_l<<"\n";
	
	outfile<<"% NUMA Data \n";
	for (int i =0;i<pr*pr;i++){
		if(infile.eof()) {
			cout<<"ERROR: file terminated too early).\n";
			exit(1);
           	}		

		getline(infile, line);
		while(!infile.eof() && line.at(0)=='%')
			getline(infile, line);
					
		sscanf(line.c_str(), "%d %d %d", &n_pr1, &n_pr2, &n_cost);
		outfile<<n_pr1<< "\t"<<n_pr2<<"\t"<<n_cost<<"\n";
	}
	infile.close();
	
	 

}

int main(int argc, char *argv[]) {

	if (argc != 2) {
		printf("\n usage: %s <input file> \n\n", argv[0]);
		exit(-1);
	}
	string file_input= string(argv[1]);
	ifstream infile(file_input);
	
	fstream outfile;
	string s1 = "source3_weights_cluster_" + file_input;
	outfile.open("source3_weights_cluster_" + file_input, ios::out);
	
	if(!infile.is_open()){
		cout<<"Error opening the file:	"<<file_input<<"\n";
		exit(1);
	}
	if(!outfile.is_open()){
		cout<<"Error creating the outputfile:	"<<s1<<"\n";
		exit(1);
	}
	
	hyperDag a;
	
	a.readHyperDag(infile, outfile);
	a.transform_hyperdag_to_dag();
	a.scheduling(outfile);
	
	outfile.close();
	
	return 0;
	
	
}
