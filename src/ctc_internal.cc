#include "ctc_internal.h"
#include <exception>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <cmath>

namespace CTCFunctions{
	CTCLoss::CTCLoss(){}
	//CTCLoss::~CTCLoss(){}; //bad implimentation


	
}

void FindingSimpleSccs(const SimpleDAG& graph, SimpleDAG& sccs)
{
	if (graph.size() == 0)
	{
		fprintf(stdout, "Graph is empty\n");
		return;
	}
	sccs.resize(1);//bad choice
	for (int i = 0;i < graph.size();i++)
	{
		if (graph[i].size() == 0) //find basic scc
		{
			IntVector& this_scc = sccs.back();
			this_scc.push_back(i);
			sccs.resize(sccs.size()+1);
		}
	}
}


int NumofNodes(const SimpleDAG& graph)
{
	if(graph.size() == 0)
		return -1;
	int min_index = INT_MAX, max_index = -1;
	for (int i=0;i<graph.size();i++)
	{
		for (int j=0;j<graph[i].size();j++)
		{
			if (min_index > graph[i][j] && graph[i][j] != -1)
			{
				min_index = graph[i][j];
			}
			if (max_index < graph[i][j] && graph[i][j] != -1)
			{
				min_index = graph[i][j];
			}

		}
	}
	return std::abs(max_index-min_index+1);
}


void PrintLog(const SimpleDAG& mat)
{
	if(mat.size() == 0)
		return;
	for (int i=0;i<mat.size();i++)
	{
		fprintf(stdout, "The %dth scc : \n", i);
		for (int j=0;j<mat[i].size();j++)
		{
			std::cout << mat[i][j] << " ";
		}
		std::cout << "\n";
	}
}


void DfsByNode(const SimpleDAG& graph, int node_index, IntVector& find_times, IntVector& completed_times)
{

}

void DfsVisit(const SimpleDAG& graph, IntVector& find_times, IntVector& completed_times)
{
	if (graph.size() == 0)
	{
		return;
	}
	//////////////////////////////////////////////////////////////////////////
	int node_num = NumofNodes(graph);
	find_times.resize(node_num);
	completed_times.resize(node_num);
	for (int i=0;i<graph.size();i++)
	{
		const IntVector& succs = graph[i];
		for (int j=0;j<succs.size();j++)
		{
			DfsByNode(graph, succs[j], find_times, completed_times);
		}
	}
}

void DoSccFinding(const std::string& directed_graph){
	using namespace CTCFunctions;
	SimpleDAG graph, sccs, dag_graph;
	ReadGraph(directed_graph, graph);
	ConvertRawGraph2LinkTable(graph, dag_graph);
	Cnode* node_graph = ConvertLinkTable2CnodeGraph(graph);
	PrintLog(graph);
	FindingSimpleSccs(dag_graph, sccs);
	PrintLog(sccs);
}

int64_t ReadGraph(const std::string& graph_file, std::vector<std::vector<int> >& linked_table)
{
	using namespace boost;
	int64_t node_size  = -1;
	
	try
	{
		std::ifstream graph_strm(graph_file.c_str(), ios::in);
		if (graph_strm.is_open())
		{
			std::vector<std::string> nodes;
			std::string line_of_nodes;
			int node_index = 0;
			linked_table.resize(1);
			while(std::getline(graph_strm, line_of_nodes, '\n')){
				if (line_of_nodes.empty() || line_of_nodes[0] == '#')
				{
					continue;
				}
				split(nodes, line_of_nodes, is_any_of(" \t"), token_compress_off);
				std::vector<int>& curr_layer = linked_table.back();
				curr_layer.resize(nodes.size());
				for (int i=0;i<curr_layer.size();i++)
				{
					curr_layer[i] = std::atoi(nodes[i].c_str());
				}
				linked_table.resize(linked_table.size()+1);
			}
		}
		graph_strm.close();
	}
	catch (const std::exception* e)
	{
		std::cout << e->what() << std::endl;
	}
	


	return node_size;
}

void FindSccs(CTCFunctions::Cnode* graph, int node_size){
	if (! graph || node_size <= 0)
	{
		std::cout << "Graph is not exist!\n";
		return;
	}
	if (node_size == 1)
	{
		fprintf(stdout, "The only scc is : %d\n", graph->uindex);
		return;
	}
	std::vector<std::vector<int> > sccs;
	std::vector<std::vector<int> > graph_indices;
	std::vector<int> findTime(node_size, -1), completedTime(node_size, -1);
	
}


std::vector<float> GenerateRandomParams(int64_t rows, int64_t cols)
{
	CTCFunctions::RandomGenerator<float> randomer(rows,cols);
	randomer.Init();
	std::vector<float> params = randomer.GetRandom(0,rows, 0, cols);
	return params;
}

void GenerateRandomRational(const std::string ouput_file, int64_t rows, int64_t cols)
{
	CTCFunctions::RandomGenerator<float> randomer(rows,cols);
	randomer.SetOutputFile(ouput_file);
	randomer.Init();
	std::vector<float> params = randomer.GetRandom(2,5,2,6);
	randomer.Write(randomer.GetOutputFile(), false);
}
void DoCTC_Loss_Computation(int argc, char* argv[])
{
	CTCFunctions::CTCLoss ctc_computation(100);
	ctc_computation.InitTrellis();
}


CTCFunctions::Cnode* ConvertLinkTable2CnodeGraph(const SimpleDAG& link_table)
{
	using namespace CTCFunctions;
	Cnode* node = nullptr;
	if (link_table.empty())
	{
		return node;
	}
	int node_num = NumofNodes(link_table);
	node = new CTCFunctions::Cnode[node_num](); 
	Cnode* next_unused = node;
	for (int i=0;link_table[i].size() > 0 && i<link_table.size();i++)
	{
		int k = 0;
		int root_node_index = link_table[i][k];
		Cnode& curr_root_node = node[root_node_index];
		for (int j=1;link_table[i].size() > 0 && j<link_table[i].size();j++)
		{
			int curr_node_index = link_table[i][j];
			if(curr_node_index != -1)
				curr_root_node.succs.push_back(&(node[curr_node_index]));
		}
	}
	return node;
}


void FindSccsFromNodeGraph(const CTCFunctions::Cnode* graph, SimpleDAG& sccs)
{
	if (graph == 0)
	{
		std::cout << "No graph used here\n";
		return;
	}
	SimpleDAG tmp_sccs;

}


void FindSccsWithTarjan(const SimpleDAG& dag_graph, SimpleDAG& unordered_sccs)
{
	if (dag_graph.size() == 0)
	{
		return;
	}
	unordered_sccs.resize(dag_graph.size()); // this is not the best choice
}

void ConvertRawGraph2LinkTable(const SimpleDAG& raw_graph, SimpleDAG& link_table)
{
	int node_num = NumofNodes(raw_graph);
	link_table.resize(node_num);
	for (int i=0;raw_graph[i].size() > 0 && i<raw_graph.size();i++)
	{
		int root_node_index = raw_graph[i][0];
		for (int j=1;raw_graph[i].size() > 0 && j<raw_graph[i].size();j++)
		{
			int curr_node_index = raw_graph[i][j];
			if(curr_node_index != -1)
				link_table[root_node_index].push_back(curr_node_index);
		}
	}
}