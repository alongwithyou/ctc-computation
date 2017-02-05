#ifndef _CTC_INTERNAL_H_
#define _CTC_INTERNAL_H_

#include <vector>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cstdint>
#include <cmath>
#include <fstream>
#include <ctime>
#include <algorithm>
#include <random> //This is very useful include
// some boost lib include
#include <boost/xpressive/xpressive.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
using namespace std;



typedef std::vector<std::vector<int> > SimpleDAG;
typedef std::vector<int> IntVector;
namespace CTCFunctions{
	template <typename DataT>
	class RandomGenerator{
	public:
		RandomGenerator(int64_t rows, int64_t cols, int64_t stride=100):
		  rows_(rows), cols_(cols), stride_(stride){
			  stride_ = RoundUp2NearestTwoExponential(cols_);
			  assert(rows_ >= 0 || 
				  cols_ >= 0 ||
				  stride_ >= 0);


			  int64_t mem_len(rows_*stride_);
			  data_ = new DataT[mem_len];
			  assert(data_ != 0);
		  };
		  ~RandomGenerator(){
			  if (data_)
			  {
				  delete [] data_;
				  data_ = 0;
			  }
		  };

		  void SetOutputFile(const std::string& output_file){
			 log_file_name_.assign(output_file);
		  }
		  const std::string GetOutputFile() const {
			  return log_file_name_;
		  }
		void Init(){
			std::mt19937 rng(time(0));
			int64_t total_number = rows_*cols_;
			std::vector<DataT> tmp_result;
			tmp_result.resize(total_number);
			for (int i = 0;i < total_number;i++)
			{
				tmp_result[i] = rng()/(double)(INT_MAX);
			}
			for (int j=0;j<rows_;j++)
			{
				for (int k=0;k<cols_;k++)
				{
					data_[stride_*j+k]= tmp_result[cols_*j+k];
				}
			}
		};
		std::vector<DataT> GetRandom(int64_t row_index, int64_t row_offset,
			int64_t col_index, int64_t col_offset){
				std::vector<DataT> results;
				results.resize(0);
				if (row_index > rows_ || col_index > cols_ ||
					(row_index+row_offset) > rows_ ||
					(col_index+row_offset) > cols_)
				{
					return results;
				}
				int64_t total_random_num = row_offset*col_offset;
				assert(total_random_num >= 0);
				results.resize(total_random_num);
				for (int i=row_index;i<(row_index+row_offset) && 
					i<rows_ ;i++)
				{
					for (int j = col_index; j<(col_index+col_offset) && j<cols_;j++)
					{
						int64_t dst_index = (i-row_index)*col_offset+j-col_index;
						int64_t src_index = i*stride_+j;
						results[dst_index] = data_[src_index];
					}
				}
				return results;

		};
		void Destrtoy(){
			if (data_)
			{
				delete [] data_;
				data_ = 0;
			}
			rows_ = -1;
			cols_ = -1;
			stride_ = -1;
		};
		inline int64_t Rows(){
			return rows_;
		}
		inline int64_t Cols(){
			return cols_;
		}
		inline void ResetRows(int64_t rows){
			rows_ = rows;
		}
		inline void ResetCols(int64_t cols){
			cols_ = cols;
			if (stride_<0)
			{
				stride_ = RoundUp2NearestTwoExponential(cols_);
			}
		}
		typedef struct _SubMatrix{
			DataT* sub_data_;
			int64_t rows_;
			int64_t rows_offset_;
			int64_t cols_;
			int64_t cols_offset_;
			int64_t stride_;
			_SubMatrix(int64_t rows, int64_t rows_offset, int64_t cols, 
				int64_t cols_offset, int64_t stride):rows_(rows),
			rows_offset_(rows_offset),cols_(cols), cols_offset_(cols_offset),
			stride_(stride)
			{
			}

		}SubMatrix;
		SubMatrix& GetSubData(int64_t rows, int64_t rows_offset, int64_t cols, 
			int64_t cols_offset, int64_t stride){
				SubMatrix sub_mat(rows, rows_offset, cols, cols_offset, stride);
				assert(rows >= 0 || rows_offset >= 0 || 
					cols >= 0 || cols_offset >= 0);
				if (rows > rows_ || cols > cols_)
				{
					std::cout << "Cols and Rows set error\n";
					return sub_mat;
				}

				return sub_mat;
				



		}

		int64_t Write(const std::string& output_file, bool is_binary){
			if (output_file.empty())
			{
				std::cout << "Output file name is empty\n";
				return -1;
			}
			std::ofstream ostrm(output_file.c_str(), 
				ios_base::out);
			int64_t succ_num = 0;
			if (ostrm.is_open())
			{
				ostrm << "[ \n";
				for (int i=0;i<rows_;i++)
				{

					if (is_binary)
					{
						ostrm.write(reinterpret_cast<char*>(data_), sizeof(DataT)*cols_);
						ostrm << "\n";
					}else{
						for (int c = 0;c<rows_;c++)
						{
							for (int k = 0;k<cols_;k++)
							{
								ostrm << data_[c*stride_+k] << " ";
							}
							ostrm << std::endl;
						}
					}
					succ_num += cols_;
				}
				ostrm << "\n]\n";
			}

			ostrm.close();
			return succ_num;
		}

		int64_t Read(const std::string& input_file, bool is_binary){
			using namespace boost;
			int64_t succ_read = 0;
			if (input_file.empty())
			{
				std::cout << "File name is empty, can not be read\n";
				"System exits\n";
				return -1;
			}
			Destrtoy();
			std::vector<DataT> new_vector;
			std::ifstream istrm(input_file.c_str(), ios_base::in | ios_base::binary = is_binary);
			if (istrm.is_open())
			{
				int64_t succ_data_line = 0;
				std::string param_line;
				std::vector<std::string> params;
				while(std::getline(istrm, param_line, '\n'))
				{
					if (param_line.empty())
						continue;
					if (ExpectToken(param_line, "[") || ExpectToken(param_line, "]"))
					{
						continue;
					}
					
				}
				split(params, param_line, is_any_of(" \t"), token_compress_off);
				if (params.size())
				{
					
					new_vector.resize(new_vector.size()+params.size());
					for ( int i=0;i<params.size();i++)
					{
						new_vector.push_back(lexical_cast<DataT>(params[i]));
					}
					succ_data_line++;
					if (Cols() < 0)
					{
						ResetCols(params.size());
					}
				}

				if (Rows()<0)
				{
					ResetRows(succ_data_line);
				}
				Init();

			}

			istrm.close();
			return succ_read;
		}
	protected:
		void CopyVectorFrom(const std::vector<DataT>& src_v){
			if (src_v.size() != Rows()*Cols())
			{
				std::cout << "Memory allocated error\n";
				return;
			}
			for (int i=0;i<src_v.size();i++)
			{
				;
			}
		}
		std::ofstream& Write2Stream(std::ofstream& ostrm, 
			const std::vector<DataT>& data, bool is_binary = false){
			if (is_binary == false)
			{
				for (int i=0;i<data.size();i++)
				{
					ostrm << data[i] << " ";
				}
				ostrm << std::endl;
			}
			return ostrm;
		}
		bool ExpectToken(const std::string& token_line, const std::string& token){
			if(token_line.empty())
				return false;
			if(token.empty()){
				std::cout << "token is empty\n";return false;
			}
			std::string clean_token(token);
			TrimSymbols(clean_token);
			std::string::size_type token_pos = token_line.find_first_of(clean_token.c_str());
			if (token_pos != std::string::npos)
			{
				return true;
			}
			return false;

		}
		void TrimSymbols(std::string& token){
			if (token.empty())
			{
				return;
			}
			std::string clean_token;
			unsigned t = 0;
			while(token[t] != ' ')
				clean_token.push_back(token[t]);
			token.assign(clean_token);
		}
		SubMatrix& GetSubMatrix(const SubMatrix& sub_matrix) const
		{
			SubMatrix sub_matrix;

			return sub_matrix;
		}
		int64_t RoundUp2NearestTwoExponential(int64_t num){
			int64_t up_num =num;
			int64_t remain = up_num >> 1, exponential = 0;
			while(remain >= 1){
				remain = remain >> 1;
				exponential++;
			}
			exponential++;
			return std::pow(2.0, static_cast<double>(exponential));//here we can get the right result 
		}

	private:
		DataT* data_;
		int64_t rows_;
		int64_t cols_;
		int64_t stride_;
		std::string log_file_name_;

	};




	typedef struct CTCNode{
		int tindex;
		int uindex;
		float fweight;
		float bweight;
		bool is_log;
		CTCNode* prev;
		CTCNode* next;
		std::vector<CTCNode*> succs;
		CTCNode():tindex(-1),uindex(-1),fweight(0),bweight(0), is_log(true){
			succs.resize(1);
			prev = 0;
			next = 0;
		}

	}Cnode;

	



	class CTCLoss{
		
	public:
		CTCLoss();
		CTCLoss(int64_t target_num):targetNum_(target_num),
		end_stake_(0){
			seq_length_ = 10;
			trellis_.reserve(seq_length_);
			curr_seq_length_ = 0;
		};
		~CTCLoss(){
			if (trellis_.size())
			{
				for(int i=0;i<trellis_.size();i++){
					if (trellis_[i] != nullptr)
					{
						delete [] trellis_[i];
						trellis_[i] = 0;
					}
					
				}
			}
		};
		void SetSeqenceLength(int64_t seq_len){
			seq_length_ = seq_len;
		};
		int64_t GetSequenceLength(){
			return seq_length_;
		}

		void InitTrellis(){
			if(targetNum_ <= 0 || seq_length_ <= 0){
				std::cout << "CTC computation object has not been set \n";
				return;
			}
			trellis_.resize(seq_length_);
			int trellis_width = 2*targetNum_+1;
			for (int i=0;i<seq_length_;i++)
			{
				Cnode* tmp = new Cnode[trellis_width];
				trellis_[i] = tmp;
				Cnode* next = tmp;
				for (int j=0;j<trellis_width-1;j++)
				{
					tmp[j].next = &(tmp[j+1]);
					if(j==0){
						tmp[j].prev = 0;
					}else{
						tmp[j].prev = &(tmp[j-1]);
					}
				}
				tmp[trellis_width-1].next = end_stake_;
				tmp[trellis_width-1].prev = &(tmp[trellis_width-2]);

			}

				
		}
		float ComputeCTCObj(bool is_log = false);
		void ComputeGradient();
	protected:
		float ComputeForward();
		float ComputeBackward();

	private:
		std::vector<Cnode*> trellis_;
		Cnode* end_stake_;
		

		int targetNum_;
		int seq_length_;
		int curr_seq_length_;
	};
};//end of name space CTCFunctions


int NumofNodes(const SimpleDAG& graph);
void DfsVisit(const SimpleDAG& graph, IntVector& find_times, IntVector& completed_times);
void DfsByNode(const SimpleDAG& graph, int node_index, IntVector& find_times, IntVector& completed_times);
void PrintLog(const SimpleDAG& mat);
void FindingSimpleSccs(const SimpleDAG& graph, SimpleDAG& sccs);
void DoSccFinding(const std::string& directed_graph);
void FindSccsFromNodeGraph(const CTCFunctions::Cnode* graph, SimpleDAG& sccs);
void FindSccsWithTarjan(const SimpleDAG& dag_graph, SimpleDAG& unordered_sccs);
int64_t ReadGraph(const std::string& graph_file, std::vector<std::vector<int> >& linked_table);
void FindSccs(CTCFunctions::Cnode* graph, int node_size);
std::vector<float> GenerateRandomParams(int64_t rows, int64_t cols);
void GenerateRandomRational(const std::string ouput_file, int64_t rows, int64_t cols);
void DoCTC_Loss_Computation(int argc, char* argv[]);
CTCFunctions::Cnode* ConvertLinkTable2CnodeGraph(const SimpleDAG& link_table);
void ConvertRawGraph2LinkTable(const SimpleDAG& raw_graph, SimpleDAG& link_table);
#endif