// *********************************************************************
// Extract phase for PS candidates from complex interferograms
// ---------------------------------------------------------------------
// AUTHOR    : Andy Hooper
// ---------------------------------------------------------------------
// WRITTEN   : 11.12.2004
//
// Change History
// ==============================================
// 03/2009 MA Fix for gcc 4.3.x
// 11/2025 Mingjia Li I/O optimizations, increase throughput to decrease seek times
// ==============================================

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cstring>
#include <stdint.h>

using namespace std;

// 阈值设定：2MB
// 机械硬盘 Seek 一次约 10ms，这 10ms 足够顺序读取 1.5MB - 2MB 数据。
// 如果两个点间距小于此值，顺序读比跳跃更快。
const uint64_t READ_VS_SEEK_THRESHOLD = 2 * 1024 * 1024; 

struct PSCandidate {
    uint64_t file_offset; 
    size_t original_idx;  
};

bool compareCandidates(const PSCandidate& a, const PSCandidate& b) {
    return a.file_offset < b.file_offset;
}

int main(int argc, char *argv[]) {
    // 禁用 C++ IO 同步，提速
    std::ios_base::sync_with_stdio(false);
    
    try {
        if (argc < 2) {
            cout << "Usage: pscphase parmfile pscands.1.ij pscands.1.ph" << endl;
            return 1;
        }

        const char *ijname = (argc < 3) ? "pscands.1.ij" : argv[2];
        const char *outfilename = (argc < 4) ? "pscands.1.ph" : argv[3];

        // 1. 读取参数
        ifstream parmfile(argv[1]);
        if (!parmfile.is_open()) throw "Error opening parmfile";

        long width = 0;
        parmfile >> width;
        cout << "Width: " << width << endl;

        vector<string> ifg_filenames;
        char linebuf[1024];
        parmfile.getline(linebuf, 1024);
        while (parmfile >> linebuf) {
            ifg_filenames.push_back(string(linebuf));
        }
        parmfile.close();

        int num_files = ifg_filenames.size();
        cout << "Number of interferograms: " << num_files << endl;

        // 2. 读取并排序 PS 点
        ifstream psfile(ijname);
        if (!psfile.is_open()) throw "Error opening PS file";
        
        cout << "Reading PS candidates list..." << endl;
        
        vector<PSCandidate> candidates;
        candidates.reserve(1000000); 

        long pscid, y, x;
        size_t idx = 0;
        char buffer[1024];
        
        while (psfile >> pscid >> y >> x) {
            psfile.getline(buffer, 1024); 
            PSCandidate cand;
            cand.file_offset = ((uint64_t)y * width + x) * 8; 
            cand.original_idx = idx;
            candidates.push_back(cand);
            idx++;
        }
        psfile.close();
        
        size_t num_ps = candidates.size();
        cout << "Total PS candidates: " << num_ps << endl;

        // 关键步骤：排序
        cout << "Sorting candidates for optimized access..." << endl;
        vector<PSCandidate> sorted_candidates = candidates;
        std::sort(sorted_candidates.begin(), sorted_candidates.end(), compareCandidates);

        // 3. 智能处理
        vector<char> output_buffer(num_ps * 8);
        vector<char> smart_buffer; // 动态缓冲区
        smart_buffer.reserve(READ_VS_SEEK_THRESHOLD * 2); 

        ofstream outfile(outfilename, ios::out | ios::binary);
        if (!outfile.is_open()) throw "Error opening output file";

        for (int i = 0; i < num_files; ++i) {
            ifstream ifg(ifg_filenames[i].c_str(), ios::in | ios::binary);
            if (!ifg.is_open()) {
                cerr << "Warning: Cannot open " << ifg_filenames[i] << endl;
                continue; 
            }

            // Header 处理
            char header[32];
            long magic = 0x59a66a95;
            ifg.read(header, 32);
            if (*reinterpret_cast<long*>(header) != magic) {
                ifg.seekg(0, ios::beg);
            }

            size_t cand_idx = 0;
            
            // 智能聚合循环
            while (cand_idx < num_ps) {
                // 1. 确定当前“聚合块”的起始和结束
                // 我们尝试把一群靠得近的点合并成一次 Read
                uint64_t block_start = sorted_candidates[cand_idx].file_offset;
                uint64_t block_end = block_start + 8;
                
                size_t batch_count = 1;
                
                // 向前探查：如果下个点离得不远，就把它也包进来
                while (cand_idx + batch_count < num_ps) {
                    uint64_t next_pos = sorted_candidates[cand_idx + batch_count].file_offset;
                    uint64_t gap = next_pos - block_end;
                    
                    if (gap < READ_VS_SEEK_THRESHOLD) {
                        block_end = next_pos + 8; // 延伸块的尾部
                        batch_count++;
                    } else {
                        break; // 距离太远，断开，准备 Seek
                    }
                }

                // 2. 执行 I/O
                uint64_t read_size = block_end - block_start;
                
                // 确保缓冲区够大
                if (smart_buffer.size() < read_size) {
                    smart_buffer.resize(read_size * 1.5); // 稍微多给点避免频繁重分配
                }

                // Seek 到块的开头 (如果需要)
                ifg.seekg(block_start, ios::beg);
                
                // 一次性读入整个块 (包含多个 PS 点和中间的缝隙)
                ifg.read(smart_buffer.data(), read_size);

                // 3. 从内存中提取数据
                for (size_t k = 0; k < batch_count; ++k) {
                    const auto& cand = sorted_candidates[cand_idx + k];
                    
                    // 计算该点在 smart_buffer 中的相对位置
                    size_t offset_in_block = cand.file_offset - block_start;
                    
                    // 拷贝到最终输出 buffer
                    memcpy(&output_buffer[cand.original_idx * 8], 
                           &smart_buffer[offset_in_block], 
                           8);
                }

                // 推进指针
                cand_idx += batch_count;
            }

            ifg.close();
            
            // 写出结果
            outfile.write(output_buffer.data(), output_buffer.size());
            cout << "Processed IFG " << (i + 1) << " / " << num_files << endl;
        }

        outfile.close();
        cout << "Done." << endl;

    } catch (const char* msg) {
        cerr << "EXCEPTION: " << msg << endl;
        return 999;
    } catch (...) {
        return 999;
    }

    return 0;
}