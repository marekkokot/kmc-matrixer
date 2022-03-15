#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cinttypes>
#include <limits>
#include <memory>
#include "kmer_api.h"
#include "kmc_file.h"

using kmer_t = CKmerAPI;

class KMCFileWrapper
{
    std::unique_ptr<CKMCFile> kmc_file;
    size_t tot_kmers;
    size_t cur_kmer_no{};
    kmer_t cur;
    size_t cur_count;
    uint32_t k;
public:
    KMCFileWrapper(KMCFileWrapper&&) = default;
    KMCFileWrapper& operator=(KMCFileWrapper&&) = default;

    KMCFileWrapper(const std::string& path)
    {
        kmc_file = std::make_unique<CKMCFile>();
        if(!kmc_file->OpenForListing(path))
        {
            std::cerr << "Error: cannot open kmc database " << path << "\n";
            exit(1);
        }
        if(kmc_file->IsKMC2() == true)
        {
            std::cerr << "Error: kmc database not sorted: " << path << "\n";
            std::cerr << "Hint: KMC database may be sorted with kmc_tools transform <input> sort <output>\n";
            exit(1);
        }
        CKMCFileInfo kmc_file_info;
        kmc_file->Info(kmc_file_info);
        k = kmc_file_info.kmer_length;
        tot_kmers = kmc_file_info.total_kmers;
        
        cur = kmer_t(k);
        
        if(!Finished())
            Next();        
    }
    uint32_t GetK() const
    {
        return k;
    }
    bool Finished()
    {        
        return cur_kmer_no >= tot_kmers;
    }
    const kmer_t& First() const
    {
        return cur;
    }
    size_t FirstCount() const
    {
        return cur_count;
    }
    void Next()
    {
        ++cur_kmer_no;
        if (!kmc_file->ReadNextKmer(cur, cur_count))
        {
            std::cerr << "Error: critical, this should not happen, details: " << __FILE__ << "(" << __LINE__ << ")\n";
            exit(1);
        }
    }
    ~KMCFileWrapper() noexcept
    {
        if (kmc_file)
            kmc_file->Close();
    }
};

bool allKAreSame(const std::vector<KMCFileWrapper>& samples)
{
    if (samples.empty())
        return true;
    uint32_t k = samples.front().GetK();
    for (const auto& sample : samples)
        if (k!= sample.GetK())
            return false;
    return true;
}


int main(int argc, char**argv)
{
    
    if (argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " <file_with_sorted_kmc_dbs_per_line>\n";
        return 1;
    }
    std::ifstream in(argv[1]);
    if (!in)
    {
        std::cerr << "Error: cannot open file " << argv[1] << "\n";
        return 1;
    } 

    std::string path;
    std::vector<KMCFileWrapper> samples;
    while (std::getline(in, path))
        samples.emplace_back(path);

    if (!allKAreSame(samples))
    {
        std::cerr << "Error: each database should have the same k\n";
        return 1;
    }
    std::string str_kmer;
    while (true)
    {        
        size_t min_id = std::numeric_limits<size_t>::max();
        for (size_t i = 0 ; i < samples.size() ; ++i)
        {
            if (!samples[i].Finished())    
            {                
                if (min_id == std::numeric_limits<size_t>::max() || samples[i].First() < samples[min_id].First())
                    min_id = i;                
            }
        }
        if (min_id == std::numeric_limits<size_t>::max()) // no more k-mers
            break;

        
        //todo: print kmer as string, and tab
        auto min_kmer = samples[min_id].First();
        min_kmer.to_string(str_kmer);
        std::cout << str_kmer;
        for (auto& sample : samples)
        {
            size_t c;
            if (sample.Finished() || !(sample.First() == min_kmer))
                c = 0;
            else
            {
                c = sample.FirstCount();
                sample.Next();
            }
                
            std::cout << "\t" << c;
        }
        std::cout << "\n";
        
    }   
    return 0;
}