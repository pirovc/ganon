#include <seqan/seq_io.h>
#include <seqan/kmer/kmer_base.h>
#include <seqan/kmer/kmer_ibf.h>
#include <seqan/kmer/filtervector.h>
#include <utils/safequeue.hpp>
#include <cxxopts.hpp>
#include <mutex>
#include <vector>
#include <future>
#include <map>
#include <chrono>
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace seqan;

static const uint32_t filterMetadataSize = 256;
static const uint64_t gbInBits = 8589934592;

struct SeqBin{
    CharString id;
    Dna5String seq;
    uint64_t bin;
};

int main(int argc, char* argv[]){
    int old_argc = argc; // parser always set argc to 1
    cxxopts::Options options("ganon-update", "Ganon bloom filter update");
    options.add_options()
        ("e,seqid-bin", "Seqid bin file", cxxopts::value<std::string>())
        ("b,bloom-filter", "Bloom filter file", cxxopts::value<std::string>())
        ("o,output-file", "Alternative output file (default same as bloom-filter)", cxxopts::value<std::string>())
        ("complete", "Old and new sequences are provided for updated bins", cxxopts::value<bool>()->default_value("false"))
        ("t,threads", "Number of threads", cxxopts::value<int>()->default_value("1"))
        ("h,help", "Print help")
        ("v,version", "Show version")
        ("references", "references", cxxopts::value<std::vector<std::string>>())
    ;
    options.parse_positional({"references"});
    options.positional_help("reference.fna[.gz] [reference2.fna[.gz] ... referenceN.fna[.gz]]");

    auto args = options.parse(argc, argv);

    if(args.count("help") || old_argc<=1){
        std::cerr << options.help() << std::endl;
        return 0;
    }else if(args.count("version")){
        std::cerr << "version" << std::endl;
        return 0;
    }

    std::cerr << "seqid-bin: " << args["seqid-bin"].as<std::string>() << std::endl;
    std::cerr << "bloom-filter: " << args["bloom-filter"].as<std::string>() << std::endl;
    std::cerr << "output-file: " << args["output-file"].as<std::string>() << std::endl;
    std::cerr << "threads: " << args["threads"].as<int>() << std::endl;
    std::cerr << "references: " << std::endl;
    for (const auto& s : args["references"].as<std::vector<std::string>>()) {
        std::cerr << s << std::endl;
    }

    int threads = args["threads"].as<int>();

    auto start = std::chrono::high_resolution_clock::now();
    KmerFilter<Dna5, InterleavedBloomFilter, Uncompressed> filter;
    retrieve(filter, toCString(args["bloom-filter"].as<std::string>()));
    uint64_t number_of_bins = getNumberOfBins(filter);
    uint64_t kmer_size = getKmerSize(filter);  
    std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - start;
    std::cerr << "Loading Bloom filter: " << elapsed.count() << std::endl;

    std::map<std::string, uint64_t> bins;
    std::ifstream infile(args["seqid-bin"].as<std::string>());
    std::string seqid;
    uint64_t bin;
    std::unordered_set<uint64_t> updated_bins;
    while (infile >> seqid >> bin){
        bins[seqid] = bin;
        updated_bins.insert(bin);   
    }
    std::cerr << bins.size() << " sequences on " << updated_bins.size() <<  " updated bins (out of " << number_of_bins << " bins)" << std::endl;

    // Reset bins if complete set of sequences is provided (re-create updated bins)
    if(args.count("complete")){
        std::vector<uint32_t> ubins;
        ubins.insert(ubins.end(), updated_bins.begin(), updated_bins.end());
        clear(filter, ubins, threads);
    }

    std::mutex mtx;
    SafeQueue<SeqBin> q;
    std::vector<std::future<void>> tasks;
    bool finished = false;
   
    //Start threads async
    start = std::chrono::high_resolution_clock::now();
    for (int taskNo = 0; taskNo < threads; ++taskNo){
        tasks.emplace_back(
            std::async(std::launch::async, [=, &filter, &q, &finished, &mtx] { 
                while(true){
                    SeqBin val = q.pop();
                    if(val.id!=""){ //if not empty
                        insertKmer(filter, val.seq, val.bin, 0);
                        mtx.lock();
                        std::cerr << val.id << " -> k-mers added to bin " << val.bin << std::endl; 
                        mtx.unlock();
                    }
                    if(finished && q.empty())
                        break;
                }
            })
        ); 
    }

    // extra thread for reading the input
    tasks.emplace_back(
        std::async(std::launch::async, [=, &bins, &q, &finished, &mtx] {
            for(auto const& reference_fasta_file: args["references"].as<std::vector<std::string>>()) {
                SeqFileIn seqFileIn;
                if (!open(seqFileIn, toCString(reference_fasta_file))){
                    std::cerr << "Unable to open " << reference_fasta_file << std::endl;
                    continue;
                }
                while(!atEnd(seqFileIn)){
                    while(q.size() > (threads*10)){
                        ; //spin
                    }
                    StringSet<CharString> ids;
                    StringSet<IupacString> seqs;
                    readRecords(ids, seqs, seqFileIn, threads*5);
                    for (uint16_t readID = 0; readID < length(ids); ++readID){
                        if(length(seqs[readID])<kmer_size){ //sequence too small
                        	mtx.lock();
                            std::cerr << ids[readID] << " has sequence smaller than k-mer size" << std::endl; 
                            mtx.unlock();
                            continue;
                        }
                        std::string cid = toCString(ids[readID]);
                        std::string acc = cid.substr(0, cid.find(' '));
                        if(bins.count(acc)==0){ //not defined on the bins
                        	mtx.lock();
                            std::cerr << acc << " not defined on bins file" << std::endl; 
                            mtx.unlock();
                            continue;
                        }
                        q.push(SeqBin{acc, seqs[readID], bins[acc]});
                    }
                }
                close(seqFileIn);
            }
            finished = true;
        })
    ); 

    for (auto &&task : tasks){
        task.get();
    }
    elapsed = std::chrono::high_resolution_clock::now() - start;
    std::cerr << "Adding k-mers: " << elapsed.count() << std::endl;

    start = std::chrono::high_resolution_clock::now();
    if(args.count("output-file"))
        store(filter, toCString(args["output-file"].as<std::string>()));
    else
        store(filter, toCString(args["bloom-filter"].as<std::string>()));
    elapsed = std::chrono::high_resolution_clock::now() - start;
    std::cerr << "Saving Bloom filter: " << elapsed.count() << std::endl;

    return 0;
}