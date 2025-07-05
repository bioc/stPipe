// [[Rcpp::depends(Rhdf5lib, Rhtslib, Rcpp)]]
//' @name RunDemultiplex
//' @title Demultiplex Stereo-seq data
//' @description This function processes sequencing data for spatial transcriptomics.
//' @param read_1_fq Path to the first FASTQ file.
//' @param read_2_fq Path to the second FASTQ file.
//' @param h5_mapping Path to the HDF5 barcode mapping file.
//' @param output_fq Path to the output FASTQ file.
//' @param n_reads Number of reads to process.
//' @param bc_start Start position of barcode.
//' @param bc_len Length of barcode.
//' @param umi_start Start position of UMI.
//' @param umi_len Length of UMI.
//' @param bin_size Binning size of n * n.
//' @return None. Writes the demultiplexed FASTQ file to the specified path.
//' @export

#include <string>
#include <map>
#include <set>
#include <vector>
#include <unordered_map>
#include <utility>
#include <zlib.h>
#include <htslib/kseq.h>
#include <htslib/hts.h>
#include "H5Cpp.h"
#include "progressbar.h"
#include <Rcpp.h>

using Rcpp::Rcout;
using Rcpp::Rcerr;

#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

// progressbar.h now uses Rcerr internally

// helpers
unsigned int CoordPairToInt(unsigned int a, unsigned int b) {
    return (a << 16) + b;
}

KSEQ_INIT(gzFile, gzread);

unsigned long seq_to_int(const std::string& seq) {
    static const std::map<char, unsigned long> base_to_int_map{
        {'A',0}, {'C',1}, {'T',2}, {'G',3}
    };
    unsigned long acc = 0;
    for (int i = seq.length() - 1; i >= 0; --i) {
        acc += base_to_int_map.at(seq[i]) << (i*2);
    }
    return acc;
}

// ReadStats class
class ReadStats {
public:
    int total, clean_hit, miss, err_hit;
    ReadStats();
    double getHitAcc();
    void reportStats();
    void reportStatsR();
};

ReadStats::ReadStats()
    : total(0), clean_hit(0), miss(0), err_hit(0) {}

double ReadStats::getHitAcc() {
    return (1.0 - double(miss)/double(total));
}

void ReadStats::reportStats() {
    // switch std::cout -> Rcout
    Rcout << "Total number of reads: " << total << "\n";
    Rcout << "Clean hits: "           << clean_hit << "\n";
    Rcout << "Single Error hits: "    << err_hit   << "\n";
    Rcout << "Misses: "               << miss      << "\n";
    Rcout << "Total hits: "           << (clean_hit + err_hit) << "\n";
}

void ReadStats::reportStatsR() {
    // same for R interface
    Rcout << "Total number of reads: " << total << "\n";
    Rcout << "Clean hits: "           << clean_hit << "\n";
    Rcout << "Single Error hits: "    << err_hit   << "\n";
    Rcout << "Misses: "               << miss      << "\n";
    Rcout << "Total hits: "           << (clean_hit + err_hit) << "\n";
}

//' @export
// [[Rcpp::export]]
void RunDemultiplex(
    const char* read_1_fq_path,
    const char* read_2_fq_path,
    const char* h5_mapping_path,
    const char* output_fq_path,
    int n_reads,
    int coord_bc_start,
    int coord_bc_len,
    int umi_start,
    int umi_len,
    int bin_size = 1
) {
    gzFile fp1 = gzopen(read_1_fq_path, "r");
    gzFile fp2 = gzopen(read_2_fq_path, "r");
    gzFile fp_write = gzopen(output_fq_path, "w");
    kseq_t* seq1 = kseq_init(fp1);
    kseq_t* seq2 = kseq_init(fp2);

    H5File file(h5_mapping_path, H5F_ACC_RDONLY);
    DataSet dataset = file.openDataSet("bpMatrix_1");
    DataSpace dataspace = dataset.getSpace();
    hsize_t dims[3];
    dataset.getSpace().getSimpleExtentDims(dims, nullptr);
    unsigned long dim_x = dims[0], dim_y = dims[1], dim_z = dims[2];

    Rcout << "x:" << dim_x << " y:" << dim_y << " z:" << dim_z << "\n";

    std::vector<unsigned long> buf(dim_x*dim_y*dim_z);
    dataset.read(buf.data(), PredType::NATIVE_ULONG);

    std::unordered_map<unsigned long,unsigned int> barcode_map;
    std::set<unsigned long> dupes;
    int c=0;
    progressbar pb; // uses Rcerr
    Rcout << "building map..." << "\n";

    for (auto b : buf) {
        if (b!=0) {
            if (barcode_map.count(b)) {
                barcode_map[b]=CoordPairToInt(0,0);
                dupes.insert(b);
            } else {
                barcode_map[b]=CoordPairToInt(c%dim_x, c/dim_x);
            }
        }
        ++c;
    }
    Rcout << "done!\n";

    ReadStats stats;
    Rcout << "beginning deconvolution..." << "\n";

    for (int i=0; i < n_reads; ++i) {
        int l1 = kseq_read(seq1);
        int l2 = kseq_read(seq2);
        if (l1<0 || l2<0) break;

        std::string s = seq1->seq.s;
        std::string trimmed = s.substr(coord_bc_start, coord_bc_len);
        std::string umi = s.substr(umi_start, umi_len);
        unsigned long code = seq_to_int(trimmed);
        unsigned long masked = code;
        int mask_i=0;
        int seq_len = trimmed.size()*2;
        while (!barcode_map.count(masked) && mask_i<seq_len) {
            masked = code ^ (1<<mask_i++);
        }
        if (mask_i<seq_len) {
            if (mask_i==0) stats.clean_hit++;
            else          stats.err_hit++;
            auto coord = barcode_map[masked];
            std::string x = std::to_string((coord>>16)/bin_size);
            std::string y = std::to_string((coord&0xFFFF)/bin_size);
            std::string header = "@" + trimmed + "_" + umi + "#" + seq2->name.s
                              + ":x" + x + ":y" + y + "\n";
            std::string out = header
                            + seq2->seq.s + "\n+\n"
                            + seq2->qual.s + "\n";
            gzwrite(fp_write, out.c_str(), out.size());
        } else {
            stats.miss++;
        }
        stats.total++;
    }

    stats.reportStats();
    Rcout << "Demultiplexing completed.\n";

    kseq_destroy(seq1);
    kseq_destroy(seq2);
    gzclose(fp1);
    gzclose(fp2);
    gzclose(fp_write);
}

