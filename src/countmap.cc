#include <getopt.h>
#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <ctime>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <tuple>
#include <utility>
#include <chrono>

#include "countmap.hpp"
#include "fastaq.hpp"
#include "mapping_params.hpp"
#include "index.hpp"
#include "map.hpp"
#include "brown_minimizers.hpp"
#include "bioparser/bioparser.hpp"
#include "thread_pool/thread_pool.hpp"

// double clip_time = 0;
// double minimizer_time = 0;
// double hits_time = 0;
// double radix_time = 0;
// double candidates_time = 0;
// double check_time = 0;
// double circ_time = 0;
// double clip_time_2 = 0;
// double minimizer_time_2 = 0;
// double hits_time_2 = 0;
// double radix_time_2 = 0;
// double candidates_time_2 = 0;
// double check_time_2 = 0;
// double circ_time_2 = 0;
// double region_time = 0;
// double expand_time = 0;
// double revcomp_time = 0;
// double ksw2_time = 0;
// double process_time = 0;
// double process_time_2 = 0;
// double second_time = 0;
// double sam_time = 0;
// double unmapped_time = 0;

// Minimizer: value, position, origin
typedef std::tuple<uint64_t, uint32_t, bool> minimizer_t;
// Index: position, range
typedef std::pair<uint32_t, uint32_t> index_pos_t;
// Paired read
typedef std::pair<std::vector<std::unique_ptr<fastaq::FastAQ>>, std::vector<std::unique_ptr<fastaq::FastAQ>>> paired_reads_t;
// bioparser ptr
typedef std::unique_ptr<bioparser::Parser<fastaq::FastAQ>> parser_ptr_t;
// Paired bioparser ptr
typedef std::pair<std::unique_ptr<bioparser::Parser<fastaq::FastAQ>>, std::unique_ptr<bioparser::Parser<fastaq::FastAQ>>> paired_parser_ptr_t;

// Accepted file formats
const std::unordered_set<std::string> fasta_formats = {".fasta", ".fa", ".fasta.gz", ".fa.gz"};
const std::unordered_set<std::string> fastq_formats = {".fastq", ".fq", ".fastq.gz", ".fq.gz"};

// Sample size in bytes used for insert size inferrence
constexpr uint32_t sample_bytes = 256 * 1024 * 1024;

static struct option long_options[] = {
  {"help", no_argument, NULL, 'h'},
  {"version", no_argument, NULL, 'v'},
  {"paired", no_argument, NULL, 'p'},
  {"infer", required_argument, NULL, 'I'},
  {"all", no_argument, NULL, 'a'},
  {"match", required_argument, NULL, 'm'},
  {"mismatch", required_argument, NULL, 'M'},
  {"gap-open", required_argument, NULL, 'o'},
  {"gap-extend", required_argument, NULL, 'e'},
  {"band", required_argument, NULL, 'b'},
  {"kmer", required_argument, NULL, 'k'},
  {"kmer_s", required_argument, NULL, 'K'},
  {"window_length", required_argument, NULL, 'w'},
  {"window_length_s", required_argument, NULL, 'W'},
  {"frequency", required_argument, NULL, 'f'},
  {"insert_size", required_argument, NULL, 'i'},
  {"st_deviation", required_argument, NULL, 's'},
  {"threshold", required_argument, NULL, 'T'},
  {"threshold_s", required_argument, NULL, 'S'},
  {"threads", required_argument, NULL, 't'},
  {"batch_size", required_argument, NULL, 'B'},
  {NULL, no_argument, NULL, 0}
};


void help(void) {
  printf("CountMap - a tool for mapping short reads to a reference genome.\n\n"

         "Usage: countmap [OPTIONS] reference.fa [read.fq]\n"
         "  reference.fa - FASTA file containing reference genome\n"
         "  [read.fq]    - one or two FASTQ files containing sequencing data\n\n"

         "Supported file extensions: .fasta\n"
         "                           .fa\n"
         "                           .fastq\n"
         "                           .fq\n"
         "                           .fasta.gz\n"
         "                           .fa.gz\n"
         "                           .fastq.gz\n"
         "                           .fq.gz\n\n"
         
         "OPTIONS:\n"
         "  -h  or  --help             show a summary of available options and exit\n"
         "  -v  or  --version          show the current version number and exit\n"
         "  -p  or  --paired           paired-end mode\n"
         "  -I  or  --infer            only output inferred insert size and standard deviation\n"
         "  -a  or  --all              output all found mappings\n"
         "  -m  or  --match            <int>\n"
         "                               default: 2\n"
         "                               match value\n"
         "  -M  or  --mismatch         <int>\n"
         "                               default: 4\n"
         "                               mismatch value\n"
         "  -o  or  --gap-open         <int>\n"
         "                               default: 4\n"
         "                               gap open value\n"
         "  -e  or  --gap-extend       <int>\n"
         "                               default: 2\n"
         "                               gap extend value\n"
         "  -b  or  --band             <int>\n"
         "                               default: -1\n"
         "                               ksw2 alignment band, band < 0 => disabled\n"
         "                               option: -1  => set band to half of QLEN\n"
         "  -k  or  --kmer             <uint>\n"
         "                               default: 12\n"
         "                               constraints: largest supported is 32\n"
         "                               k-mer size\n"
         "  -K  or  --kmer_s           <uint>\n"
         "                               default: 9\n"
         "                               constraints: largest supported is 32\n"
         "                               k-mer size (second round)\n"
         "  -w  or  --window_length    <uint>\n"
         "                               default: 5\n"
         "                               length of window\n"
         "  -W  or  --window_length_s  <uint>\n"
         "                               default: 2\n"
         "                               length of window (second round)\n"
         "  -f  or  --frequency        <float>\n"
         "                               default: 0.01\n"
         "                               constraints: must be from [0, 1]\n"
         "                               number of most frequent minimizers that\n"
         "                               are not taken into account\n"
         "  -i  or  --insert_size      <uint>\n"
         "                               default: 215\n"
         "                               fragment insert size mean\n"
         "                               if set overrides automatic insert size inferrence\n"
         "  -s  or  --st_deviation     <float>\n"
         "                               default: 10.0\n"
         "                               fragment insert size standard deviation\n"
         "  -T  or  --threshold        <uint>\n"
         "                               default: 2\n"
         "                               number of hits needed in order to consider\n"
         "                               a region a candidate for mapping\n"
         "  -S  or  --threshold_s      <uint>\n"
         "                               default: 3\n"
         "                               number of hits needed in order to consider\n"
         "                               a region a candidate for mapping (second round)\n"
         "  -t  or  --threads          <uint>\n"
         "                               default: 3\n"
         "                               number of threads\n"
         "  -B  or  --batch_size       <uint>\n"
         "                               default: 256\n"
         "                               read loading batch size in MB\n"
  );
}

void version(void) {
  printf("countmap %d.%d\n",
    countmap_VERSION_MAJOR,
    countmap_VERSION_MINOR
  );
}

// Check file extension
// Args: filename   - name of file to be checked for extension
//       extensions - set of accepted extensions
// Return: extension accepted or not accepted
bool check_extension(const std::string& filename, const std::unordered_set<std::string>& extensions) {
  for (const auto& it : extensions) {
    if (filename.size() > it.size()) {
      if (filename.compare(filename.size()-it.size(), std::string::npos, it) == 0) {
        return true;
      }
    }
  }
  return false;
}

int main(int argc, char **argv) {
  int optchr;
  mapping_params_t parameters;
  parameters.all = false;
  parameters.mch = 2;
  parameters.mis = 4;
  parameters.gapo = 4;
  parameters.gape = 2;
  parameters.band = -1;
  parameters.k = 12;
  parameters.w = 5;
  parameters.k_2 = 9;
  parameters.w_2 = 2;
  parameters.f = 0.01f;
  parameters.insert_size = 215;
  parameters.threshold = 2;
  parameters.threshold_2 = 3;
  bool paired = false;
  bool infer_is = false;
  bool set_insert = false;
  uint32_t t = 3;
  uint32_t batch_size = 256 * 1024 * 1024;
  std::string cl_flags;

  while ((optchr = getopt_long(argc, argv, "hvpIam:M:o:e:b:k:K:w:W:f:i:s:T:S:t:B:", long_options, NULL)) != -1) {
    cl_flags += "-", cl_flags += optchr, cl_flags += " ";
    if (optarg != nullptr) cl_flags += optarg, cl_flags += " ";
    switch (optchr) {
      case 'h': {
        help();
        exit(0);
      }
      case 'v': {
        version();
        exit(0);
      }
      case 'p': {
        paired = true;
        break;
      }
      case 'I': {
        infer_is = true;
        break;
      }
      case 'a': {
        parameters.all = true;
        break;
      }
      case 'm': {
        parameters.mch = atoi(optarg);
        break;
      }
      case 'M': {
        parameters.mis = atoi(optarg);
        break;
      }
      case 'o': {
        parameters.gapo = atoi(optarg);
        break;
      }
      case 'e': {
        parameters.gape = atoi(optarg);
        break;
      }
      case 'b': {
        parameters.band = atoi(optarg);
        break;
      }
      case 'k': {
        parameters.k = atoi(optarg);
        break;
      }
      case 'K': {
        parameters.k_2 = atoi(optarg);
        break;
      }
      case 'w': {
        parameters.w = atoi(optarg);
        break;
      }
      case 'W': {
        parameters.w_2 = atoi(optarg);
        break;
      }
      case 'f': {
        parameters.f = atof(optarg);
        if (parameters.f < 0.0f || parameters.f > 1.0f) {
          fprintf(stderr, "[countmap] error: f must be from [0, 1].\n"); 
          exit(1); 
        }
        break;
      }
      case 'i': {
        set_insert = true;
        parameters.insert_size = atoi(optarg);
        break;
      }
      case 's': {
        parameters.sd = atof(optarg);
        break;
      }
      case 'T': {
        parameters.threshold = atoi(optarg);
        break;
      }
      case 'S': {
        parameters.threshold_2 = atoi(optarg);
        break;
      }
      case 't': {
        t = atoi(optarg);
        break;
      }
      case 'B': {
        batch_size = atoi(optarg) * 1024 * 1024;
        break;
      }
      default: {
        fprintf(stderr, "[countmap] error: Unknown option. Type %s --help for usage.\n", argv[0]);
        exit(1);
      }
    }
  }

  if (argc - optind > 3 || argc - optind < 2) {
    fprintf(stderr, "[countmap] error: Expected read(s) and reference. Use --help for usage.\n");
    exit(1);
  }
  if (argc - optind == 2 && infer_is) {
    fprintf(stderr, "[countmap] error: Expected paired reads in order to infer insert size (option -I).\n");
    exit(1);
  }
  if (argc - optind == 2 && paired) {
    fprintf(stderr, "[countmap] error: Expected paired reads in order to use paired-end mode (option -p).\n");
    exit(1);
  }


  fprintf(stderr, "[countmap-load] loading reference... ");

  std::string reference_file(argv[optind]);
  if (!check_extension(reference_file, fasta_formats)) {
      fprintf(stderr, "[countmap] error: Unsupported reference file format. Check --help for supported file formats.\n");
      exit(1);
  }
  std::vector<std::unique_ptr<fastaq::FastAQ>> reference;

  parser_ptr_t ref_parser;
  if (check_extension(reference_file, fasta_formats)) {
    ref_parser = bioparser::createParser<bioparser::FastaParser, fastaq::FastAQ>(reference_file);
  } else {
    ref_parser = bioparser::createParser<bioparser::FastqParser, fastaq::FastAQ>(reference_file);
  }
  ref_parser->parse(reference, -1);

  fprintf(stderr, "\r[countmap-load] loaded reference           \n"
                  "[countmap-index] indexing reference... ");

  auto i_start = std::chrono::steady_clock::now();
  std::shared_ptr<thread_pool::ThreadPool> thread_pool = thread_pool::createThreadPool(t);

  std::vector<minimizer_t> t_minimizers = collect_minimizers(reference[0], parameters.w, parameters.k, t);
  prep_ref(t_minimizers, parameters.f);
  std::unordered_map<uint64_t, index_pos_t> ref_index = index_ref(t_minimizers);

  std::vector<minimizer_t> t_minimizers_2 = collect_minimizers(reference[0], parameters.w_2, parameters.k_2, t);
  prep_ref(t_minimizers_2, parameters.f / 10.0);
  std::unordered_map<uint64_t, index_pos_t> ref_index_2 = index_ref(t_minimizers_2);

  fprintf(stderr, "\r[countmap-index] indexed reference        \n");

  fastaq::FastAQ::print_statistics(reference, reference_file);

  auto i_end = std::chrono::steady_clock::now();
  auto i_interval = std::chrono::duration_cast<std::chrono::duration<double>>(i_end - i_start);
  fprintf(stderr, "[countmap-index] index time: %.2f sec\n", i_interval.count());
  
  std::srand(std::time(0));

  if (infer_is && argc - optind == 3) {
    std::string reads_file1(argv[optind + 1]);
    std::string reads_file2(argv[optind + 2]);
    if (!(check_extension(reads_file1, fasta_formats) || check_extension(reads_file1, fastq_formats))
        || !(check_extension(reads_file2, fasta_formats) || check_extension(reads_file2, fastq_formats))) {
      fprintf(stderr, "[countmap] error: Unsupported paired-end reads formats. Check --help for supported file formats.\n");
      exit(1);
    }
    paired_parser_ptr_t parser;
    if (check_extension(reads_file1, fasta_formats) && check_extension(reads_file2, fasta_formats)) {
      parser.first = bioparser::createParser<bioparser::FastaParser, fastaq::FastAQ>(reads_file1);
      parser.second = bioparser::createParser<bioparser::FastaParser, fastaq::FastAQ>(reads_file2);
    } else if (check_extension(reads_file1, fastq_formats) && check_extension(reads_file2, fastq_formats)) {
      parser.first = bioparser::createParser<bioparser::FastqParser, fastaq::FastAQ>(reads_file1);
      parser.second = bioparser::createParser<bioparser::FastqParser, fastaq::FastAQ>(reads_file2);
    } else {
      fprintf(stderr, "[countmap] error: Paired-end reads formats not equal.\n");
      exit(1);
    }

    paired_reads_t paired_reads;
    parser.first->parse(paired_reads.first, sample_bytes);
    parser.second->parse(paired_reads.second, sample_bytes);

    fastaq::FastAQ::print_statistics(paired_reads.first, reads_file1);
    fastaq::FastAQ::print_statistics(paired_reads.second, reads_file2);

    infer_insert_size(ref_index, t_minimizers,
                      reference[0], paired_reads,
                      parameters);
    fprintf(stderr, "[countmap-map] inferred insert size mean, standard deviation: %u, %.2f\n", 
            parameters.insert_size, parameters.sd);
    return 0;
  }

  auto m_start = std::chrono::steady_clock::now();

  if (argc - optind == 3) {
    if (paired) {
      fprintf(stderr, "[countmap] paired-end mode\n");
    }
    else {
      fprintf(stderr, "[countmap] single-end mode\n");
    }

    std::string reads_file1(argv[optind + 1]);
    std::string reads_file2(argv[optind + 2]);
    if (!(check_extension(reads_file1, fasta_formats) || check_extension(reads_file1, fastq_formats))
        || !(check_extension(reads_file2, fasta_formats) || check_extension(reads_file2, fastq_formats))) {
      fprintf(stderr, "[countmap] error: Unsupported paired-end reads formats. Check --help for supported file formats.\n");
      exit(1);
    }
    paired_parser_ptr_t parser;
    if (check_extension(reads_file1, fasta_formats) && check_extension(reads_file2, fasta_formats)) {
      parser.first = bioparser::createParser<bioparser::FastaParser, fastaq::FastAQ>(reads_file1);
      parser.second = bioparser::createParser<bioparser::FastaParser, fastaq::FastAQ>(reads_file2);
    } else if (check_extension(reads_file1, fastq_formats) && check_extension(reads_file2, fastq_formats)) {
      parser.first = bioparser::createParser<bioparser::FastqParser, fastaq::FastAQ>(reads_file1);
      parser.second = bioparser::createParser<bioparser::FastqParser, fastaq::FastAQ>(reads_file2);
    } else {
      fprintf(stderr, "[countmap] error: Paired-end reads formats not equal.\n");
      exit(1);
    }

    printf("@HD\tVN:1.6\n"
           "@SQ\tSN:%s\tLN:%lu\n"
           "@PG\tID:countmap\tPN:countmap\tVN:%d.%d\tCL:%s %s%s %s %s\n", 
           reference[0]->name.c_str(), reference[0]->sequence.size(),
           countmap_VERSION_MAJOR, countmap_VERSION_MINOR, 
           argv[0], cl_flags.c_str(), 
           argv[optind], argv[optind + 1], argv[optind + 2]);

    uint32_t batch_num = 1;
    while (true) {
      fprintf(stderr, "[countmap] batch number %u\n", batch_num++);
      fprintf(stderr, "[countmap-load] loading paired-end reads... ");

      paired_reads_t paired_reads;
      bool status1 = parser.first->parse(paired_reads.first, batch_size);
      bool status2 = parser.second->parse(paired_reads.second, batch_size);
      
      fprintf(stderr, "\r[countmap-load] loaded paired-end reads        \n");
      
      fastaq::stats pr1_stats = fastaq::FastAQ::print_statistics(paired_reads.first, reads_file1);
      fastaq::stats pr2_stats = fastaq::FastAQ::print_statistics(paired_reads.second, reads_file2);

      if (paired) {
        if (pr1_stats.num != pr2_stats.num) {
          fprintf(stderr, "[countmap] error: Paired-end read files must have equal number of reads (pairs).\n");
          exit(1);
        }
        if (pr1_stats.max - pr1_stats.min > 0 || pr2_stats.max - pr2_stats.min > 0) {
          fprintf(stderr, "[countmap] warning: Reads are not of fixed size.\n");
        }
        if (!set_insert) {
          infer_insert_size(ref_index, t_minimizers,
                            reference[0], paired_reads,
                            parameters);
          fprintf(stderr, "[countmap-map] inferred insert size mean, standard deviation: %u, %.2f\n", 
                  parameters.insert_size, parameters.sd);
        }
      }
      
      auto c_start = std::chrono::steady_clock::now();
      std::vector<std::future<std::string>> thread_futures;
      for (unsigned int tasks = 0; tasks < t - 1; ++tasks) {
        thread_futures.emplace_back(thread_pool->submit(paired ? map_paired : map_as_single, 
                std::ref(ref_index), std::ref(ref_index_2), std::ref(t_minimizers), std::ref(t_minimizers_2), 
                std::ref(reference[0]), std::ref(paired_reads),
                std::ref(parameters), tasks * paired_reads.first.size() / t, (tasks + 1) * paired_reads.first.size() / t));  
      }
      thread_futures.emplace_back(thread_pool->submit(paired ? map_paired : map_as_single, 
                std::ref(ref_index), std::ref(ref_index_2), std::ref(t_minimizers), std::ref(t_minimizers_2), 
                std::ref(reference[0]), std::ref(paired_reads),
                std::ref(parameters), (t - 1) * paired_reads.first.size() / t, paired_reads.first.size()));
      
      double print_time = 0;
      for (auto& it : thread_futures) {
        it.wait();
        auto p_start = std::chrono::steady_clock::now();
        printf("%s", it.get().c_str());
        auto p_end = std::chrono::steady_clock::now();
        auto p_interval = std::chrono::duration_cast<std::chrono::duration<double>>(p_end - p_start);
        print_time += p_interval.count();
      }
      fprintf(stderr, "[countmap-map] print time: %.2f sec\n", print_time);

      auto c_end = std::chrono::steady_clock::now();
      auto c_interval = std::chrono::duration_cast<std::chrono::duration<double>>(c_end - c_start);
      fprintf(stderr, "[countmap-map] batch time: %.2f sec\n", c_interval.count());
      
      if (status1 == false || status2 == false) {
        break;
      }
    }
  } else {
    std::string reads_file(argv[optind + 1]);
    if (!(check_extension(reads_file, fasta_formats) || check_extension(reads_file, fastq_formats))) {
      fprintf(stderr, "[countmap] error: Unsupported format. Check --help for supported file formats.\n");
      exit(1);
    }

    parser_ptr_t parser;
    if (check_extension(reads_file, fasta_formats)) {
      parser = bioparser::createParser<bioparser::FastaParser, fastaq::FastAQ>(reads_file);
    } else {
      parser = bioparser::createParser<bioparser::FastqParser, fastaq::FastAQ>(reads_file);
    }

    printf("@HD\tVN:1.6\n"
           "@SQ\tSN:%s\tLN:%lu\n"
           "@PG\tID:countmap\tPN:countmap\tCL:%s %s%s %s\n", 
           reference[0]->name.c_str(), reference[0]->sequence.size(), 
           argv[0], cl_flags.c_str(), 
           argv[optind], argv[optind + 1]);

    uint32_t batch_num = 1;
    while(true) {
      fprintf(stderr, "[countmap] batch number %u\n", batch_num++);
      fprintf(stderr, "[countmap-load] loading reads... ");


      std::vector<std::unique_ptr<fastaq::FastAQ>> reads;
      bool status = parser->parse(reads, batch_size);

      fprintf(stderr, "\r[countmap-load] loaded reads        \n");

      fastaq::FastAQ::print_statistics(reads, reads_file);

      auto c_start = std::chrono::steady_clock::now();
      std::vector<std::future<std::string>> thread_futures;
      for (unsigned int tasks = 0; tasks < t - 1; ++tasks) {
        thread_futures.emplace_back(thread_pool->submit(map_single, std::ref(ref_index), std::ref(t_minimizers),
                std::ref(reference[0]), std::ref(reads),
                std::ref(parameters), tasks * reads.size() / t, (tasks + 1) * reads.size() / t));  
      }
      thread_futures.emplace_back(thread_pool->submit(map_single, std::ref(ref_index), std::ref(t_minimizers),
                std::ref(reference[0]), std::ref(reads),
                std::ref(parameters), (t - 1) * reads.size() / t, reads.size()));
      
      for (auto& it : thread_futures) {
        it.wait();
        printf("%s", it.get().c_str());
      }

      auto c_end = std::chrono::steady_clock::now();
      auto c_interval = std::chrono::duration_cast<std::chrono::duration<double>>(c_end - c_start);
      fprintf(stderr, "[countmap-map] batch time: %.2f sec\n", c_interval.count());

      if (status == false) {
        break;
      }
    }
  }

  auto m_end = std::chrono::steady_clock::now();
  auto m_interval = std::chrono::duration_cast<std::chrono::duration<double>>(m_end - m_start);
  fprintf(stderr, "[countmap-map] mapping time: %.2f sec\n", m_interval.count());

  // std::cerr << "Clip        " << clip_time << std::endl;
  // std::cerr << "            " << clip_time / m_interval.count() * 100 << "%" << std::endl;
  // std::cerr << "Clip2       " << clip_time_2 << std::endl;
  // std::cerr << "            " << clip_time_2 / m_interval.count() * 100 << "%" << std::endl;
  // std::cerr << "Minimizer   " << minimizer_time << std::endl;
  // std::cerr << "            " << minimizer_time / m_interval.count() * 100 << "%" << std::endl;
  // std::cerr << "Minimizer2  " << minimizer_time_2 << std::endl;
  // std::cerr << "            " << minimizer_time_2 / m_interval.count() * 100 << "%" << std::endl;
  // std::cerr << "Hits        " << hits_time << std::endl;
  // std::cerr << "            " << hits_time / m_interval.count() * 100 << "%" << std::endl;
  // std::cerr << "Hits2       " << hits_time_2 << std::endl;
  // std::cerr << "            " << hits_time_2 / m_interval.count() * 100 << "%" << std::endl;
  // std::cerr << "Radix       " << radix_time << std::endl;
  // std::cerr << "            " << radix_time / m_interval.count() * 100 << "%" << std::endl;
  // std::cerr << "Radix2      " << radix_time_2 << std::endl;
  // std::cerr << "            " << radix_time_2 / m_interval.count() * 100 << "%" << std::endl;
  // std::cerr << "Candidates  " << candidates_time << std::endl;
  // std::cerr << "            " << candidates_time / m_interval.count() * 100 << "%" << std::endl;
  // std::cerr << "Candidates2 " << candidates_time_2 << std::endl;
  // std::cerr << "            " << candidates_time_2 / m_interval.count() * 100 << "%" << std::endl;
  // std::cerr << "Check       " << check_time << std::endl;
  // std::cerr << "            " << check_time / m_interval.count() * 100 << "%" << std::endl;
  // std::cerr << "Check2      " << check_time_2 << std::endl;
  // std::cerr << "            " << check_time_2 / m_interval.count() * 100 << "%" << std::endl;
  // std::cerr << "Circ        " << circ_time << std::endl;
  // std::cerr << "            " << circ_time / m_interval.count() * 100 << "%" << std::endl;
  // std::cerr << "Circ2       " << circ_time_2 << std::endl;
  // std::cerr << "            " << circ_time_2 / m_interval.count() * 100 << "%" << std::endl;
  // std::cerr << "Region      " << region_time << std::endl;
  // std::cerr << "            " << region_time / m_interval.count() * 100 << "%" << std::endl;
  // std::cerr << "Expand      " << expand_time << std::endl;
  // std::cerr << "            " << expand_time / m_interval.count() * 100 << "%" << std::endl;
  // std::cerr << "Revcom      " << revcomp_time << std::endl;
  // std::cerr << "            " << revcomp_time / m_interval.count() * 100 << "%" << std::endl;
  // std::cerr << "  KSW2        " << ksw2_time << std::endl;
  // std::cerr << "              " << ksw2_time / m_interval.count() * 100 << "%" << std::endl;
  // std::cerr << "Process     " << process_time << std::endl;
  // std::cerr << "            " << process_time / m_interval.count() * 100 << "%" << std::endl;
  // std::cerr << "Process2    " << process_time_2 << std::endl;
  // std::cerr << "            " << process_time_2 / m_interval.count() * 100 << "%" << std::endl;
  // std::cerr << "SAM         " << sam_time << std::endl;
  // std::cerr << "            " << sam_time / m_interval.count() * 100 << "%" << std::endl;
  // std::cerr << "Unmapped    " << unmapped_time << std::endl;
  // std::cerr << "            " << unmapped_time / m_interval.count() * 100 << "%" << std::endl;

  // double first = clip_time + minimizer_time + hits_time + radix_time + candidates_time + check_time + circ_time + process_time;
  // double second = clip_time_2 + minimizer_time_2 + hits_time_2 + radix_time_2 + candidates_time_2 + check_time_2 + circ_time_2 + process_time_2;
  // std::cerr << "First       " << (first / m_interval.count() * 100) << "%" << std::endl;
  // std::cerr << "Second      " << (second / m_interval.count() * 100) << "%" << std::endl;

  // double prep_time = clip_time + minimizer_time + hits_time + radix_time + candidates_time + check_time + region_time + expand_time + revcomp_time;
  // std::cerr << "Preparation " << prep_time << std::endl;
  // std::cerr << "            " << prep_time / m_interval.count() * 100 << "%" << std::endl;
  // std::cerr << "Alignment   " << ksw2_time << std::endl;
  // std::cerr << "            " << clip_time / m_interval.count() * 100 << "%" << std::endl;


  return 0;
}