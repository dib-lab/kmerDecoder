#include <string>
#include <iostream>
#include <vector>
#include <set>
#include <list>

/* 
--------------------------------------------------------
                        InputModule:Parent
--------------------------------------------------------
*/

class KDecoder{

public:
    std::list<std::string> kmers;

public:
    virtual std::list<std::string>* getKmers(std::string &seq) = 0;
};


/* 
--------------------------------------------------------
                        Default Kmers
--------------------------------------------------------
*/


class Kmers : public KDecoder{

private:
    int kSize;

public:
    Kmers(int kSize){
        this->kSize = kSize;
    }
    std::list<std::string>* getKmers(std::string &seq);
};


/* 
--------------------------------------------------------
                        Skipmers
--------------------------------------------------------
*/

class Skipmers : public KDecoder
{
private:
  int m, n, k;
  int S;

public:
  Skipmers() {}
  Skipmers(uint8_t m, uint8_t n, uint8_t k)
  {
    if (n < 1 or n < m or k < m or k > 31 or k % m != 0)
    {
      std::cout << "Error: invalid skip-mer shape! m= " << m << " n=" << n << " k= " << k << std::endl
                << "Conditions: 0 < m <= n, k <= 31 , k must multiple of m." << std::endl;

      exit(1);
    }

    this->m = m;
    this->n = n;
    this->k = k;
    this->S = k + ((k - 1) / m) * (n - m);
  }
  std::list<std::string>* getKmers(std::string &x);
  virtual ~Skipmers() {}
};


/* 
--------------------------------------------------------
                        Minimizers
--------------------------------------------------------
*/

typedef struct mkmh_minimizer
{
  uint64_t pos;
  uint32_t length;
  std::string seq;
  bool operator<(const mkmh_minimizer &rhs) const { return seq < rhs.seq; };
} mkmh_minimizer;


class Minimzers : public KDecoder
{
private:
  int k, w;
  struct mkmh_kmer_list_t
  {
    char **kmers;
    int length;
    int k;

    mkmh_kmer_list_t(){

    };

    mkmh_kmer_list_t(int length, int k)
    {
      length = length;
      k = k;
      kmers = new char *[length];
    };

    ~mkmh_kmer_list_t()
    {
      for (int i = 0; i < length; ++i)
      {
        delete[] kmers[i];
      }
      delete[] kmers;
    };
  };

protected:
  std::vector<mkmh_minimizer> kmer_tuples(std::string seq, int k);
  mkmh_kmer_list_t kmerize(char *seq, int seq_len, int k);
  std::vector<std::string> kmerize(std::string seq, int k);
  void kmerize(char *seq, const int &seq_len, const int &k, char **kmers, int &kmer_num);
  template <typename T>
  std::vector<T> v_set(std::vector<T> kmers);

public:
  Minimzers() {}
  Minimzers(int k, int w)
  {
    this->k = k;
    this->w = w;
  }

  std::vector<mkmh_minimizer> getMinimizers(std::string &seq);
  std::list<std::string>* getKmers(std::string &x);

  virtual ~Minimzers(){};
};