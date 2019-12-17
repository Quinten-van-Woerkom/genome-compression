#include <vector>
#include "dna.h"

struct matchlocation {
    std::size_t size;
    std::size_t index1, index2;
    matchlocation(int size, int index1, int index2) : size(size), index1(index1), index2(index2) {}
};

class jumps_SW
{
    protected:
    int n, m, F, L, E;

    public:

    jumps_SW(const std::vector<dna>& seq1,const std::vector<dna>& seq2, int offset, int L, int E) : L(L), E(E) {
        n = seq1.size();
        m = seq2.size();
        F = 1 + offset;
    }

    std::vector<matchlocation> match(const std::vector<dna>& seq1,const std::vector<dna>& seq2) {
        std::vector<matchlocation> matchloc;
        if (F<=n) {

            int errorloc=0;
            int cord0 = n-F;

            std::size_t I;
            if (F<m) { I = m; }
            else { I = F; }

            std::size_t i=L-1;

            while (i<I) {
                auto j = 0;
                int error = 0;

                while (j<L && error<=E) {
                    if (seq1[cord0+i-j]!=seq2[i-j]) {
                        error += 1;
                        errorloc = i-j;
                        break;
                    }
                    j++;
                }
                j++;


                while (j<L && error <= E) {
                    if (seq1[cord0+i-j] != seq2[i-j]) {
                        error++;
                    }
                    j++;
                }

                if (error>E) {
                    i = errorloc + L;
                }
                else { 
                    std::size_t index1 = cord0 -i+j-1;
                    std::size_t index2 = i-j+1;

                    while (error<=E && i<I-1) {
                        i ++;
                        if (seq1[cord0+1] != seq2[i]) {
                            error++;
                        }
                    }
                    std::size_t size = i-index2;
                    i += L;
                    matchloc.push_back(matchlocation(index1, index2, size));
                 }
            }
        }

        else if (F<m+n) {
            int cord1 = F-n;
            int errorloc=0;
            int I;
            if (F<m) { I = m; }
            else { I = F; }

            auto i = L-1;

            while (i<I) {
                int j = 0;
                int error = 0;

                while (j<L && error <= E) {
                    if (seq1[i-j] != seq2[cord1+i-j]) {
                        error++;
                        errorloc = i-j;
                        break;
                    }
                    j++;
                }
                j++;

                while (j<L && error <= E) {
                    if (seq1[i-j] != seq2[cord1 + i-j]) {
                        error++;
                    }
                    j++;
                }

                if (error>E) {
                    i = errorloc+L;
                }

                else {
                    std::size_t index1 = i-j+1;
                    std::size_t index2 = cord1 + i-j+1;

                    while (error<=E && i<I-1) {
                        i++;
                        if (seq1[i]!=seq2[cord1+i]) {error++;}
                        }
                    std::size_t size = i-index1;
                    matchloc.push_back(matchlocation(index1, index2, size));
                    }
                }
            }
        return matchloc;
    }
};