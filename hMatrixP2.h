#ifndef H_MATRIX_P2
#define H_MATRIX_P2

#include <vector>
#include <complex>

class HMatrixP2{
public:
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class IndexPartitioning{
        std::vector<unsigned int> data;

    public:
        IndexPartitioning(std::vector<unsigned int> _data);
        IndexPartitioning(unsigned int size);

        unsigned int total() const;
        unsigned int levels() const;
        std::vector<unsigned int> at_level(unsigned int level) const;
    };
    //\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//
private:
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class Kernel{
        std::vector<unsigned int> left_indices;
        std::vector<unsigned int> right_indices;
        std::vector<std::complex<double> > data;
    public:
        Kernel(std::vector<unsigned int> _left_indices, std::vector<unsigned int> _right_indices, const std::complex<double>* matrix);

        void apply(std::complex<double>* left, const std::complex<double>* right) const;
        long applyCount() const;
    };
    //\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class R1Element{
        std::vector<unsigned int> left_indices;
        std::vector<unsigned int> right_indices;
        std::vector<std::complex<double> > left_data;
        std::vector<std::complex<double> > right_data;

        R1Element(std::vector<unsigned int> _left_indices, std::vector<unsigned int> _right_indices, std::complex<double>* _left_data, std::complex<double>* _right_data);
    public:
        static std::vector<R1Element*> full_decomposition_at_level(const std::complex<double>* matrix, std::vector<unsigned int> left_indices, std::vector<unsigned int> right_indices);

        void apply(std::complex<double>* left, const std::complex<double>* right) const;
        long applyCount() const;
    };
    //\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//

    Kernel* kernel;
    std::vector<R1Element*> elements;

public:
    HMatrixP2(const std::complex<double>* matrix, IndexPartitioning p_left, IndexPartitioning p_right);

    void apply(std::complex<double>* left, const std::complex<double>* right) const;
    long applyCount() const;
};

#endif //H_MATRIX_P2
