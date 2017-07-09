#include "hMatrixP2.h"

#include <iostream>
#include <vector>
#include <complex>

#include "lapacke.h"

static void svDecompose(std::complex<double>* A, unsigned int rows, unsigned int cols, double* s, std::complex<double>* U, std::complex<double>* Vadj){
        LAPACKE_zgesdd(
            /* order */ LAPACK_COL_MAJOR,
            /* jobz  */ 'A',
            /* m     */ rows,
            /* n     */ cols,
            /* A     */ A,
            /* LDA   */ rows,
            /* s     */ s,
            /* U     */ U,
            /* LDU   */ rows,
            /* VT    */ Vadj,
            /* LDVT  */ cols
    );
}

HMatrixP2::IndexPartitioning::IndexPartitioning(std::vector<unsigned int> _data): data(_data){}

HMatrixP2::IndexPartitioning::IndexPartitioning(unsigned int size){
    std::vector<unsigned int> current;
    std::vector<unsigned int> next;

    current.push_back(size);

    while(true){
        bool finish=false;
        for(unsigned int i=0; i<current.size(); i++){
            if(current[i]/2<2){
                finish=true;
                break;
            }

            next.push_back(current[i]/2);
            next.push_back(current[i]-current[i]/2);
        }

        if(finish) break;

        current=next;
        next.clear();
    }

    data=current;
}

unsigned int HMatrixP2::IndexPartitioning::total() const{
    unsigned int res=0;
    for(unsigned int i=0; i<data.size(); i++) res+=data[i];
    return res;
}

unsigned int HMatrixP2::IndexPartitioning::levels() const{
    unsigned int levels=0;
    for(;std::pow(2,levels)<data.size();levels++);
    return levels;
}

std::vector<unsigned int> HMatrixP2::IndexPartitioning::at_level(unsigned int level) const{
    unsigned int width = std::pow(2,level);
    std::vector<unsigned int> res;
    for(unsigned int i=0; i<data.size(); i++){
        if(i%width==0) res.push_back(0);
        res.back()+=data[i];
    }

    return res;
}

HMatrixP2::Kernel::Kernel(std::vector<unsigned int> _left_indices, std::vector<unsigned int> _right_indices, const std::complex<double>* matrix):
    left_indices(_left_indices),
    right_indices(_right_indices){

    unsigned int block_counter=0;
    unsigned int left_block=0;
    unsigned int right_block=0;
    unsigned int left=0;
    unsigned int right=0;

    unsigned int left_size=0;
    for(unsigned int i=0; i<left_indices.size(); i++) left_size+=left_indices[i];

    for(unsigned int i=0;;i++){
        data.push_back(matrix[left+right*left_size]);

        if(left_block==left_indices[block_counter]-1){
            if(right_block==right_indices[block_counter]-1){
                block_counter++;
                left_block=0;
                right_block=0;
                
                left++;
                right++;

                if(block_counter>=left_indices.size() or block_counter>=right_indices.size()) break;
            }else{
                left_block=0;
                right_block++;
                
                left-=left_indices[block_counter]-1;
                right++;
            }
        }else{
            left_block++;
            left++;
        }
    }
}

void HMatrixP2::Kernel::apply(std::complex<double>* left_vec, const std::complex<double>* right_vec) const{

    unsigned int block_counter=0;
    unsigned int left_block=0;
    unsigned int right_block=0;
    unsigned int left=0;
    unsigned int right=0;


    for(unsigned int i=0;i<data.size();i++){
        left_vec[left]+=data[i]*right_vec[right];

        if(left_block==left_indices[block_counter]-1){
            if(right_block==right_indices[block_counter]-1){
                block_counter++;
                left_block=0;
                right_block=0;
                
                left++;
                right++;
            }else{
                left_block=0;
                right_block++;
                
                left-=left_indices[block_counter]-1;
                right++;
            }
        }else{
            left_block++;
            left++;
        }
    }
}

long HMatrixP2::Kernel::applyCount() const{ return data.size(); }

HMatrixP2::R1Element::R1Element(std::vector<unsigned int> _left_indices, std::vector<unsigned int> _right_indices, std::complex<double>* _left_data, std::complex<double>* _right_data):
    left_indices(_left_indices),
    right_indices(_right_indices){

    unsigned int left_size=0;
    unsigned int right_size=0;
    for(unsigned int i=0; i<left_indices.size(); i++) left_size+=left_indices[i];
    for(unsigned int i=0; i<right_indices.size(); i++) right_size+=right_indices[i];

    for(unsigned int i=0; i<left_size; i++) left_data.push_back(_left_data[i]);
    for(unsigned int i=0; i<right_size; i++) right_data.push_back(_right_data[i]);
}

std::vector<HMatrixP2::R1Element*> HMatrixP2::R1Element::full_decomposition_at_level(const std::complex<double>* matrix, std::vector<unsigned int> left_indices, std::vector<unsigned int> right_indices){
    std::vector<std::vector<std::complex<double> > > left_data;
    std::vector<std::vector<std::complex<double> > > right_data;

    for(unsigned int i=0; i<left_indices.size(); i++){
        unsigned int j=i+1;
        if(i%2==1) j=i-1;

        unsigned int left_start=0;
        unsigned int right_start=0;
        unsigned int rows = left_indices[i];
        unsigned int cols = right_indices[i];

        for(unsigned int k=0; k<i; k++) left_start+=left_indices[k];
        for(unsigned int k=0; k<j; k++) right_start+=right_indices[k];


        std::vector<std::complex<double> > submatrix;
        unsigned int left_size=0;
        for(unsigned int k=0; k<left_indices.size(); k++) left_size+=left_indices[k];

        for(unsigned int right=right_start; right<right_start+cols; right++){
            for(unsigned int left=left_start; left<left_start+rows; left++){
                submatrix.push_back(matrix[left + right*left_size]);
            }
        }

        std::vector<std::complex<double> > U(rows*rows);
        std::vector<std::complex<double> > Vadj(cols*cols);
        std::vector<double> s(std::min(rows, cols));

        svDecompose(submatrix.data(), rows, cols, s.data(), U.data(), Vadj.data());

        // Store U*sigma in left_data and Vadj in right_data
        for(unsigned int k=0; k<1/*s.size()*/; k++){
            if(left_data.size()<=k) left_data.push_back(std::vector<std::complex<double> >());
            if(right_data.size()<=k) right_data.push_back(std::vector<std::complex<double> >());
            for(unsigned int l=0; l<rows; l++){
                left_data[k].push_back(U[k*rows+l]*s[k]);
            }
            for(unsigned int l=0; l<cols; l++){
                right_data[k].push_back(Vadj[l*rows+k]);
            }
        }
    }

    std::vector<R1Element*> res;

    for(unsigned int k=0; k<left_data.size(); k++){
        res.push_back(new R1Element(left_indices, right_indices, left_data[k].data(), right_data[k].data()));
    }

    return res;
}

void HMatrixP2::R1Element::apply(std::complex<double>* left, const std::complex<double>* right) const{
    unsigned int right_counter=0;
    unsigned int left_counter=0;
    unsigned int left_vec_counter=0;

    for(unsigned int i=0; i<left_indices.size(); i++){
        unsigned int j;
        if(i%2==0) j=i+1;
        else j=i-1;
       
        unsigned int right_vec_counter=0;
        for(unsigned int j_=0; j_<j; j_++) right_vec_counter+=right_indices[j_];

        std::complex<double> prod=0.;
        for(unsigned int k=0; k<right_indices[j]; k++) prod+=right_data[right_counter++]*right[right_vec_counter++];
        for(unsigned int k=0; k<left_indices[i]; k++) left[left_vec_counter++]+=prod*left_data[left_counter++];
    }
}

long HMatrixP2::R1Element::applyCount() const{ return left_data.size()+right_data.size(); }


HMatrixP2::HMatrixP2(const std::complex<double>* matrix, IndexPartitioning p_left, IndexPartitioning p_right){
    kernel = new Kernel(p_left.at_level(0),p_right.at_level(0),matrix);
    for(unsigned int i=0; i<p_left.levels(); i++){
        std::vector<R1Element*> dat = R1Element::full_decomposition_at_level(matrix, p_left.at_level(i), p_right.at_level(i));
        elements.insert(elements.end(), dat.begin(), dat.end());
    }
}

void HMatrixP2::apply(std::complex<double>* left, const std::complex<double>* right) const{
    kernel->apply(left, right);
    for(unsigned int i=0; i<elements.size(); i++) elements[i]->apply(left, right);
}

long HMatrixP2::applyCount() const{
    long res=kernel->applyCount();
    for(unsigned int i=0; i<elements.size(); i++) res+=elements[i]->applyCount();
    return res;
}
